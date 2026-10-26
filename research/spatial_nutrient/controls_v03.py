"""Registered isolated lifecycle controls for the prospective v0.3 study.

These are Own-model qualification fixtures, not biological validation. They use
the production experiment scheduler, uptake, stochastic fission, death and mass
ledgers. All descendants are observed to the full 600 s horizon. See
PROTOCOL_V03.md, isolated controls; DOCUMENTATION.md, lifecycle qualification.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
import copy
import gzip
import hashlib
import io
import json
import math
import multiprocessing
import os
from pathlib import Path
import random
import traceback

from .inputs import seed_streams, scoped_python_random, _tuples
from .runner import ROOT, atomic_json, canonical, code_identity, digest, runtime_identity

PLAN_PATH = ROOT / 'research/spatial_nutrient/study_plan_v03.json'
BASE_PATH = ROOT / 'research/spatial_nutrient/protocol_v02.json'
CONDITIONS = ('no_feed', 'abundant_local_feed')
ROLES = ('focal', 'resident')
LIFESPANS = (.4, 1.0, 1.6)
SEEDS = tuple(range(6000, 6020))


def fixture_specs():
    """Return the entire prespecified denominator in deterministic order."""
    return [dict(fixture_id=f'{condition}__{role}__L{str(lifespan).replace(".", "p")}__seed{seed}',
                 condition=condition, role=role, lifespan_multiplier=lifespan, seed=seed)
            for seed in SEEDS for condition in CONDITIONS for role in ROLES
            for lifespan in LIFESPANS]


def _validate_plan(plan, base):
    c = plan['controls']
    exact = {'seeds': list(SEEDS), 'lifespan_multipliers': list(LIFESPANS),
             'size': 1, 'initial_reserve_pg': 2.5, 'horizon_s': 600, 'dt_s': .1,
             'movement': 0, 'maintenance_aging': False, 'daughter_size_policy': 'fixed',
             'environmental_and_age_mortality': False}
    for name, value in exact.items():
        if c[name] != value:
            raise ValueError(f'Unsupported prespecified control setting: {name}')
    if digest(base) != plan['base_configuration']['canonical_sha256']:
        raise ValueError('Historical base protocol digest mismatch')
    if plan['study_version'] != '0.3':
        raise ValueError('This harness implements only the v0.3 controls')


def _resolved_config(plan, base, spec):
    c = plan['controls']
    species = copy.deepcopy(next(s for s in base['species'] if s['role'] == spec['role']))
    # Neutral pH and oxygen-indifference give envf=1, making the production
    # environmental death probability exactly zero; base_aging=0 removes age
    # mortality without patching the lifecycle code or its RNG sequence.
    species.update(init=1, speed=c['movement'], aging=0, so=0, ph=7)
    species.pop('role')
    environment = copy.deepcopy(base['environment_tick_start_clamp'])
    environment.update(ph=7, drop_frequency=0, phage_inflow=0)
    return {
        'species': species, 'environment': environment,
        'geometry': copy.deepcopy(base['geometry']),
        'capacity_cells': base['capacity_cells'], 'steps': 6000,
        'dt_s': c['dt_s'], 'horizon_s': c['horizon_s'],
        'founder': {'size': c['size'], 'initial_reserve_pg': c['initial_reserve_pg'],
                    'age': 0, 'lifespan_multiplier': spec['lifespan_multiplier']},
        'daughter_lifespan': 'production constructor draws; no age-dependent maintenance or mortality',
        'abundant_node_pg': 100, 'substrate_id': 0,
        'abundance_rule': 'one node at each live cell; top up to 100 pg before every tick; retain all nodes',
        'ledger': {'maintenance_aging': False, 'daughter_size_policy': 'fixed',
                   'recycling': base['lifecycle']['recycling'],
                   'environment_clamp': environment,
                   'washout_rate': base['neutralized_radial_parameters']['LUMEN_FOOD_WASHOUT_RATE_PER_SEC'],
                   'motility_factor': 1, 'source_angle': 0,
                   'max_food_nodes': None, 'mass_tolerance': base['execution']['mass_tolerance']}}


def _make_species(config, role):
    from bacteria.core_sim import BacteriaSpecies
    sp = BacteriaSpecies(**config['species'], color=(230, 80, 80) if role == 'focal' else (80, 130, 230))
    sp.population_role = role
    sp.maintenance_aging = False
    sp.daughter_size_policy = 'fixed'
    return sp


def _make_input(plan, base, spec):
    from bacteria.core_sim import Bacterium
    config = _resolved_config(plan, base, spec)
    streams = seed_streams(spec['seed'])
    with scoped_python_random(streams['founders']['python_state']):
        sp = _make_species(config, spec['role'])
        cell = Bacterium(sp, *config['geometry']['center'], size=config['founder']['size'],
                         initial_food=config['founder']['initial_reserve_pg'])
        cell.agent_id = 0
        cell.age = 0
        cell.lifespan_variance = spec['lifespan_multiplier']
        founder = {k: copy.deepcopy(v) for k, v in vars(cell).items() if k != 'species'}
    return json.loads(canonical({'fixture': spec, 'resolved_config': config,
                                 'streams': streams, 'founder_state': founder}))


def build_control(inputs, event_sink=None):
    """Restore a saved single-founder fixture; caller owns process RNG scope."""
    from bacteria.core_sim import Bacterium, Simulation
    config, spec = inputs['resolved_config'], inputs['fixture']
    sp = _make_species(config, spec['role'])
    sim = Simulation([sp], copy.deepcopy(config['environment']), config['capacity_cells'],
                     steps=config['steps'], physics_backend='python', ecm_backend='python')
    geom = config['geometry']
    sim._ensure_petri_geometry(geom['world_width_um'], geom['world_height_um'])
    sim.petri_center = tuple(geom['center'])
    sim.petri_radius = math.hypot(geom['world_width_um'], geom['world_height_um']) / 2 + geom['padding_um']
    cell = Bacterium(sp, *sim.petri_center, size=config['founder']['size'],
                     initial_food=config['founder']['initial_reserve_pg'])
    for key, value in inputs['founder_state'].items():
        setattr(cell, key, copy.deepcopy(value))
    sp.members = [cell]
    sim._populated = True
    sim.enable_experiment(config['ledger'], inputs['streams']['scheduler']['seed'], event_sink)
    random.setstate(_tuples(inputs['streams']['ecology']['python_state']))
    sim.running = True
    return sim


def replenish_abundance(sim, nodes_by_agent, amount_pg=100):
    """Maintain the declared per-cell node without accumulating duplicate boluses.

    Existing nodes remain in the external pool. Every top-up is a positive
    external mass transfer and is recorded before begin_tick, exactly as a new
    scheduled bolus is. No nutrient is removed or reset without accounting.
    """
    live_nodes = {id(node) for node in sim.food_nodes}
    for cell in sorted((b for sp in sim.species for b in sp.members), key=lambda b: b.agent_id):
        node = nodes_by_agent.get(cell.agent_id)
        if node is None or id(node) not in live_nodes:
            sim.add_external_food(cell.x, cell.y, amount_pg, substrate_id=0,
                                  event_id=f'control:{sim.step_idx}:{cell.agent_id}')
            node = nodes_by_agent[cell.agent_id] = sim.food_nodes[-1]
        else:
            if node.x != cell.x or node.y != cell.y or node.substrate_id != 0:
                raise ValueError('Immobilized abundance fixture moved relative to its food node')
            addition = max(0.0, amount_pg - node.amount)
            node.amount += addition
            if addition:
                sim.experiment.record('feed', amount_pg=addition, x=node.x, y=node.y,
                                      substrate_id=0,
                                      event_id=f'control:{sim.step_idx}:{cell.agent_id}',
                                      replenished_agent_id=cell.agent_id)
        if node.amount != amount_pg:
            raise ValueError('Abundance node is not exactly the declared 100 pg before tick')


@contextmanager
def _gzip_jsonlines(path):
    # Fixed header makes identical retained logs byte-identical on replay.
    with Path(path).open('wb') as raw:
        try:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, compresslevel=1, mtime=0) as zipped:
                with io.TextIOWrapper(zipped, encoding='utf-8') as text:
                    yield lambda value: text.write(canonical(value) + '\n')
        finally:
            # Preserve complete gzip footers for failed attempts as well.
            raw.flush()
            os.fsync(raw.fileno())


def _file_metadata(path):
    path = Path(path)
    with path.open('rb') as handle:
        sha = hashlib.file_digest(handle, 'sha256').hexdigest()
    return {'sha256': sha, 'bytes': path.stat().st_size}


def register_controls(output, *, plan_path=PLAN_PATH, base_path=BASE_PATH):
    """Register all 240 exact inputs before any control is executed."""
    output = Path(output)
    plan, base = json.loads(Path(plan_path).read_text()), json.loads(Path(base_path).read_text())
    _validate_plan(plan, base)
    output.mkdir(parents=True, exist_ok=False)
    manifest = {'schema_version': 1, 'cohort': 'isolated_controls_v03',
                'study_plan': plan, 'study_plan_sha256': digest(plan),
                'base_protocol': base, 'base_protocol_sha256': digest(base),
                'code': code_identity(), 'runtime': runtime_identity(), 'fixtures': []}
    for spec in fixture_specs():
        inputs = _make_input(plan, base, spec)
        atomic_json(output / 'inputs' / (spec['fixture_id'] + '.json'), inputs)
        manifest['fixtures'].append({**spec, 'input_sha256': digest(inputs)})
    manifest['manifest_sha256'] = digest(manifest)
    atomic_json(output / 'manifest.json', manifest)
    return manifest


def _read_manifest(output, *, check_runtime=True):
    manifest = json.loads((Path(output) / 'manifest.json').read_text())
    saved = manifest.pop('manifest_sha256')
    if digest(manifest) != saved:
        raise ValueError('Controls manifest hash mismatch')
    manifest['manifest_sha256'] = saved
    _validate_plan(manifest['study_plan'], manifest['base_protocol'])
    if digest(manifest['study_plan']) != manifest['study_plan_sha256'] or digest(manifest['base_protocol']) != manifest['base_protocol_sha256']:
        raise ValueError('Controls configuration hash mismatch')
    specs = [{k: row[k] for k in fixture_specs()[0]} for row in manifest['fixtures']]
    if specs != fixture_specs():
        raise ValueError('Controls denominator must be the entire240 prespecified fixtures')
    if check_runtime and (manifest['code']['sha256'] != code_identity()['sha256'] or manifest['runtime'] != runtime_identity()):
        raise ValueError('Controls executable/runtime changed after registration')
    return manifest


def _read_input(output, manifest, row):
    inputs = json.loads((Path(output) / 'inputs' / (row['fixture_id'] + '.json')).read_text())
    if digest(inputs) != row['input_sha256']:
        raise ValueError('Controls input hash mismatch')
    expected = {key: row[key] for key in fixture_specs()[0]}
    if inputs['fixture'] != expected or inputs['resolved_config'] != _resolved_config(manifest['study_plan'], manifest['base_protocol'], expected):
        raise ValueError('Controls input differs from prespecified fixture')
    if inputs != _make_input(manifest['study_plan'], manifest['base_protocol'], expected):
        raise ValueError('Controls founder or random streams differ from prespecified fixture')
    return inputs


def _outcome_passes(result):
    if result['status'] != 'completed' or result['completed_ticks'] != 6000 or not result['ledger_closed']:
        return False
    if result['deaths_by_cause'].get('capacity', 0) or result['deaths_by_cause'].get('environment_or_age', 0):
        return False
    if result['fixture']['condition'] == 'no_feed':
        return (result['births'] == 0 and result['final_count'] == 0 and
                result['deaths_by_cause'] == {'starvation': 1} and
                result['first_starvation_s'] is not None and result['first_starvation_s'] < 600 and
                result['cumulative']['feed_pg'] == 0)
    return (result['births'] > 0 and result['deaths'] == 0 and
            result['first_birth_s'] is not None and result['first_birth_s'] < 600 and
            result['abundance_ticks'] == 6000 and result['abundance_checks'] >= 6000 and
            result['abundance_failures'] == 0 and result['final_count'] == result['births'] + 1)


def _run_fixture(output, fixture_id):
    output = Path(output)
    manifest = _read_manifest(output)
    row = next(f for f in manifest['fixtures'] if f['fixture_id'] == fixture_id)
    inputs = _read_input(output, manifest, row)
    target = output / 'fixtures' / fixture_id
    target.mkdir(parents=True, exist_ok=False)  # Attempts are never overwritten.
    result = {'fixture': inputs['fixture'], 'input_sha256': row['input_sha256'],
              'manifest_sha256': manifest['manifest_sha256'], 'status': 'running',
              'completed_ticks': 0, 'ledger_closed': False, 'births': 0, 'deaths': 0,
              'first_birth_s': None, 'first_starvation_s': None, 'deaths_by_cause': {},
              'abundance_ticks': 0, 'abundance_checks': 0, 'abundance_failures': 0,
              'max_abs_mass_residual_pg': 0, 'max_abs_reserve_residual_pg': 0}
    atomic_json(target / 'started.json', {'fixture_id': fixture_id, 'input_sha256': row['input_sha256'],
                                         'manifest_sha256': manifest['manifest_sha256']})
    sim = None
    try:
        config = inputs['resolved_config']
        with _gzip_jsonlines(target / 'events.jsonl.gz') as write_event, _gzip_jsonlines(target / 'trajectory.jsonl.gz') as write_tick:
            def observe(event):
                write_event(event)
                time = (event['tick'] + 1) * config['dt_s']
                if event['kind'] == 'birth' and result['first_birth_s'] is None:
                    result['first_birth_s'] = time
                if event['kind'] == 'death':
                    cause = event['cause']
                    result['deaths_by_cause'][cause] = result['deaths_by_cause'].get(cause, 0) + 1
                    if cause == 'starvation' and result['first_starvation_s'] is None:
                        result['first_starvation_s'] = time
            with scoped_python_random(inputs['streams']['ecology']['python_state']):
                sim = build_control(inputs, observe)
                nodes = {}
                for _ in range(config['steps']):
                    if inputs['fixture']['condition'] == 'abundant_local_feed':
                        count = len(sim.species[0].members)
                        try:
                            replenish_abundance(sim, nodes, config['abundant_node_pg'])
                        except Exception:
                            result['abundance_failures'] += 1
                            raise
                        result['abundance_checks'] += count
                        result['abundance_ticks'] += 1
                    sim.update(config['geometry']['world_width_um'], config['geometry']['world_height_um'], config['dt_s'])
                    tick = sim.experiment.last_tick
                    write_tick(tick)
                    result['completed_ticks'] += 1
                    result['max_abs_mass_residual_pg'] = max(result['max_abs_mass_residual_pg'], abs(tick['mass_residual_pg']))
                    result['max_abs_reserve_residual_pg'] = max(result['max_abs_reserve_residual_pg'], *(abs(p['reserve_residual_pg']) for p in tick['populations'].values()))
                atomic_json(target / 'rng_final.json', {'ecology': random.getstate(), 'scheduler': sim.experiment.scheduler.getstate()})
        result.update(status='completed', ledger_closed=True)
    except Exception as exc:
        result.update(status='failed', error=f'{type(exc).__name__}: {exc}', traceback=traceback.format_exc())
    if sim is not None:
        result.update(cumulative=sim.experiment.cumulative, final_count=len(sim.species[0].members),
                      births=sim.experiment.cumulative['births'], deaths=sim.experiment.cumulative['deaths'],
                      temperature_factor=sim.species[0].temperature_factor(config['environment']['temperature']),
                      final_external_mass_pg=math.fsum(n.amount for n in sim.food_nodes),
                      final_reserve_pg=math.fsum(b.food_level for b in sim.species[0].members))
    result['passed'] = _outcome_passes(result)
    result['artifacts'] = {p.name: _file_metadata(p) for p in sorted(target.iterdir()) if p.is_file()}
    result['result_sha256'] = digest(result)
    atomic_json(target / 'result.json', result)
    return result


def _read_result(output, manifest, row):
    target = Path(output) / 'fixtures' / row['fixture_id']
    result = json.loads((target / 'result.json').read_text())
    saved = result.pop('result_sha256')
    if digest(result) != saved:
        raise ValueError('Control result hash mismatch')
    result['result_sha256'] = saved
    if result['manifest_sha256'] != manifest['manifest_sha256'] or result['input_sha256'] != row['input_sha256'] or result['fixture'] != {key: row[key] for key in fixture_specs()[0]}:
        raise ValueError('Control result registration mismatch')
    required = {'started.json', 'events.jsonl.gz', 'trajectory.jsonl.gz'}
    if result['status'] == 'completed':
        required.add('rng_final.json')
    if set(result['artifacts']) != required:
        raise ValueError('Control source artifact set mismatch')
    for name, metadata in result['artifacts'].items():
        if _file_metadata(target / name) != metadata:
            raise ValueError(f'Control source artifact hash mismatch: {name}')
    if result['status'] == 'completed':
        inputs = _read_input(output, manifest, row)
        _verify_source_outcomes(target, inputs, result)
    if result['passed'] != _outcome_passes(result):
        raise ValueError('Control outcome gate mismatch')
    return result


def _verify_source_outcomes(target, inputs, result):
    """Reconstruct essential outcomes from events and independently close ticks.

    A summary is not evidence for its own gate: even a rehashed summary must
    agree with the retained transfer/lifecycle records and population snapshots.
    Reading one tick at a time keeps the full-horizon audit memory bounded.
    """
    config = inputs['resolved_config']
    name, dt = config['species']['name'], config['dt_s']
    tolerance = config['ledger']['mass_tolerance']
    transfer_keys = {'feed': 'feed_pg', 'recycling': 'recycled_pg', 'uptake': 'uptake_pg',
                     'washout': 'washout_pg', 'prune': 'pruned_pg', 'maintenance': 'maintenance_pg',
                     'fission': 'fission_cost_pg', 'death_reserve': 'death_reserve_pg'}
    cumulative = {key: 0.0 for key in (*transfer_keys.values(), 'reserve_gain_pg',
                                      'maintenance_demand_pg', 'unmet_maintenance_pg')}
    cumulative.update(births=0, deaths=0)
    summary = dict(completed_ticks=0, ledger_closed=True, first_birth_s=None,
                   first_starvation_s=None, deaths_by_cause={}, abundance_ticks=0,
                   abundance_checks=0, abundance_failures=0,
                   max_abs_mass_residual_pg=0, max_abs_reserve_residual_pg=0)
    live_ids, last_mass, last_reserve, elapsed = {0}, 0.0, config['founder']['initial_reserve_pg'], 0.0
    with gzip.open(target / 'events.jsonl.gz', 'rt') as event_file, gzip.open(target / 'trajectory.jsonl.gz', 'rt') as tick_file:
        events = (json.loads(line) for line in event_file)
        event = next(events, None)
        for line in tick_file:
            tick = json.loads(line)
            tick_index = summary['completed_ticks']
            if tick['tick'] != tick_index:
                raise ValueError('Controls source ticks are missing, duplicated, or out of order')
            start_ids, before = set(live_ids), dict(cumulative)
            gains, losses, feeds, births, deaths = [], [], [], 0, {}
            while event is not None and event['tick'] == tick_index:
                kind = event['kind']
                if kind in transfer_keys:
                    amount = event['amount_pg']
                    if not math.isfinite(amount) or amount < 0:
                        raise ValueError('Controls source contains invalid transfer')
                    cumulative[transfer_keys[kind]] += amount
                if kind == 'feed':
                    fields = event['event_id'].split(':')
                    if len(fields) != 3 or fields[:2] != ['control', str(tick_index)] or event['substrate_id'] != 0:
                        raise ValueError('Controls source feed is not the declared local replenishment')
                    feeds.append(int(fields[2]))
                if kind == 'maintenance':
                    demand = event['demand_pg']
                    if not math.isfinite(demand) or demand < amount or event['age_multiplier'] != 1:
                        raise ValueError('Controls source maintenance differs from the declared policy')
                    cumulative['maintenance_demand_pg'] += demand
                    cumulative['unmet_maintenance_pg'] += max(0.0, demand - amount)
                if kind == 'uptake':
                    gain = event['reserve_gain_pg']
                    cumulative['reserve_gain_pg'] += gain
                    gains.append(gain)
                if kind in ('maintenance', 'fission', 'death_reserve'):
                    losses.append(amount)
                if kind == 'birth':
                    if event['agent_id'] in live_ids or event['parent_id'] not in live_ids or event['base_size'] != 1 or event['capacity_pg'] != 100:
                        raise ValueError('Controls source birth violates lineage/fixed-size policy')
                    live_ids.add(event['agent_id'])
                    births += 1
                    cumulative['births'] += 1
                    if summary['first_birth_s'] is None:
                        summary['first_birth_s'] = (tick_index + 1) * dt
                if kind == 'death':
                    if event['agent_id'] not in live_ids:
                        raise ValueError('Controls source death has no living cell')
                    live_ids.remove(event['agent_id'])
                    cause = event['cause']
                    deaths[cause] = deaths.get(cause, 0) + 1
                    summary['deaths_by_cause'][cause] = summary['deaths_by_cause'].get(cause, 0) + 1
                    cumulative['deaths'] += 1
                    if cause == 'starvation' and summary['first_starvation_s'] is None:
                        summary['first_starvation_s'] = (tick_index + 1) * dt
                event = next(events, None)
            if event is not None and event['tick'] < tick_index + 1:
                raise ValueError('Controls source events are out of order')
            delta = {key: value - before[key] for key, value in cumulative.items()}
            if any(tick[key] != value for key, value in delta.items()):
                raise ValueError('Controls source tick transfers differ from events')
            population = tick['populations'][name]
            if (set(tick['populations']) != {name} or population['start'] != len(start_ids)
                    or population['end'] != len(live_ids) or population['births'] != births
                    or population['deaths'] != deaths):
                raise ValueError('Controls source population counts differ from events')
            snapshot = tick['exposure']['populations'][name]
            if (snapshot['agent_ids'] != sorted(live_ids)
                    or math.fsum(snapshot['reserves_pg']) != population['reserve_pg']
                    or any(not math.isfinite(value) or value < 0 or value > 100 + tolerance for value in snapshot['reserves_pg'])):
                raise ValueError('Controls source live-cell snapshot differs from population ledger')
            mass_expected = (last_mass + delta['feed_pg'] + delta['recycled_pg']
                             - delta['uptake_pg'] - delta['washout_pg'] - delta['pruned_pg'])
            mass_residual = tick['mass_end_pg'] - mass_expected
            reserve_expected = last_reserve + math.fsum(gains) - math.fsum(losses)
            reserve_residual = population['reserve_pg'] - reserve_expected
            if (tick['mass_start_pg'] != last_mass or not math.isfinite(tick['mass_end_pg'])
                    or tick['mass_residual_pg'] != mass_residual
                    or abs(mass_residual) > tolerance * max(1, last_mass + delta['feed_pg'] + delta['recycled_pg'])
                    or population['reserve_residual_pg'] != reserve_residual
                    or abs(reserve_residual) > tolerance * max(1, last_reserve + math.fsum(gains))):
                raise ValueError('Controls source mass/reserve ledger does not close')
            if inputs['fixture']['condition'] == 'abundant_local_feed':
                if len(feeds) != len(start_ids) or set(feeds) != start_ids:
                    raise ValueError('Controls source lacks per-live-cell replenishment')
                # With no deaths, one100 pg node per cell accounts for the entire
                # external pool before washout. Dead-cell nodes are retained.
                if not before['deaths'] and abs(last_mass + delta['feed_pg'] - 100 * len(start_ids)) > tolerance * max(1, 100 * len(start_ids)):
                    raise ValueError('Controls source pre-tick abundance does not match100 pg per cell')
                summary['abundance_ticks'] += 1
                summary['abundance_checks'] += len(start_ids)
            elif feeds:
                raise ValueError('Controls source no-feed fixture contains external input')
            elapsed += dt
            if tick['time_seconds'] != elapsed:
                raise ValueError('Controls source elapsed time differs from the prescribed clock')
            last_mass, last_reserve = tick['mass_end_pg'], population['reserve_pg']
            summary['completed_ticks'] += 1
            summary['max_abs_mass_residual_pg'] = max(summary['max_abs_mass_residual_pg'], abs(mass_residual))
            summary['max_abs_reserve_residual_pg'] = max(summary['max_abs_reserve_residual_pg'], abs(reserve_residual))
        if event is not None:
            raise ValueError('Controls source has events beyond the retained completed horizon')
    summary.update(cumulative=cumulative, final_count=len(live_ids), births=cumulative['births'],
                   deaths=cumulative['deaths'], final_external_mass_pg=last_mass,
                   final_reserve_pg=last_reserve)
    for key, value in summary.items():
        if result.get(key) != value:
            raise ValueError(f'Controls result summary differs from retained sources: {key}')


def _make_report(output, manifest):
    results = []
    for row in manifest['fixtures']:
        _read_input(output, manifest, row)
        results.append(_read_result(output, manifest, row))
    passed = sum(r['passed'] for r in results)
    report = {'schema_version': 1, 'cohort': 'isolated_controls_v03',
              'status': 'passed' if passed == 240 else 'failed', 'required_fixtures': 240,
              'passed_fixtures': passed, 'failed_fixtures': 240 - passed,
              'manifest_sha256': manifest['manifest_sha256'],
              'study_plan_sha256': manifest['study_plan_sha256'],
              'base_protocol_sha256': manifest['base_protocol_sha256'],
              'code': manifest['code'], 'runtime': manifest['runtime'], 'fixture_results': results}
    report['report_sha256'] = digest(report)
    return report


def run_controls(output, workers=1):
    """Execute the complete registration, retaining failures without replacement."""
    output = Path(output)
    manifest = _read_manifest(output)
    # Already completed cases may be inspected/resumed after an interruption;
    # partial attempts are retained and block this execution, never silently reset.
    pending = []
    for row in manifest['fixtures']:
        target = output / 'fixtures' / row['fixture_id']
        if (target / 'result.json').exists():
            _read_result(output, manifest, row)
        elif target.exists():
            raise ValueError(f'Incomplete retained control attempt: {target}')
        else:
            pending.append(row['fixture_id'])
    if type(workers) is not int or workers < 1:
        raise ValueError('workers must be a positive integer')
    if workers == 1:
        for fixture_id in pending:
            result = _run_fixture(output, fixture_id)
            print(f'{fixture_id}: {result["status"]}, passed={result["passed"]}', flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as pool:
            futures = [pool.submit(_run_fixture, str(output), fixture_id) for fixture_id in pending]
            for future in as_completed(futures):
                result = future.result()
                print(f'{result["fixture"]["fixture_id"]}: {result["status"]}, passed={result["passed"]}', flush=True)
    # Reject files edited during the run as well as before registration.
    _read_manifest(output)
    report = _make_report(output, manifest)
    atomic_json(output / 'report.json', report)
    return report


def verify_controls(path, *, code_sha256=None, runtime=None, plan_sha256=None):
    """Verify all retained sources, the exact denominator, and the aggregate gate."""
    output = Path(path)
    if output.is_file():
        output = output.parent
    manifest = _read_manifest(output, check_runtime=False)
    if code_sha256 is not None and manifest['code']['sha256'] != code_sha256:
        raise ValueError('Controls evidence executable mismatch')
    if runtime is not None and manifest['runtime'] != runtime:
        raise ValueError('Controls evidence runtime mismatch')
    if plan_sha256 is not None and manifest['study_plan_sha256'] != plan_sha256:
        raise ValueError('Controls evidence study-plan mismatch')
    report = _make_report(output, manifest)
    if json.loads((output / 'report.json').read_text()) != report:
        raise ValueError('Controls report differs from retained sources')
    if report['status'] != 'passed':
        raise ValueError('Controls gate failed; all 240 fixtures must pass before the grid')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('register', 'run', 'verify'))
    parser.add_argument('output', type=Path)
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    if args.command == 'register':
        value = register_controls(args.output)
        print(value['manifest_sha256'])
    elif args.command == 'run':
        value = run_controls(args.output, args.workers)
        print(canonical({k: value[k] for k in ('status', 'passed_fixtures', 'failed_fixtures', 'report_sha256')}))
    else:
        print(verify_controls(args.output)['report_sha256'])


if __name__ == '__main__':
    main()
