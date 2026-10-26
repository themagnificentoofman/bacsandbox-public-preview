"""Registered, append-only paired execution for the spatial nutrient experiment.

Every attempt is retained; retries restart the same saved inputs, never a new seed.
Only subprocess workers mutate the legacy engine's process-global controls.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import ExitStack, contextmanager
from fractions import Fraction
import hashlib
import gzip
import io
import importlib.metadata
import json
import multiprocessing
import os
from pathlib import Path
import platform
import random
import sys
import traceback

from .protocol import load_protocol, validate_protocol
from .inputs import make_pair_inputs, arm_inputs, assert_matched_arms, tick_schedule, build_simulation

ROOT = Path(__file__).resolve().parents[2]


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    with temp.open('w') as handle:
        json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp, path)


def code_identity():
    """Hash executable sources/assets, not mutable outputs or commit messages."""
    paths = sorted([*ROOT.joinpath('bacteria').rglob('*.py'),
                    *ROOT.joinpath('bacteria/models').rglob('*.json'),
                    *ROOT.joinpath('research/spatial_nutrient').glob('*.py')])
    files = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if not p.name.startswith('test_')}
    return {'sha256': digest(files), 'files': files}


def runtime_identity():
    packages = {}
    for name in ('numpy', 'pygame-ce', 'scipy'):
        packages[name] = importlib.metadata.version(name)
    return {'python': platform.python_version(), 'implementation': platform.python_implementation(),
            'platform': platform.platform(), 'packages': packages, 'backend': 'python'}


def registered_seeds(protocol, cohort):
    if cohort == 'pilot':
        return protocol['replication']['pilot_seeds']
    if cohort == 'main':
        return protocol['replication']['main_seeds']
    if cohort in ('numerical_coarse', 'numerical_fine'):
        return protocol['replication']['numerical_seeds']
    raise ValueError('diagnostic and sensitivity cohorts require explicit registered seeds')


def _verify_main_freeze(freeze, protocol, code, runtime):
    if freeze is None:
        raise ValueError('main execution requires a successful qualification freeze')
    frozen = json.loads(Path(freeze).read_text()) if isinstance(freeze, (str, Path)) else freeze
    study_v03 = any(a['change'].startswith('study_v03:') for a in protocol.get('amendments', []))
    if study_v03 and frozen.get('study_version') != '0.3':
        raise ValueError('v0.3 main execution requires a v0.3 qualification freeze')
    if frozen.get('study_version') == '0.3':
        from .qualification_v03 import verify_freeze
    else:
        from .qualification_run import verify_freeze
    return verify_freeze(frozen, protocol_sha256=digest(protocol),
                         code_sha256=code['sha256'], runtime=runtime)


def register(output, protocol, cohort, seeds=None, *, freeze=None):
    """Save the complete denominator, both arms and execution order before any run."""
    output = Path(output)
    protocol = load_protocol(protocol) if isinstance(protocol, (str, Path)) or protocol is None else protocol
    protocol = validate_protocol(protocol)
    seeds = list(registered_seeds(protocol, cohort) if seeds is None else seeds)
    if not seeds or len(set(seeds)) != len(seeds) or any(type(s) is not int or s < 0 for s in seeds):
        raise ValueError('seeds must be unique nonnegative integers')
    if cohort in ('pilot', 'main') and seeds != registered_seeds(protocol, cohort):
        raise ValueError('pilot/main denominator must equal the entire registered cohort')
    main_seeds = set(protocol['replication']['main_seeds'])
    if cohort != 'main' and main_seeds.intersection(seeds):
        raise ValueError('fresh main seeds cannot be used for qualification or sensitivity')
    code, runtime = code_identity(), runtime_identity()
    if cohort == 'main':
        freeze = _verify_main_freeze(freeze, protocol, code, runtime)
    output.mkdir(parents=True, exist_ok=False)
    manifest = {'schema_version': 1, 'cohort': cohort, 'protocol': protocol,
                'protocol_sha256': digest(protocol), 'code': code, 'runtime': runtime,
                'seeds': seeds, 'pairs': [], 'freeze': freeze}
    # A partially written registration has no final manifest and is never runnable.
    atomic_json(output / 'registration_started.json', {'seeds': seeds, 'cohort': cohort})
    for seed in seeds:
        pair = make_pair_inputs(protocol, seed)
        u, h = arm_inputs(pair, 'U'), arm_inputs(pair, 'H')
        assert_matched_arms(u, h)
        pair_dir = output / 'inputs' / str(seed)
        pair_dir.mkdir(parents=True)
        for arm, inputs in [('U', u), ('H', h)]:
            with _compressed_jsonl(pair_dir / f'{arm}.json.gz') as handle:
                handle.write(canonical(inputs) + '\n')
        # Execution order is not the ecological RNG stream or a treatment input.
        order = ['U', 'H']
        random.Random(pair['streams']['execution_order']['seed']).shuffle(order)
        manifest['pairs'].append({'seed': seed, 'order': order,
                                  'inputs': {'U': digest(u), 'H': digest(h)}})
    manifest['manifest_sha256'] = digest(manifest)
    atomic_json(output / 'manifest.json', manifest)
    return manifest


def read_manifest(output, check_runtime=True):
    manifest = json.loads((Path(output) / 'manifest.json').read_text())
    expected = manifest.pop('manifest_sha256')
    if digest(manifest) != expected:
        raise ValueError('manifest hash mismatch')
    manifest['manifest_sha256'] = expected
    if check_runtime:
        if manifest['code']['sha256'] != code_identity()['sha256']:
            raise ValueError('code changed; register a new version and rerun the cohort')
        if manifest['runtime'] != runtime_identity():
            raise ValueError('runtime changed; use the registered runtime')
    return manifest


def _json_state(value):
    if isinstance(value, tuple):
        return [_json_state(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _json_state(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_state(v) for v in value]
    if hasattr(value, 'item'):
        return value.item()
    return value


def _last_counts(sim):
    return {sp.name: len(sp.members) for sp in sim.species}


def artifact_path(directory, logical_name):
    """Resolve one retained format, rejecting ambiguous plain/compressed copies."""
    directory = Path(directory)
    paths = [directory / logical_name]
    if logical_name.endswith(('.jsonl', '.json')):
        paths.append(directory / (logical_name + '.gz'))
    found = [p for p in paths if p.exists()]
    if len(found) != 1:
        raise ValueError(f'Expected one retained artifact for {logical_name}, found {len(found)}')
    return found[0]


def open_jsonl(path):
    """Read lossless compressed evidence and historical plain JSONL alike."""
    path = Path(path)
    return gzip.open(path, 'rt', encoding='utf-8') if path.suffix == '.gz' else path.open(encoding='utf-8')


def input_path(output, seed, arm):
    return artifact_path(Path(output) / 'inputs' / str(seed), f'{arm}.json')


def read_inputs(output, seed, arm):
    with open_jsonl(input_path(output, seed, arm)) as handle:
        return json.load(handle)


@contextmanager
def _compressed_jsonl(path):
    # Fixed header metadata makes identical ecological replays byte-identical.
    # fsync follows the gzip footer, before any completion record can be written.
    with Path(path).open('wb') as raw:
        try:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=3) as compressed:
                with io.TextIOWrapper(compressed, encoding='utf-8', newline='\n') as handle:
                    yield handle
        finally:
            raw.flush()
            os.fsync(raw.fileno())


def artifact_manifest(directory):
    result = {}
    for name in ('events.jsonl', 'trajectory.jsonl', 'rng_final.json'):
        path = artifact_path(directory, name)
        hasher = hashlib.sha256()
        with path.open('rb') as handle:
            for chunk in iter(lambda: handle.read(1024*1024), b''):
                hasher.update(chunk)
        result[path.name] = {'sha256':hasher.hexdigest(), 'size_bytes':path.stat().st_size}
    return result


def validate_completed(record, attempt_dir, manifest, pair_meta, seed, arm):
    if (record.get('status') != 'completed' or record.get('validated') is not True
            or record.get('seed') != seed or record.get('arm') != arm
            or record.get('input_sha256') != pair_meta['inputs'][arm]
            or record.get('manifest_sha256') != manifest['manifest_sha256']
            or record.get('outputs') != artifact_manifest(attempt_dir)):
        raise ValueError('Completed result or retained output integrity mismatch')


def recover_interrupted(output, seed, arm):
    """Release only a dead local worker's lock, retaining explicit recovery evidence."""
    manifest = read_manifest(output)
    if seed not in manifest['seeds'] or arm not in ('U','H'):
        raise ValueError('unregistered arm')
    arm_dir = Path(output)/'runs'/str(seed)/arm
    lock = arm_dir/'active.lock'
    if not lock.exists():
        return False
    owner = json.loads(lock.read_text())
    if owner.get('host') != platform.node() or type(owner.get('pid')) is not int or owner['pid'] <= 0:
        raise RuntimeError('Cannot establish that this lock belongs to a dead local worker')
    try:
        os.kill(owner['pid'],0)
    except ProcessLookupError:
        recoveries = list(arm_dir.glob('recovery-*.json'))
        atomic_json(arm_dir/f'recovery-{len(recoveries)+1:03d}.json',
                    {'previous_lock':owner,'reason':'worker no longer exists; same saved input retry required'})
        lock.unlink()
        return True
    raise RuntimeError('Worker is still alive; cannot recover its lock')


def execute_arm(output, seed, arm, retry=False):
    """Worker entry point. A fresh process owns each arm, including each retry."""
    output = Path(output)
    manifest = read_manifest(output)
    if manifest['cohort'] == 'main':
        _verify_main_freeze(manifest['freeze'], manifest['protocol'], manifest['code'], manifest['runtime'])
    if (manifest['cohort'].startswith('diagnostic_v03_')
            or any('study_v03:' in str(a) for a in manifest['protocol'].get('amendments', []))):
        from .study_v03 import verify_admission
        verify_admission(output, manifest)
    pair_meta = next((p for p in manifest['pairs'] if p['seed'] == seed), None)
    if pair_meta is None or arm not in ('U', 'H'):
        raise ValueError('unregistered arm')
    inputs = read_inputs(output, seed, arm)
    if digest(inputs) != pair_meta['inputs'][arm]:
        raise ValueError('saved inputs changed')
    arm_dir = output / 'runs' / str(seed) / arm
    arm_dir.mkdir(parents=True, exist_ok=True)
    # Atomic exclusive claim prevents simultaneous retries/duplicate writers.
    lock = arm_dir / 'active.lock'
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        raise RuntimeError(f'arm {seed}/{arm} is already claimed; inspect stale worker before recovery')
    os.write(descriptor, canonical({'pid':os.getpid(), 'host':platform.node()}).encode())
    os.close(descriptor)
    attempt_dir = None
    try:
        attempts = sorted(arm_dir.glob('attempt-*'))
        previous = [json.loads((p / 'result.json').read_text()) for p in attempts if (p / 'result.json').exists()]
        completed = [r for r in previous if r.get('status') == 'completed']
        if len(completed) > 1:
            raise ValueError('Duplicate completed attempts')
        if completed:
            saved = completed[0]
            saved_dir = arm_dir / f"attempt-{saved['attempt']:03d}"
            validate_completed(saved, saved_dir, manifest, pair_meta, seed, arm)
            marker = {'attempt':saved['attempt'], 'result_sha256':digest(saved)}
            marker_path = arm_dir/'completed.json'
            if marker_path.exists() and json.loads(marker_path.read_text()) != marker:
                raise ValueError('Completion marker mismatch')
            if not marker_path.exists():
                if not retry:
                    raise RuntimeError('Interrupted completion marker requires explicit same-input retry')
                atomic_json(marker_path, marker)
            return saved
        if attempts and not retry:
            raise RuntimeError(f'{seed}/{arm} has a failed/interrupted attempt; explicit same-input retry required')
        attempt_dir = arm_dir / f'attempt-{len(attempts)+1:03d}'
        attempt_dir.mkdir(exist_ok=False)
        protocol = manifest['protocol']
        dt = protocol['clock']['dt_model_seconds']
        horizon = protocol['clock']['horizon_model_seconds']
        steps = int(Fraction(str(horizon)) / Fraction(str(dt)))
        record = {'seed': seed, 'arm': arm, 'cohort': manifest['cohort'], 'attempt': len(attempts)+1,
                  'status': 'started', 'validated': False, 'reached_time_s': 0.0,
                  'manifest_sha256': manifest['manifest_sha256'], 'input_sha256': digest(inputs),
                  'protocol_sha256': manifest['protocol_sha256'], 'code_sha256': manifest['code']['sha256']}
        atomic_json(attempt_dir / 'started.json', record)
        reached = 0.0
        sim = None
        try:
            with ExitStack() as stack:
                events = stack.enter_context(_compressed_jsonl(attempt_dir / 'events.jsonl.gz'))
                trajectory = stack.enter_context(_compressed_jsonl(attempt_dir / 'trajectory.jsonl.gz'))
                def event_sink(event):
                    events.write(canonical(_json_state(event)) + '\n')
                sim = build_simulation(inputs, event_sink=event_sink)
                feed_at = tick_schedule(inputs['feed_schedule'], dt)
                focal_name = next(s['name'] for s in protocol['species'] if s['role'] == 'focal')
                focal = next(s for s in sim.species if s.name == focal_name)
                first_extinction = None
                tail_min = None
                max_mass_residual = 0.0
                width, height = protocol['geometry']['world_width_um'], protocol['geometry']['world_height_um']
                sample_stride = max(1, int(Fraction('1') / Fraction(str(dt))))
                for tick in range(steps):
                    for event in feed_at.get(tick, []):
                        for node in event['nodes']:
                            sim.add_external_food(node['x'], node['y'], node['mass_pg'], node['substrate_id'],
                                                  event_id=node['node_id'])
                    sim.update(width, height, dt)
                    reached = float((tick+1) * Fraction(str(dt)))
                    n_focal = len(focal.members)
                    if n_focal == 0 and first_extinction is None:
                        first_extinction = reached
                    if reached >= 0.8 * horizon:
                        tail_min = n_focal if tail_min is None else min(tail_min, n_focal)
                    last_tick = _json_state(sim.experiment.last_tick)
                    max_mass_residual = max(max_mass_residual, abs(last_tick.get('mass_residual_pg', 0.0)))
                    if (tick+1) % sample_stride == 0 or tick+1 == steps:
                        sample = {'time_s': reached, 'counts': _last_counts(sim), 'tick': last_tick}
                        if hasattr(sim.experiment, 'spatial_snapshot'):
                            sample['exposure'] = _json_state(sim.experiment.spatial_snapshot(sim))
                        trajectory.write(canonical(sample) + '\n')
                        trajectory.flush()
                        events.flush()
                        atomic_json(attempt_dir / 'progress.json', {**record, 'reached_time_s': reached,
                                                                 'counts': _last_counts(sim)})
                stack.close()
                ledger = _json_state(sim.experiment.cumulative)
                record.update(status='completed', validated=True, reached_time_s=reached,
                              focal_count=len(focal.members), counts=_last_counts(sim),
                              initial_focal_count=next(s['init'] for s in protocol['species'] if s['role']=='focal'),
                              tail_min_focal_count=tail_min,
                              extinction_time_s=first_extinction,
                              established={str(t): tail_min >= t for t in (1,5,10)},
                              exposure={'populations': last_tick['population_exposure'],
                                        'early_window_seconds': [0.0, 25.0],
                                        'early_observed_seconds': last_tick['early_exposure_time_seconds']},
                              diagnostics={'max_mass_residual_pg': max_mass_residual,
                                           'food_balance_residual_pg': max_mass_residual,
                                           'exposure_contrast': last_tick['exposure_contrast'],
                                           'alive_weighted_exposure_contrast': last_tick['alive_weighted_exposure_contrast'],
                                           'early_exposure_contrast': last_tick['early_exposure_contrast'],
                                           'alive_cell_seconds': last_tick['alive_cell_seconds']}, ledger=ledger)
                # Numerical comparison also receives named scalar diagnostics without discarding the detailed ledger.
                record['diagnostics'].update({k: float(v) for k,v in ledger.items() if isinstance(v,(int,float))})
                atomic_json(attempt_dir / 'rng_final.json', {'ecology': _json_state(random.getstate()),
                           'scheduler': _json_state(sim.experiment.scheduler_rng.getstate())})
                record['outputs'] = artifact_manifest(attempt_dir)
                atomic_json(attempt_dir / 'result.json', record)
                # This marker follows validation and complete, durable outputs only.
                atomic_json(arm_dir / 'completed.json', {'attempt': record['attempt'], 'result_sha256': digest(record)})
                return record
        except BaseException as exc:
            record.update(status='failed', validated=False, reached_time_s=reached,
                          error_type=type(exc).__name__, error=str(exc), traceback=traceback.format_exc())
            if sim is not None:
                record['counts'] = _last_counts(sim)
                record['ledger'] = _json_state(getattr(sim.experiment, 'cumulative', {}))
            atomic_json(attempt_dir / 'result.json', record)
            return record
    finally:
        lock.unlink(missing_ok=True)


def _execute_pair(output, seed, order, retry):
    # This dispatcher starts a new process for every arm, so both arms begin with clean module state.
    results=[]
    for arm in order:
        with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context('spawn')) as executor:
            results.append(executor.submit(execute_arm, output, seed, arm, retry).result())
    return results


def run(output, workers=1, retry=False):
    manifest = read_manifest(output)
    if type(workers) is not int or workers < 1:
        raise ValueError('workers must be a positive integer')
    results=[]
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as executor:
        tasks={executor.submit(_execute_pair,str(output),p['seed'],p['order'],retry):p['seed'] for p in manifest['pairs']}
        for task in as_completed(tasks):
            result=task.result()
            results.extend(result)
            print(canonical({'seed':tasks[task], 'statuses':[r['status'] for r in result]}),flush=True)
    return results


def collect(output):
    """Resolve one completed attempt per registered arm; retain failures in attempt folders."""
    output=Path(output)
    manifest=read_manifest(output,check_runtime=False)
    records=[]
    expected_dirs={str(s) for s in manifest['seeds']}
    actual_dirs={p.name for p in (output/'runs').glob('*') if p.is_dir()}
    if actual_dirs - expected_dirs:
        raise ValueError('unregistered run directories')
    for pair in manifest['pairs']:
        for arm in ('U','H'):
            arm_dir=output/'runs'/str(pair['seed'])/arm
            marker=arm_dir/'completed.json'
            if not marker.exists():
                raise ValueError(f'missing completed arm {pair["seed"]}/{arm}; failures are not extinction')
            completion=json.loads(marker.read_text())
            path=arm_dir/f'attempt-{completion["attempt"]:03d}'/'result.json'
            record=json.loads(path.read_text())
            if digest(record)!=completion['result_sha256'] or record['input_sha256']!=pair['inputs'][arm]:
                raise ValueError('result/input hash mismatch')
            validate_completed(record, path.parent, manifest, pair, pair['seed'], arm)
            completed=[p for p in arm_dir.glob('attempt-*/result.json') if json.loads(p.read_text())['status']=='completed']
            if len(completed)!=1:
                raise ValueError('duplicate completed attempts')
            records.append(record)
    return manifest,records


def analyze(output):
    from .analysis import analyze_pairs
    manifest,records=collect(output)
    p=manifest['protocol']
    result=analyze_pairs(records,manifest['seeds'],cohort=manifest['cohort'],
                         horizon_s=p['clock']['horizon_model_seconds'],dt_s=p['clock']['dt_model_seconds'],
                         bootstrap_resamples=p['analysis']['bootstrap_resamples'],
                         bootstrap_seed=p['analysis']['bootstrap_seed'],confidence=p['analysis']['confidence_level'])
    result['manifest_sha256']=manifest['manifest_sha256']
    atomic_json(Path(output)/'analysis.json',result)
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    reg=sub.add_parser('register');reg.add_argument('output');reg.add_argument('--protocol');reg.add_argument('--cohort',required=True)
    reg.add_argument('--seeds',type=int,nargs='+');reg.add_argument('--freeze')
    runp=sub.add_parser('run');runp.add_argument('output');runp.add_argument('--workers',type=int,default=1);runp.add_argument('--retry',action='store_true')
    ana=sub.add_parser('analyze');ana.add_argument('output')
    recover=sub.add_parser('recover');recover.add_argument('output');recover.add_argument('--seed',type=int,required=True);recover.add_argument('--arm',choices=['U','H'],required=True)
    args=parser.parse_args(argv)
    if args.command=='register':
        result=register(args.output,load_protocol(args.protocol),args.cohort,args.seeds,freeze=args.freeze)
        print(canonical({'registered_pairs':len(result['pairs']),'manifest_sha256':result['manifest_sha256']}))
    elif args.command=='run':
        results=run(args.output,args.workers,args.retry)
        if any(r['status']!='completed' for r in results):return 1
    elif args.command=='recover':
        print(canonical({'recovered':recover_interrupted(args.output,args.seed,args.arm)}))
    else:print(json.dumps(analyze(args.output),indent=2))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
