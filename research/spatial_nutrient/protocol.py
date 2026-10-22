"""Strict, versioned inputs for the synthetic legacy-pg experiment.

All numerical settings are Own-model choices, not measured species parameters.
See DOCUMENTATION.md, spatial nutrient experiment. Validation deliberately rejects
unknown fields and silent legacy clamping before a registered run starts.
"""
from __future__ import annotations

import copy
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

DEFAULT_PROTOCOL = Path(__file__).with_name('protocol_v02.json')


def canonical_json(value):
    """Portable deterministic JSON; nonfinite numerical values are never records."""
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical_json(value).encode('utf-8')).hexdigest()


def _number(value, name, minimum=None, maximum=None, integer=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{name} must be a finite number')
    if not math.isfinite(value) or (integer and type(value) is not int):
        raise ValueError(f'{name} must be finite' + (' and integral' if integer else ''))
    if minimum is not None and value < minimum:
        raise ValueError(f'{name} must be >= {minimum}')
    if maximum is not None and value > maximum:
        raise ValueError(f'{name} must be <= {maximum}')


def _shape(value, template, path='protocol'):
    if isinstance(template, dict):
        if not isinstance(value, dict):
            raise ValueError(f'{path} must be an object')
        extra, missing = set(value)-set(template), set(template)-set(value)
        if extra or missing:
            raise ValueError(f'{path}: unknown fields {sorted(extra)}; missing fields {sorted(missing)}')
        for key in template:
            _shape(value[key], template[key], f'{path}.{key}')
    elif isinstance(template, list):
        if not isinstance(value, list):
            raise ValueError(f'{path} must be a list')
        if template:
            for i, item in enumerate(value):
                _shape(item, template[0], f'{path}[{i}]')
    elif isinstance(template, bool):
        if type(value) is not bool:
            raise ValueError(f'{path} must be boolean')
    elif isinstance(template, (int, float)):
        _number(value, path)
    elif isinstance(template, str) and not isinstance(value, str):
        raise ValueError(f'{path} must be text')


def exact_steps(duration, dt, name='duration'):
    """Convert physical seconds to ticks without floating-point remainder tests."""
    _number(duration, name, 0)
    _number(dt, 'dt', 0)
    if dt == 0:
        raise ValueError('dt must be positive')
    ratio = Fraction(str(duration)) / Fraction(str(dt))
    if ratio.denominator != 1:
        raise ValueError(f'{name}={duration} is not an exact multiple of dt={dt}')
    return ratio.numerator


def validate_protocol(protocol):
    """Return an independent normalized configuration; never fill missing settings."""
    template = json.loads(DEFAULT_PROTOCOL.read_text())
    p = copy.deepcopy(protocol)
    _shape(p, template)
    canonical_json(p)
    for key, expected in [('protocol_id','SPATIAL-NUTRIENT-001'), ('version','0.2'),
                          ('engine','legacy_pg_experiment'), ('backend','python')]:
        if p[key] != expected:
            raise ValueError(f'{key} must be {expected!r}')
    if p['status'] not in ('qualification_required', 'frozen'):
        raise ValueError('status must be qualification_required or frozen')
    if len(p['species']) != 2 or {s['role'] for s in p['species']} != {'focal','resident'}:
        raise ValueError('Exactly one focal and one resident species are required')
    if len({s['name'] for s in p['species']}) != 2:
        raise ValueError('Species names must be distinct')
    p['species'].sort(key=lambda s: s['role'])
    disabled = ('aging','so','t6ss_activity','bacteriocin_production','immunity_factor',
                'predation_thresh','acetate_secretion','acetate_uptake','persister_switch',
                'phage_resistance','hgt_competence','aff_glycan','aff_mdcf')
    unit_fields = ('greediness','t6ss_activity','bacteriocin_production','immunity_factor',
                   'rep_thresh','share_thresh','predation_thresh','acetate_secretion',
                   'acetate_uptake','persister_switch','persister_enter_thresh',
                   'persister_exit_thresh','phage_resistance','hgt_competence')
    for s in p['species']:
        for key, value in s.items():
            if key not in ('name','role'):
                _number(value, f"species.{s['role']}.{key}")
        for key in disabled:
            if s[key] != 0:
                raise ValueError(f'species.{s["role"]}.{key} must be disabled (zero)')
        for key in unit_fields:
            _number(s[key], key, 0, 1)
        for key in ('rep','food_capacity','metabolic_rate','consumption_speed','speed','sensor_radius','aff_calorie'):
            _number(s[key], key, 0)
        _number(s['init'],'init',1,integer=True)
        _number(s['base_size'],'base_size',0.1)
        _number(s['ph'],'ph',0,14)
        _number(s['baranyi_q0'],'baranyi_q0',1e-4,1e3)
        if s['temp_max'] < s['temp_min'] + 1:
            raise ValueError('temp_max must be at least temp_min + 1')
        if s['persister_exit_thresh'] < s['persister_enter_thresh'] + 0.05:
            raise ValueError('persister_exit_thresh must be at least enter threshold + .05')
        if s['share_thresh'] != 1 or s['sf'] != 1:
            raise ValueError('Sharing must be disabled (share_thresh=1); sf must equal 1')
    f = p['founders']
    _number(f['size'],'founders.size',0.1)
    _number(f['age'],'founders.age',0)
    _number(f['initial_food_pg'],'founders.initial_food_pg',0)
    for s in p['species']:
        if s['base_size'] != f['size']:
            raise ValueError('Species base_size must match the declared fixed founder size')
        capacity = s['food_capacity'] * f['size']**3
        if f['initial_food_pg'] > capacity:
            raise ValueError('Founder reserves exceed capacity')
        if capacity <= 40.0*f['size']**3 or s['rep_thresh'] >= 1 or s['rep'] <= 0:
            raise ValueError('Lifecycle infeasible: capacity, fission cost, and threshold cannot enable division')
    g = p['geometry']
    for key in ('world_width_um','world_height_um'):
        _number(g[key],key,1)
    _number(g['padding_um'],'padding_um',0)
    if g['center'] != [0.0,0.0] or g['padding_um'] != 100.0:
        raise ValueError('Geometry must match legacy origin and 100 um padding')
    c = p['clock']
    _number(c['horizon_model_seconds'],'horizon_model_seconds',0)
    if c['horizon_model_seconds'] <= 0:
        raise ValueError('horizon must be positive')
    steps = exact_steps(c['horizon_model_seconds'],c['dt_model_seconds'],'horizon')
    _number(c['steps'],'steps',1,integer=True)
    if c['steps'] != steps:
        raise ValueError('steps does not match the exact physical clock')
    _number(p['capacity_cells'],'capacity_cells',1,integer=True)
    if p['capacity_cells'] < sum(s['init'] for s in p['species']):
        raise ValueError('Founders exceed carrying capacity')
    feed = p['feeding']
    if feed['nodes_per_event'] != 4 or feed['substrate_id'] != 0:
        raise ValueError('Experiment requires four generic-calorie nodes per feed event')
    if feed['first_event_s'] != 0 or feed['automatic_drop_frequency'] != 0:
        raise ValueError('Scheduled feed starts at zero; automatic feeding must be disabled')
    _number(feed['mass_pg_each'],'mass_pg_each',0)
    _number(feed['contrast_nodes'],'contrast_nodes',0,2,integer=True)
    if feed['interval_s'] <= 0:
        raise ValueError('feed interval must be positive')
    exact_steps(feed['interval_s'],c['dt_model_seconds'],'feed interval')
    exact_steps(c['horizon_model_seconds'],feed['interval_s'],'horizon/feed interval')
    radial = p['neutralized_radial_parameters']
    for key in ('MUCUS_MOTILITY_FACTOR','MUCUS_DROP_AMOUNT_FACTOR','LUMEN_COMPETITION_CULL_WEIGHT'):
        if radial[key] != 1:
            raise ValueError(f'{key} must equal one in the nutrient-only experiment')
    for key in ('LUMEN_FOOD_WASHOUT_RATE_PER_SEC','MUCUS_FOOD_WASHOUT_RATE_PER_SEC'):
        _number(radial[key],key,0)
    if radial['LUMEN_FOOD_WASHOUT_RATE_PER_SEC'] != radial['MUCUS_FOOD_WASHOUT_RATE_PER_SEC']:
        raise ValueError('Washout must be spatially neutral')
    if radial['LUMEN_FOOD_WASHOUT_RATE_PER_SEC']*c['dt_model_seconds'] > 1:
        raise ValueError('Washout Euler step can remove more than available food')
    env = p['environment_tick_start_clamp']
    if env['phage_inflow'] != 0:
        raise ValueError('Phage inflow must be disabled')
    _number(env['ph'],'environment.ph',0,14)
    _number(env['o2'],'environment.o2',0,1)
    _number(env['co2'],'environment.co2',0,1)
    life=p['lifecycle']
    if life['daughter_size_policy'] not in ('variable','fixed') or life['overflow_policy'] != 'return_to_parent':
        raise ValueError('Unsupported daughter size or overflow policy')
    ex=p['execution']
    if ex['scheduler'] != 'global_random_agent' or ex['birth_timing'] != 'next_tick' or ex['feed_timing'] != 'before_washout_and_uptake':
        raise ValueError('Unsupported execution policy')
    if ex['max_food_nodes'] is not None:
        _number(ex['max_food_nodes'],'max_food_nodes',1,integer=True)
    _number(ex['mass_tolerance'],'mass_tolerance',0)
    if ex['mass_tolerance']==0:
        raise ValueError('mass_tolerance must be positive')
    seen=set()
    for group in ('pilot_seeds','main_seeds','numerical_seeds','numerical_extension_seeds'):
        seeds=p['replication'][group]
        if not seeds or len(set(seeds)) != len(seeds):
            raise ValueError(f'{group} must be nonempty and unique')
        for seed in seeds:
            _number(seed,group,0,integer=True)
        if seen.intersection(seeds):
            raise ValueError('Pilot, main, and numerical cohorts must be disjoint')
        seen.update(seeds)
    _number(p['replication']['analysis_seed'],'analysis_seed',0,integer=True)
    analysis=p['analysis']
    if analysis['bootstrap_unit']!='whole pair' or analysis['bootstrap_method']!='percentile':
        raise ValueError('Paired percentile resampling is required')
    _number(analysis['bootstrap_resamples'],'bootstrap_resamples',1,integer=True)
    _number(analysis['bootstrap_seed'],'bootstrap_seed',0,integer=True)
    _number(analysis['confidence_level'],'confidence_level',0,1)
    if not 0 < analysis['confidence_level'] < 1:
        raise ValueError('confidence_level must be between zero and one')
    gates=p['numerical_gates']
    exact_steps(c['horizon_model_seconds'],gates['refined_dt'],'refined horizon')
    exact_steps(feed['interval_s'],gates['refined_dt'],'refined feed interval')
    if gates['refined_dt'] > c['dt_model_seconds']:
        raise ValueError('refined_dt must not be coarser than execution dt')
    for key in ('effect_tolerance','confidence_level','max_pruned_fraction_of_sources','min_realized_food_contrast'):
        _number(gates[key],key,0,1)
    if gates['cumulative_pair_counts'] != [50,200,1000] or gates['initial_pairs'] != 50 or gates['seed_start'] != 2000:
        raise ValueError('Numerical qualification uses prespecified cumulative looks at 50, 200, and 1000 pairs starting seed 2000')
    if p['replication']['numerical_seeds'] != list(range(2000,2050)) or p['replication']['numerical_extension_seeds'] != list(range(2050,3000)):
        raise ValueError('Numerical seed registry does not match the staged qualification plan')
    if p['replication']['analysis_seed'] != analysis['bootstrap_seed']:
        raise ValueError('Analysis seed declarations disagree')
    # These fields describe the implemented analysis; changing them requires a
    # versioned implementation amendment rather than silently unused settings.
    for key in ('primary','extinction','censoring','initial_zero','secondary_establishment','secondary_threshold_checks'):
        if p['outcomes'][key] != template['outcomes'][key]:
            raise ValueError(f'Unsupported outcome definition: {key}')
    for key in ('bootstrap_unit','paired_risk_difference','conservative_interval','per_arm_interval','optional_test','failed_runs'):
        if analysis[key] != template['analysis'][key]:
            raise ValueError(f'Unsupported analysis definition: {key}')
    return p


def load_protocol(path=None):
    return validate_protocol(json.loads(Path(path or DEFAULT_PROTOCOL).read_text()))
