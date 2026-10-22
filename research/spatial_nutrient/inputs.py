"""Saved founder states and exact-time paired nutrient-only inputs.

Randomness for founders, source geometry, scheduling, ecology, and analysis is
separate. The paired arms share imposed inputs, not necessarily corresponding
endogenous random events after populations diverge.
"""
from __future__ import annotations

import copy
from fractions import Fraction
import math
import random
from contextlib import contextmanager

import numpy as np

from .protocol import canonical_json, digest, exact_steps, validate_protocol

STREAM_NAMES = ('founders', 'source', 'scheduler', 'ecology', 'analysis', 'execution_order')


def _tuples(value):
    return tuple(_tuples(v) for v in value) if isinstance(value, list) else value


@contextmanager
def scoped_python_random(state):
    """Legacy constructors use module random; restore caller state even on failure.

    This is a process-local boundary and must not be used from concurrent threads.
    """
    previous = random.getstate()
    random.setstate(_tuples(state))
    try:
        yield
    finally:
        random.setstate(previous)


def seed_streams(seed):
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError('pair seed must be a nonnegative integer')
    children = np.random.SeedSequence(seed).spawn(len(STREAM_NAMES))
    streams = {}
    for name, child in zip(STREAM_NAMES, children):
        value = int(child.generate_state(1, dtype=np.uint64)[0])
        generator = np.random.default_rng(child)
        streams[name] = {
            'seed': value,
            'seed_sequence': {'entropy': seed, 'spawn_key': list(child.spawn_key)},
            'numpy_state': copy.deepcopy(generator.bit_generator.state),
            'python_state': random.Random(value).getstate(),
        }
    # JSON round trip makes tuple-based Python RNG state explicitly serializable.
    import json
    return json.loads(canonical_json(streams))


def _generator(stream):
    gen=np.random.Generator(np.random.PCG64())
    gen.bit_generator.state=copy.deepcopy(stream['numpy_state'])
    return gen


def realized_geometry(protocol):
    g=protocol['geometry']
    return {**copy.deepcopy(g), 'radius': math.hypot(g['world_width_um'],g['world_height_um'])/2+g['padding_um']}


def _species(protocol, order=('focal','resident')):
    from bacteria.core_sim import BacteriaSpecies
    specs={s['role']:s for s in protocol['species']}
    if len(order)!=2 or set(order)!=set(specs):
        raise ValueError('species_order must contain focal and resident once each')
    colors={'focal':(230,80,80),'resident':(80,130,230)}
    result=[]
    for role in order:
        kwargs={k:v for k,v in specs[role].items() if k!='role'}
        sp=BacteriaSpecies(**kwargs,color=colors[role])
        sp.population_role=role
        sp.maintenance_aging=protocol['lifecycle']['maintenance_aging']
        sp.daughter_size_policy=protocol['lifecycle']['daughter_size_policy']
        result.append(sp)
    return result


def founder_snapshot(protocol, geometry, stream):
    """Construct all positions before random labels; no source angle enters here."""
    from bacteria.core_sim import Bacterium
    species={sp.population_role:sp for sp in _species(protocol)}
    n=sum(s['init'] for s in protocol['species'])
    n_focal=next(s['init'] for s in protocol['species'] if s['role']=='focal')
    rng=_generator(stream)
    radius=geometry['radius']*np.sqrt(rng.random(n))
    angles=rng.uniform(0,2*math.pi,n)
    positions=list(zip(geometry['center'][0]+radius*np.cos(angles),geometry['center'][1]+radius*np.sin(angles)))
    focal=set(map(int,rng.choice(n,size=n_focal,replace=False)))
    snapshots=[]
    with scoped_python_random(stream['python_state']):
        for agent_id,(x,y) in enumerate(positions):
            role='focal' if agent_id in focal else 'resident'
            b=Bacterium(species[role],float(x),float(y),size=protocol['founders']['size'],
                        initial_food=protocol['founders']['initial_food_pg'])
            b.age=protocol['founders']['age']
            b.agent_id=agent_id
            state=copy.deepcopy({k:v for k,v in b.__dict__.items() if k!='species'})
            snapshots.append({'agent_id':agent_id,'population':role,'state':state})
    # Strict serialization rejects hidden live objects or nonfinite states.
    canonical_json(snapshots)
    return snapshots


def _feed_schedule(protocol, geometry, stream):
    rng=_generator(stream)
    angle=float(rng.uniform(0,2*math.pi))
    cfg=protocol['feeding']
    count=exact_steps(protocol['clock']['horizon_model_seconds'],cfg['interval_s'],'feed batches')
    schedule=[]
    for index in range(count):
        nodes=[]
        for slot in range(4):
            r=geometry['radius']*math.sqrt(float(rng.random()))
            theta=angle+float(rng.uniform(-math.pi/2,math.pi/2))+(0 if slot<2 else math.pi)
            nodes.append({'node_id':index*4+slot,
                          'x':geometry['center'][0]+r*math.cos(theta),
                          'y':geometry['center'][1]+r*math.sin(theta),
                          'mass_pg':cfg['mass_pg_each'],'substrate_id':cfg['substrate_id']})
        order=rng.permutation(4)
        physical=Fraction(index)*Fraction(str(cfg['interval_s']))
        schedule.append({'event_index':index,'time_s':float(physical),
                         'nodes':[nodes[int(i)] for i in order]})
    return angle,schedule


def make_pair_inputs(protocol, seed):
    p=validate_protocol(protocol)
    streams=seed_streams(seed)
    geometry=realized_geometry(p)
    founders=founder_snapshot(p,geometry,streams['founders'])
    angle,schedule=_feed_schedule(p,geometry,streams['source'])
    pair={'pair_id':seed,'protocol':p,'geometry':geometry,'founders':founders,
          'feed_schedule':schedule,'streams':streams,'orientation_rad':angle}
    pair['input_hash']=digest(pair)
    return pair


def _reflect(x,y,angle,center):
    nx,ny=math.cos(angle),math.sin(angle)
    dot=(x-center[0])*nx+(y-center[1])*ny
    return x-2*dot*nx,y-2*dot*ny


def arm_inputs(pair, arm):
    """Derive both treatments from one saved input object, without new RNG draws."""
    if arm not in ('U','H'):
        raise ValueError('arm must be U or H')
    source={k:v for k,v in pair.items() if k!='input_hash'}
    if pair.get('input_hash')!=digest(source):
        raise ValueError('Pair input hash does not match saved contents')
    out=copy.deepcopy(source)
    out['arm']=arm
    if arm=='H':
        moved=out['protocol']['feeding']['contrast_nodes']
        for event in out['feed_schedule']:
            for node in event['nodes']:
                if 2 <= node['node_id']%4 < 2+moved:
                    node['x'],node['y']=_reflect(node['x'],node['y'],out['orientation_rad'],out['geometry']['center'])
    out['input_hash']=digest(out)
    return out


def _check_arm(arm):
    keys={'pair_id','protocol','geometry','founders','feed_schedule','streams','orientation_rad','arm','input_hash'}
    if set(arm)!=keys:
        raise ValueError(f'Unexpected arm input fields: {set(arm)^keys}')
    if arm['input_hash']!=digest({k:v for k,v in arm.items() if k!='input_hash'}):
        raise ValueError('Arm input hash does not match contents')
    validate_protocol(arm['protocol'])
    if arm['geometry'] != realized_geometry(arm['protocol']):
        raise ValueError('Realized geometry differs from protocol')
    if arm['arm'] not in ('U','H'):
        raise ValueError('Invalid arm label')
    if not isinstance(arm['pair_id'],int) or isinstance(arm['pair_id'],bool) or arm['pair_id']<0:
        raise ValueError('Invalid pair id')
    if arm['streams']!=seed_streams(arm['pair_id']):
        raise ValueError('Saved streams differ from registered seed')
    canonical_json(arm)
    expected_pair=make_pair_inputs(arm['protocol'],arm['pair_id'])
    reconstructed=arm_inputs(expected_pair,arm['arm'])
    if canonical_json(reconstructed)!=canonical_json(arm):
        raise ValueError('Saved arm inputs differ from declared protocol and registered streams')
    expected=sum(s['init'] for s in arm['protocol']['species'])
    founders=arm['founders']
    if len(founders)!=expected or [f['agent_id'] for f in founders]!=list(range(expected)):
        raise ValueError('Founder identities must be complete and ordered')
    for role in ('focal','resident'):
        target=next(s['init'] for s in arm['protocol']['species'] if s['role']==role)
        if sum(f['population']==role for f in founders)!=target:
            raise ValueError('Founder population labels do not match protocol')
    for f in founders:
        if f['state'].get('agent_id')!=f['agent_id']:
            raise ValueError('Founder stable identity mismatch')
    schedule=arm['feed_schedule']
    batches=exact_steps(arm['protocol']['clock']['horizon_model_seconds'],arm['protocol']['feeding']['interval_s'])
    if len(schedule)!=batches or [e['event_index'] for e in schedule]!=list(range(batches)):
        raise ValueError('Feed event manifest is incomplete')
    for index,event in enumerate(schedule):
        if event['time_s']!=float(Fraction(index)*Fraction(str(arm['protocol']['feeding']['interval_s']))):
            raise ValueError('Feed physical event time changed')
        if len(event['nodes'])!=4 or {n['node_id'] for n in event['nodes']}!=set(range(index*4,index*4+4)):
            raise ValueError('Feed node identities changed')
        for node in event['nodes']:
            if node['mass_pg']!=arm['protocol']['feeding']['mass_pg_each'] or node['substrate_id']!=arm['protocol']['feeding']['substrate_id']:
                raise ValueError('Feed dose or composition differs from protocol')
            if math.hypot(node['x']-arm['geometry']['center'][0],node['y']-arm['geometry']['center'][1])>arm['geometry']['radius']*(1+1e-14):
                raise ValueError('Food node lies outside the domain')


def assert_matched_arms(uniform, heterogeneous):
    """Permit only declared half-disc coordinate reflections, never post-run matching."""
    _check_arm(uniform);_check_arm(heterogeneous)
    if uniform['arm']!='U' or heterogeneous['arm']!='H':
        raise ValueError('Expected U then H arm')
    u=copy.deepcopy(uniform);h=copy.deepcopy(heterogeneous)
    for value in (u,h):
        del value['arm'];del value['input_hash']
    # Transform the proposed U coordinates into exactly the permitted H coordinates.
    moved=u['protocol']['feeding']['contrast_nodes']
    for event in u['feed_schedule']:
        for node in event['nodes']:
            if 2 <= node['node_id']%4 < 2+moved:
                node['x'],node['y']=_reflect(node['x'],node['y'],u['orientation_rad'],u['geometry']['center'])
    if canonical_json(u)!=canonical_json(h):
        def differing(a,b,path='inputs'):
            if isinstance(a,dict) and isinstance(b,dict):
                for key in sorted(set(a)|set(b)):
                    if key not in a or key not in b:return f'{path}.{key}'
                    found=differing(a[key],b[key],f'{path}.{key}')
                    if found:return found
                return None
            if isinstance(a,list) and isinstance(b,list):
                if len(a)!=len(b):return path+'.length'
                for i,(av,bv) in enumerate(zip(a,b)):
                    found=differing(av,bv,f'{path}[{i}]')
                    if found:return found
                return None
            return path if a!=b else None
        raise ValueError(f'Forbidden between-arm difference at {differing(u,h)}')
    return True


def tick_schedule(schedule,dt):
    """Map precomputed physical feed times to exact ticks, retaining event identities."""
    result={}
    for event in schedule:
        tick=exact_steps(event['time_s'],dt,'feed time')
        result.setdefault(tick,[]).append(copy.deepcopy(event))
    return result


def restore_founders(sim,snapshot):
    """Restore every saved field without invoking constructors or taking RNG draws."""
    from bacteria.core_sim import Bacterium
    roles={sp.population_role:sp for sp in sim.species}
    for sp in sim.species:
        sp.members=[]
    for record in snapshot:
        sp=roles[record['population']]
        b=Bacterium.__new__(Bacterium)
        b.__dict__.update(copy.deepcopy(record['state']))
        # JSON arrays replace tuples; preserve the constructor's immutable trait type.
        if 'substrate_affinity' in b.__dict__:
            b.substrate_affinity=tuple(b.substrate_affinity)
        b.species=sp
        sp.members.append(b)
    sim._populated=True


def build_simulation(arm,event_sink=None,species_order=('focal','resident')):
    """Build a Python-only simulation from saved inputs. Feeding is runner-owned.

    The legacy ecological RNG is process-global; the runner must execute each arm
    in an isolated process (or sequentially without interleaving simulations).
    """
    from bacteria.core_sim import Simulation
    from bacteria.constants import PETRI_RADIUS_PADDING_UM
    from bacteria.experiment_support import validate_lifecycle
    _check_arm(arm)
    p=arm['protocol']
    if PETRI_RADIUS_PADDING_UM!=p['geometry']['padding_um']:
        raise ValueError('Runtime geometry padding differs from pinned protocol')
    spp=_species(p,species_order)
    for sp in spp:
        validate_lifecycle(sp,p['clock']['dt_model_seconds'],p['founders']['initial_food_pg'],p['founders']['size'])
    env=copy.deepcopy(p['environment_tick_start_clamp'])
    env['drop_frequency']=0
    sim=Simulation(spp,env,p['capacity_cells'],steps=p['clock']['steps'],seed=None,
                   physics_backend='python',ecm_backend='python')
    sim.petri_center=tuple(arm['geometry']['center'])
    sim.petri_radius=arm['geometry']['radius']
    restore_founders(sim,arm['founders'])
    config={
        'maintenance_aging':p['lifecycle']['maintenance_aging'],
        'daughter_size_policy':p['lifecycle']['daughter_size_policy'],
        'recycling':p['lifecycle']['recycling'],
        'washout_rate':p['neutralized_radial_parameters']['LUMEN_FOOD_WASHOUT_RATE_PER_SEC'],
        'motility_factor':p['neutralized_radial_parameters']['MUCUS_MOTILITY_FACTOR'],
        'source_angle':arm['orientation_rad'],
        'mass_tolerance':p['execution']['mass_tolerance'],
        'max_food_nodes':p['execution']['max_food_nodes'],
        'environment_clamp':p['environment_tick_start_clamp'],
    }
    sim.enable_experiment(config,arm['streams']['scheduler']['seed'],event_sink=event_sink)
    random.setstate(_tuples(arm['streams']['ecology']['python_state']))
    sim.running=True
    return sim
