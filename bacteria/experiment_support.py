"""Opt-in integrity support for declared synthetic spatial experiments.

All lifecycle policies here are Own (BacSandbox model) assumptions, not measured
microbial parameters. See DOCUMENTATION.md, spatial experiment qualification.
"""
import math

from .constants import FISSION_COST_PG_PER_UM3


def validate_lifecycle(species, dt, initial_food, size=1.0):
    """Reject unreachable division before a birth-and-death experiment starts.

Storage and initial reserve are independent pg glucose-equivalent values.
This establishes algebraic feasibility; the isolated fed-cell qualification
must separately exercise the actual uptake and division code.
"""
    values = (dt, size, initial_food, species.food_capacity,
              species.metabolic_rate, species.rep_thresh, species.base_rep)
    if any(isinstance(v, bool) or not math.isfinite(v) for v in values):
        raise ValueError('Lifecycle parameters must be finite numbers')
    if dt <= 0 or size <= 0 or species.metabolic_rate < 0 or species.base_rep <= 0:
        raise ValueError('Lifecycle requires positive dt, size and reproduction rate')
    capacity = species.food_capacity * size ** 3
    cost = FISSION_COST_PG_PER_UM3 * size ** 3
    survival = species.metabolic_rate * size ** 2.25 * dt
    threshold = max(capacity * species.rep_thresh, cost + 2 * survival)
    if not 0 <= initial_food <= capacity:
        raise ValueError('Initial reserve must be within capacity')
    if not 0 <= species.rep_thresh < 1 or not capacity > threshold:
        raise ValueError('Division unreachable within storage capacity')
    return {'capacity_pg': capacity, 'initial_food_pg': initial_food,
            'fission_cost_pg': cost, 'minimum_division_reserve_pg': threshold,
            'strictly_greater_required': True}


class ExperimentIntegrityError(RuntimeError):
    """A registered run failed an explicit numerical or configuration contract."""


class ExperimentLedger:
    """Event-time external-mass, reserve and population accounting.

    Events are streamed to ``event_sink`` when supplied; only one tick is kept
    in memory. Substrate mass and internal glucose-equivalent reserves have
    separate identities. Necromass is an explicit external-pool source, not a
    claim that this abstraction conserves all cellular elemental mass.
    """
    def __init__(self, sim, config, scheduler_seed, event_sink=None):
        import random
        self.sim, self.config, self.event_sink = sim, dict(config), event_sink
        self.scheduler = random.Random(scheduler_seed)
        self.scheduler_rng = self.scheduler
        self.exposure_contrast_integral = 0.0
        self.alive_weighted_contrast_integral = 0.0
        self.alive_cell_seconds = 0.0
        self.elapsed_seconds = 0.0
        self.cumulative = {k: 0.0 for k in ('feed_pg', 'recycled_pg', 'uptake_pg',
            'washout_pg', 'pruned_pg', 'reserve_gain_pg', 'maintenance_pg',
            'fission_cost_pg', 'death_reserve_pg', 'maintenance_demand_pg',
            'unmet_maintenance_pg')}
        self.cumulative.update(births=0, deaths=0)
        self.events, self.last_events, self.last_tick = [], [], None
        self.population_events = {}
        self.residence = {}
        self.last_mass = sum(n.amount for n in sim.food_nodes)
        self.last_reserves = self.reserves()
        self.last_counts = self.counts()
        self.start_cumulative = dict(self.cumulative)
        self.next_agent_id = 0
        seen = set()
        for sp in sorted(sim.species, key=lambda s: s.name):
            for b in sp.members:
                aid = getattr(b, 'agent_id', None)
                if aid is None:
                    while self.next_agent_id in seen:
                        self.next_agent_id += 1
                    aid = b.agent_id = self.next_agent_id
                if type(aid) is not int or aid < 0 or aid in seen:
                    raise ValueError('Experiment requires unique nonnegative integer agent_id')
                seen.add(aid)
                self.next_agent_id = max(self.next_agent_id, aid + 1)
        self.angle = float(config.get('source_angle', 0))
        if not math.isfinite(self.angle):
            raise ValueError('source_angle must be finite')
        if config.get('motility_factor', 1) != 1:
            raise ValueError('Experimental radial motility factor must be one')
        self.tolerance = float(config.get('mass_tolerance', 1e-9))
        if not math.isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError('mass_tolerance must be finite and positive')
        self.washout_rate = float(config.get('washout_rate', .02))
        self.max_food_nodes = config.get('max_food_nodes')
        self.environment_clamp = dict(config.get('environment_clamp', {}))
        if not math.isfinite(self.washout_rate) or self.washout_rate < 0:
            raise ValueError('washout_rate must be finite and nonnegative')
        if self.max_food_nodes is not None and (type(self.max_food_nodes) is not int or self.max_food_nodes < 1):
            raise ValueError('max_food_nodes must be positive integer or None')
        if not isinstance(config.get('recycling', True), bool):
            raise ValueError('recycling must be boolean')
        self.recycling = config.get('recycling', True)
        for sp in sim.species:
            if any(getattr(sp, name) != 0 for name in (
                't6ss_activity', 'bacteriocin_production', 'acetate_secretion',
                'acetate_uptake', 'persister_switch', 'hgt_competence', 'predation_thresh')):
                raise ValueError('Experimental ledger requires disabled interaction pathways')
            if sp.share_thresh != 1:
                raise ValueError('Experimental ledger requires share_thresh=1')
            if any(b.infected for b in sp.members):
                raise ValueError('Experimental ledger requires uninfected founders')
            sp.maintenance_aging = config.get('maintenance_aging', sp.maintenance_aging)
            sp.daughter_size_policy = config.get('daughter_size_policy', sp.daughter_size_policy)
            if not isinstance(sp.maintenance_aging, bool):
                raise ValueError('maintenance_aging must be boolean')
            if sp.daughter_size_policy not in ('variable', 'fixed'):
                raise ValueError('daughter_size_policy must be variable or fixed')
        if sim.drop_frequency != 0 or sim.env.get('phage_inflow', 0) != 0:
            raise ValueError('Experiment requires scheduled-only feeding and no phage inflow')

    def counts(self):
        return {s.name: len(s.members) for s in self.sim.species}

    def reserves(self):
        return {s.name: math.fsum(b.food_level for b in s.members) for s in self.sim.species}

    def region(self, x, y):
        cx, cy = self.sim.petri_center
        return 'rich' if (x - cx) * math.cos(self.angle) + (y - cy) * math.sin(self.angle) >= 0 else 'poor'

    def record(self, kind, agent=None, amount_pg=None, **fields):
        event = {'tick': self.sim.step_idx, 'kind': kind, **fields}
        if agent is not None:
            event.update(agent_id=agent.agent_id, population=agent.species.name,
                         x=agent.x, y=agent.y, region=self.region(agent.x, agent.y))
        if amount_pg is not None:
            if not math.isfinite(amount_pg) or amount_pg < 0:
                raise ExperimentIntegrityError('Nonfinite or negative resource transfer')
            event['amount_pg'] = amount_pg
        keys = {'feed':'feed_pg', 'recycling':'recycled_pg', 'uptake':'uptake_pg',
                'washout':'washout_pg', 'prune':'pruned_pg', 'maintenance':'maintenance_pg',
                'fission':'fission_cost_pg', 'death_reserve':'death_reserve_pg'}
        if kind in keys:
            self.cumulative[keys[kind]] += amount_pg
        if kind == 'maintenance':
            demand = fields.get('demand_pg', amount_pg)
            if not math.isfinite(demand) or demand < amount_pg:
                raise ExperimentIntegrityError('Invalid maintenance demand')
            unmet = max(0.0, demand - amount_pg)
            self.cumulative['maintenance_demand_pg'] += demand
            self.cumulative['unmet_maintenance_pg'] += unmet
            event['unmet_maintenance_pg'] = unmet
        if kind == 'uptake':
            self.cumulative['reserve_gain_pg'] += fields['reserve_gain_pg']
            agent.experiment_uptake_pg = getattr(agent, 'experiment_uptake_pg', 0.0) + amount_pg
        if kind in ('birth', 'death'):
            self.cumulative['births' if kind == 'birth' else 'deaths'] += 1
        self.events.append(event)
        if self.event_sink is not None:
            self.event_sink(event)

    def begin_tick(self, dt):
        if not math.isfinite(dt) or dt <= 0 or self.washout_rate * dt > 1:
            raise ExperimentIntegrityError('Invalid dt or Euler washout outside positivity range')
        self.sim.env.update(self.environment_clamp)
        if self.sim.drop_frequency != 0 or self.sim.env.get('phage_inflow', 0) != 0:
            raise ExperimentIntegrityError('Disabled automatic input mechanism was enabled')
        self.tick_start_mass = self.last_mass
        self.tick_start_reserves = dict(self.last_reserves)
        self.tick_start_counts = dict(self.last_counts)
        for sp in self.sim.species:
            for b in sp.members:
                r = self.region(b.x, b.y)
                k = sp.name + ':' + r
                self.residence[k] = self.residence.get(k, 0.0) + dt
                attribute = 'time_in_' + r + '_seconds'
                setattr(b, attribute, getattr(b, attribute, 0.0) + dt)

    def finish_tick(self, dt):
        delta = {k: self.cumulative[k] - self.start_cumulative[k] for k in self.cumulative}
        mass = math.fsum(n.amount for n in self.sim.food_nodes)
        expected = (self.tick_start_mass + delta['feed_pg'] + delta['recycled_pg']
                    - delta['uptake_pg'] - delta['washout_pg'] - delta['pruned_pg'])
        scale = max(1.0, self.tick_start_mass + delta['feed_pg'] + delta['recycled_pg'])
        if not math.isfinite(mass) or abs(mass - expected) > self.tolerance * scale:
            raise ExperimentIntegrityError('External food mass ledger does not reconcile')
        if any(not math.isfinite(n.amount) or n.amount < 0 for n in self.sim.food_nodes):
            raise ExperimentIntegrityError('Invalid food node amount')
        counts, reserves = self.counts(), self.reserves()
        population = {}
        for sp in self.sim.species:
            events = [e for e in self.events if e.get('population') == sp.name]
            births = sum(e['kind'] == 'birth' for e in events)
            deaths = {}
            for e in events:
                if e['kind'] == 'death':
                    deaths[e['cause']] = deaths.get(e['cause'], 0) + 1
            if counts[sp.name] != self.tick_start_counts[sp.name] + births - sum(deaths.values()):
                raise ExperimentIntegrityError('Population ledger does not reconcile')
            gain = math.fsum(e.get('reserve_gain_pg', 0) for e in events if e['kind'] == 'uptake')
            loss = math.fsum(e['amount_pg'] for e in events if e['kind'] in ('maintenance', 'fission', 'death_reserve'))
            reserve_expected = self.tick_start_reserves[sp.name] + gain - loss
            if abs(reserves[sp.name] - reserve_expected) > self.tolerance * max(1, self.tick_start_reserves[sp.name] + gain):
                raise ExperimentIntegrityError('Internal reserve ledger does not reconcile')
            if any(not math.isfinite(b.food_level) or b.food_level < 0 or b.food_level > b.food_capacity + self.tolerance for b in sp.members):
                raise ExperimentIntegrityError('Invalid internal reserve')
            population[sp.name] = {'start': self.tick_start_counts[sp.name], 'end': counts[sp.name],
                'births': births, 'deaths': deaths, 'reserve_pg': reserves[sp.name],
                'reserve_residual_pg': reserves[sp.name] - reserve_expected}
        exposure = self.exposure()
        rich, poor = exposure['food_half_pg']['rich'], exposure['food_half_pg']['poor']
        contrast = (rich - poor) / (rich + poor) if rich + poor else 0.0
        self.exposure_contrast_integral += dt * contrast
        # This weights field contrast by living tick-start abundance; it still
        # does not measure each agent's local exposure (see uptake/location events).
        alive_weight = dt * sum(self.tick_start_counts.values())
        self.alive_weighted_contrast_integral += alive_weight * contrast
        self.alive_cell_seconds += alive_weight
        self.elapsed_seconds += dt
        self.last_tick = {'tick': self.sim.step_idx, 'time_seconds': self.elapsed_seconds,
            'mass_start_pg': self.tick_start_mass, 'mass_end_pg': mass,
            'mass_residual_pg': mass - expected, **delta, 'populations': population,
            'exposure': exposure, 'exposure_contrast': self.exposure_contrast_integral / self.elapsed_seconds,
            'alive_weighted_exposure_contrast': (self.alive_weighted_contrast_integral / self.alive_cell_seconds
                                                if self.alive_cell_seconds else None),
            'alive_cell_seconds': self.alive_cell_seconds}
        self.last_mass, self.last_reserves, self.last_counts = mass, reserves, counts
        self.start_cumulative = dict(self.cumulative)
        self.last_events, self.events = self.events, []
        return self.last_tick

    def exposure(self):
        halves = {'rich': 0.0, 'poor': 0.0}
        bins = [0.0] * 16
        cx, cy = self.sim.petri_center
        radius = self.sim.petri_radius
        for node in self.sim.food_nodes:
            halves[self.region(node.x, node.y)] += node.amount
            ix = max(0, min(3, int((node.x - cx + radius) / (2 * radius) * 4)))
            iy = max(0, min(3, int((node.y - cy + radius) / (2 * radius) * 4)))
            bins[4 * iy + ix] += node.amount
        populations = {}
        for sp in self.sim.species:
            occupancy = {'rich': 0, 'poor': 0}
            for b in sp.members:
                occupancy[self.region(b.x, b.y)] += 1
            populations[sp.name] = {'occupancy': occupancy,
                'agent_ids': [b.agent_id for b in sorted(sp.members, key=lambda b:b.agent_id)],
                'reserves_pg': [b.food_level for b in sorted(sp.members, key=lambda b:b.agent_id)],
                'cumulative_uptake_pg': [getattr(b, 'experiment_uptake_pg', 0.0) for b in sorted(sp.members, key=lambda b:b.agent_id)]}
        return {'food_half_pg': halves, 'food_bins_4x4_pg': bins,
                'populations': populations, 'residence_cell_seconds': dict(self.residence)}


def spatial_snapshot(sim):
    """Common-time spatial readout; source orientation comes from the experiment."""
    if sim.experiment is None:
        raise ValueError('No experiment enabled')
    return sim.experiment.exposure()
