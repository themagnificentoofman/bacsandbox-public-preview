import random
import unittest
from unittest.mock import patch

from bacteria.core_sim import Bacterium, BacteriaSpecies, FoodNode
from bacteria.experiment_support import validate_lifecycle


def species(name='focal', **kw):
    parameters = dict(name=name, init=1, rep=1, aging=0, so=0, sf=1,
                      speed=0, ph=7, color=(0, 0, 0), food_capacity=100,
                      metabolic_rate=.1, consumption_speed=10,
                      daughter_size_policy='fixed', t6ss_activity=0,
                      bacteriocin_production=0, acetate_secretion=0, acetate_uptake=0,
                      persister_switch=0, hgt_competence=0, predation_thresh=0, share_thresh=1)
    parameters.update(kw)
    return BacteriaSpecies(**parameters)


class LifecycleQualificationTests(unittest.TestCase):
    def test_original_capacity_fails_and_amendment_passes(self):
        with self.assertRaisesRegex(ValueError, 'unreachable'):
            validate_lifecycle(species(food_capacity=10), .1, 2.5)
        result = validate_lifecycle(species(), .1, 2.5)
        self.assertEqual(result['initial_food_pg'], 2.5)
        self.assertEqual(result['capacity_pg'], 100)

    def test_fed_isolated_cell_really_divides(self):
        random.seed(314)
        sp = species()
        b = Bacterium(sp, 0, 0, initial_food=2.5)
        food = FoodNode(0, 0, 1000)
        child = None
        for _ in range(200):
            b.x = b.y = 0
            b.step([food], (0, 0, 1000), .1)
            child = b.attempt_fission(1, 1, sp.rep_thresh, .1)
            if child is not None:
                break
        self.assertIsNotNone(child)
        self.assertGreater(child.food_level, 0)

    def test_variable_daughter_capacity_overflow_returns_to_parent(self):
        sp = species(daughter_size_policy='variable')
        b = Bacterium(sp, 0, 0, initial_food=90)
        with patch('bacteria.core_sim.random.gauss', return_value=.2):
            child = b.attempt_fission(1, 1, .8, .1)
        self.assertIsNotNone(child)
        self.assertLessEqual(child.food_level, child.food_capacity)
        self.assertLessEqual(b.food_level, b.food_capacity)
        self.assertAlmostEqual(b.food_level + child.food_level, 50)

    def test_maintenance_aging_is_independent_of_mortality_aging(self):
        for enabled, expected in [(True, .1 * (1 + 601 / 180)), (False, .1)]:
            b = Bacterium(species(maintenance_aging=enabled), 0, 0, initial_food=10)
            b.age, b.lifespan_variance = 600, 1
            b.step([], (0, 0, 1000), 1)
            self.assertAlmostEqual(10 - b.food_level, expected)


if __name__ == '__main__':
    unittest.main()


from bacteria.core_sim import Simulation
from bacteria.experiment_support import ExperimentIntegrityError


def simulation(populations, seed=123, config=None, cap=5000):
    sim = Simulation(populations, {'o2': .21, 'co2': 0, 'ph': 7, 'temperature': 37,
                                  'drop_frequency': 0, 'phage_inflow': 0}, cap,
                     seed=seed, physics_backend='python', ecm_backend='python')
    sim._ensure_petri_geometry(800, 600)
    sim._populated = True
    sim.running = True
    for i, sp in enumerate(populations):
        sp.share_thresh = 1
        sp.members = [Bacterium(sp, *sim.petri_center, initial_food=2.5)]
        sp.members[0].agent_id = i
    sim.enable_experiment(config or {'washout_rate': 0, 'maintenance_aging': False}, 678)
    return sim


class ExperimentalAccountingTests(unittest.TestCase):
    def test_no_feed_cell_starves_and_population_reconciles(self):
        sim = simulation([species()])
        sim.species[0].members[0].food_level = .01
        sim.experiment.last_reserves = sim.experiment.reserves()
        sim.update(800, 600, 1)
        self.assertEqual(len(sim.species[0].members), 0)
        tick = sim.experiment.last_tick
        self.assertEqual(tick['populations']['focal']['deaths'], {'starvation': 1})
        self.assertAlmostEqual(tick['mass_residual_pg'], 0)
        self.assertEqual(sim.species[0].extinction_time, 1)

    def test_uptake_then_cull_remains_in_ledger(self):
        sim = simulation([species()], cap=0)
        sim.add_external_food(*sim.petri_center, 10)
        sim.update(800, 600, .1)
        self.assertEqual(sim.stats['total'][-1], 0)
        self.assertGreater(sim.experiment.cumulative['uptake_pg'], 0)
        self.assertAlmostEqual(sim.stats['consumed_calorie'][-1], sim.experiment.cumulative['uptake_pg'])
        self.assertEqual(sim.experiment.last_tick['populations']['focal']['deaths'], {'capacity': 1})
        self.assertAlmostEqual(sim.experiment.last_tick['mass_residual_pg'], 0)

    def test_node_limit_never_silently_deletes_food(self):
        sim = simulation([species()])
        for _ in range(1001):
            sim.add_external_food(*sim.petri_center, 1)
        sim.update(800, 600, .001)
        self.assertEqual(len(sim.food_nodes), 1001)
        self.assertAlmostEqual(sim.experiment.last_tick['mass_residual_pg'], 0, places=8)
        limited = simulation([species()], config={'max_food_nodes': 2, 'washout_rate': 0})
        for _ in range(3):
            limited.add_external_food(*limited.petri_center, 1)
        with self.assertRaisesRegex(ExperimentIntegrityError, 'resource limit'):
            limited.update(800, 600, .1)
        self.assertEqual(len(limited.food_nodes), 3)

    def test_container_order_does_not_change_dynamics(self):
        import copy
        def run(reverse):
            sim = simulation([species('F'), species('R')])
            if reverse:
                sim.species.reverse()
            sim.add_external_food(*sim.petri_center, 1)
            for _ in range(5):
                sim.update(800, 600, .1)
            agents = sorted((b for s in sim.species for b in s.members), key=lambda b:b.agent_id)
            return [(b.agent_id, b.x, b.y, b.food_level, b.age) for b in agents], sim.experiment.cumulative
        self.assertEqual(run(False), run(True))

    def test_arbitrary_neutral_labels_preserve_total_dynamics(self):
        def run(split):
            sim = simulation([species('A'), species('B')])
            if not split:
                a, b = sim.species
                for cell in b.members:
                    cell.species = a
                a.members.extend(b.members)
                b.members = []
                sim.species = [a]
                a.species_list = [a]
                sim.experiment.last_counts = sim.experiment.counts()
                sim.experiment.last_reserves = sim.experiment.reserves()
            sim.add_external_food(*sim.petri_center, 1)
            for _ in range(5):
                sim.update(800, 600, .1)
            return sorted((b.agent_id,b.x,b.y,b.food_level) for sp in sim.species for b in sp.members)
        self.assertEqual(run(True), run(False))

    def test_scheduler_gives_both_populations_first_access(self):
        winners = set()
        for seed in range(20):
            sim = simulation([species('F', metabolic_rate=0), species('R', metabolic_rate=0)])
            sim.experiment.scheduler.seed(seed)
            sim.add_external_food(*sim.petri_center, .01)
            sim.update(800, 600, .1)
            first = next(e for e in sim.experiment.last_events if e['kind'] == 'uptake')
            winners.add(first['population'])
        self.assertEqual(winners, {'F', 'R'})

    def test_daughter_not_updated_on_birth_tick(self):
        sim = simulation([species(rep=100, baranyi_q0=100)])
        b = sim.species[0].members[0]
        b.food_level = 90
        sim.experiment.last_reserves = sim.experiment.reserves()
        sim.update(800, 600, .1)
        self.assertEqual(len(sim.species[0].members), 2)
        child = next(b for b in sim.species[0].members if b.agent_id != 0)
        self.assertEqual(child.age, 0)
        self.assertEqual(sim.experiment.last_tick['births'], 1)


class DisabledMechanismTests(unittest.TestCase):
    def test_disabled_mechanisms_remain_zero(self):
        sim = simulation([species()])
        sim.add_external_food(*sim.petri_center, 100)
        for _ in range(10):
            sim.update(800, 600, .1)
            self.assertEqual(sim.stats['t6ss_deaths'][-1], 0)
            self.assertEqual(sim.stats['hgt_transfers'][-1], 0)
            self.assertFalse(any(e['kind'] == 'death' and e['cause'] != 'starvation'
                                 for e in sim.experiment.last_events))
        self.assertEqual(sim.experiment.cumulative['deaths'], 0)


class WashoutAggregationTests(unittest.TestCase):
    def test_washout_groups_preserve_exact_tick_mass_and_region(self):
        sim = simulation([species(consumption_speed=0)],
                         config={'washout_rate': .02, 'source_angle': 0, 'maintenance_aging': False})
        cx, cy = sim.petri_center
        for substrate in (0, 1, 2):
            for sign in (-1, 1):
                for mass in (.125, .25, .5):
                    sim.add_external_food(cx + sign * 100, cy, mass, substrate)
        sim.update(800, 600, .1)
        events = [e for e in sim.experiment.last_events if e['kind'] == 'washout']
        self.assertEqual(len(events), 6)
        self.assertEqual({(e['substrate_id'],e['region']) for e in events},
                         {(s,r) for s in (0,1,2) for r in ('rich','poor')})
        for event in events:
            self.assertEqual(event['node_count'], 3)
            self.assertAlmostEqual(event['amount_pg'], .875 * .02 * .1)
        self.assertAlmostEqual(sim.experiment.cumulative['washout_pg'], 5.25 * .02 * .1)
        self.assertAlmostEqual(sim.experiment.last_tick['mass_end_pg'], 5.25 * (1 - .02 * .1))
        self.assertAlmostEqual(sim.experiment.last_tick['mass_residual_pg'], 0)


class ProtocolRateLifecycleTests(unittest.TestCase):
    def test_declared_uptake_rates_can_divide_with_nominal_lifespan(self):
        # Isolate the lifecycle: immobilized, constantly fed, L=1. The actual
        # .5/1 pg/s uptake rates, temperature, initial reserve and aging law are
        # retained. A certain division trial checks reachability, not waiting time.
        for uptake in (.5, 1.0):
            with self.subTest(consumption_speed=uptake):
                random.seed(21)
                sp = species(consumption_speed=uptake, maintenance_aging=True,
                             daughter_size_policy='variable')
                b = Bacterium(sp, 0, 0, initial_food=2.5)
                b.lifespan_variance = 1.0
                food = FoodNode(0, 0, 100000)
                child = None
                for _ in range(6000):
                    b.step([food], (0, 0, 1000), .1,
                           temp_factor=sp.temperature_factor(37))
                    child = b.attempt_fission(1, 1, .8, .1)
                    if child is not None:
                        break
                self.assertIsNotNone(child)
                self.assertLess(b.age, 600)
                self.assertLessEqual(child.food_level, child.food_capacity)

    def test_fast_maintenance_aging_can_preclude_focal_division_even_when_fed(self):
        # L=.4 is the supported lifespan lower bound. Record this sensitivity
        # honestly: algebraic storage feasibility does not guarantee that every
        # cell can reach the threshold during its finite physiological lifetime.
        random.seed(21)
        sp = species(consumption_speed=.5, maintenance_aging=True)
        b = Bacterium(sp, 0, 0, initial_food=2.5)
        b.lifespan_variance = .4
        food = FoodNode(0, 0, 100000)
        peak = b.food_level
        for _ in range(6000):
            b.step([food], (0, 0, 1000), .1,
                   temp_factor=sp.temperature_factor(37))
            peak = max(peak, b.food_level)
            self.assertIsNone(b.attempt_fission(1, 1, .8, .1))
        self.assertLess(peak, 80)
        self.assertGreater(peak, 50)


class MaintenanceBoundaryTests(unittest.TestCase):
    def test_partial_feeding_rescues_after_unmet_maintenance(self):
        # Preserve the historical sequence: maintenance clips reserves at zero;
        # uptake follows; positive reserves then prevent starvation. Unpaid
        # maintenance is a recorded deficit, not carried as reserve debt.
        sim = simulation([species(consumption_speed=.5)])
        b = sim.species[0].members[0]
        b.food_level = 0
        sim.experiment.last_reserves = sim.experiment.reserves()
        for _ in range(2):
            sim.add_external_food(*sim.petri_center, .01)
            sim.update(800, 600, 1)
            self.assertEqual(len(sim.species[0].members), 1)
            self.assertAlmostEqual(b.food_level, .01)
            self.assertGreater(sim.experiment.last_tick['unmet_maintenance_pg'], .08)
        ledger = sim.experiment.cumulative
        self.assertAlmostEqual(ledger['maintenance_demand_pg'],
                               ledger['maintenance_pg'] + ledger['unmet_maintenance_pg'])
        self.assertAlmostEqual(ledger['maintenance_pg'], .01)
        self.assertEqual(ledger['deaths'], 0)
        event = next(e for e in sim.experiment.last_events if e['kind'] == 'maintenance')
        self.assertAlmostEqual(event['demand_pg'], event['amount_pg'] + event['unmet_maintenance_pg'])

    def test_alive_weighted_field_contrast_stops_after_extinction(self):
        sim = simulation([species(consumption_speed=0)])
        b = sim.species[0].members[0]
        b.food_level = 0
        sim.experiment.last_reserves = sim.experiment.reserves()
        cx, cy = sim.petri_center
        sim.add_external_food(cx + 100, cy, 10)
        sim.update(800, 600, 1)
        self.assertEqual(len(sim.species[0].members), 0)
        first = sim.experiment.last_tick['alive_weighted_exposure_contrast']
        sim.add_external_food(cx - 100, cy, 100)
        sim.update(800, 600, 1)
        tick = sim.experiment.last_tick
        self.assertEqual(tick['alive_cell_seconds'], 1)
        self.assertEqual(tick['alive_weighted_exposure_contrast'], first)
        self.assertLess(tick['exposure_contrast'], first)
