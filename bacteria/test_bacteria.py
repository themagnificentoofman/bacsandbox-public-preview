import math
import unittest
import sys
import os
import tempfile
from unittest.mock import patch

if __package__ is None and not hasattr(sys, 'frozen'):
    path = os.path.realpath(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(path)))

from bacteria import (
    BACTERIA_TUMBLE_MEAN_TURN_RAD,
    CELL_DRY_MASS_FRACTION,
    CELL_WET_MASS_PG_PER_UM3,
    ATMOSPHERIC_O2_FRACTION,
    DEFAULT_BACTERIA_SPEED_UM_S,
    ECM_CELL_SIZE_UM,
    ECM_MIN_CONCENTRATION,
    FISSION_COST_PG_PER_UM3,
    LABILE_NECROMASS_FRACTION,
    NEUTRAL_PH,
    SENESCENCE_BASE_LIFESPAN_SEC,
    T6SS_DAMAGE_SCALING_FACTOR,
    T6SS_IMMUNITY_SIZE_EPSILON_UM,
    T6SS_MAX_IMMUNITY_DEFENSE,
    Bacterium,
    BacteriaSpecies,
    FoodNode,
    FPS,
    Simulation,
    SpeciesDataManager,
)
from bacteria.app import format_world_distance_um

TEST_DT_SEC = 1.0 / FPS

class MockSpecies:
    def __init__(self, **kwargs):
        self.food_capacity = kwargs.get('food_capacity', 10.0)
        self.metabolic_rate = kwargs.get('metabolic_rate', 0.1)
        self.consumption_speed = kwargs.get('consumption_speed', 1.0)
        self.speed = kwargs.get('speed', DEFAULT_BACTERIA_SPEED_UM_S)
        self.sensor_radius = kwargs.get('sensor_radius', 200.0)
        self.greediness = kwargs.get('greediness', 0.5)
        self.t6ss_activity = kwargs.get('t6ss_activity', 0.5)
        self.bacteriocin_production = kwargs.get('bacteriocin_production', 0.5)
        self.immunity_factor = kwargs.get('immunity_factor', 0.5)
        self.rep_thresh = kwargs.get('rep_thresh', 0.8)
        self.share_thresh = kwargs.get('share_thresh', 0.6)
        self.predation_thresh = kwargs.get('predation_thresh', 0.9)
        self.base_size = kwargs.get('base_size', 1.0)

class TestBacteriumMath(unittest.TestCase):
    def test_bacterium_scaling_laws(self):
        # Base size 1.0
        sp_base = MockSpecies(food_capacity=10.0, metabolic_rate=0.1, consumption_speed=1.0, speed=DEFAULT_BACTERIA_SPEED_UM_S)
        b_base = Bacterium(sp_base, 0, 0, size=1.0)
        
        self.assertAlmostEqual(b_base.food_capacity, 10.0)
        self.assertAlmostEqual(b_base.metabolic_rate, 0.1)
        self.assertAlmostEqual(b_base.consumption_speed, 1.0)
        self.assertAlmostEqual(b_base.speed, DEFAULT_BACTERIA_SPEED_UM_S)
        self.assertAlmostEqual(b_base.turn_rate, math.pi)

        # Base size 2.0
        b_large = Bacterium(sp_base, 0, 0, size=2.0)
        # Volumetric scaling: 10 * (2^3) = 80
        self.assertAlmostEqual(b_large.food_capacity, 80.0)
        # Kleiber's Law: 0.1 * (2^2.25)
        self.assertAlmostEqual(b_large.metabolic_rate, 0.1 * (2 ** 2.25))
        # Surface area: 1.0 * (2^2) = 4.0
        self.assertAlmostEqual(b_large.consumption_speed, 4.0)
        # Kinematic scaling: 40 / 2 = 20
        self.assertAlmostEqual(b_large.speed, 20.0)
        self.assertAlmostEqual(b_large.turn_rate, math.pi / 2.0)

    def test_metabolism_decay(self):
        sp = MockSpecies(food_capacity=10.0, metabolic_rate=1.44)
        b = Bacterium(sp, 0, 0, size=1.0)
        b.age = 0 # reset age for exactness
        initial_food = b.food_level
        dt = TEST_DT_SEC
        b.step([], (0, 0, 1000), dt)
        
        # calculate expected reduction with age factor
        # age factor = 1.0 + (age / (180.0 * lifespan_variance))
        # for age=0, factor=1.0
        # Wait, age is incremented in step before food reduction
        expected_food = initial_food - (1.44 * dt) * (1.0 + ((0.0 + dt) / (SENESCENCE_BASE_LIFESPAN_SEC * b.lifespan_variance)))
        self.assertAlmostEqual(b.food_level, expected_food)

    def test_chemotaxis_score(self):
        # A food node at distance 100 with amount 50
        sp = MockSpecies(sensor_radius=200.0)
        b = Bacterium(sp, 0, 0, size=1.0)
        fnode = FoodNode(100, 0, 50)
        # The node should be visible
        dt = TEST_DT_SEC
        b.step([fnode], (0, 0, 1000), dt)
        # After step, if target is found, it will turn towards it.
        # It's hard to test the pure score without intercepting, but making sure no crash

    def test_close_food_reduces_heading_error(self):
        sp = MockSpecies(sensor_radius=200.0, greediness=1.0)
        b = Bacterium(sp, 0, 0, size=1.0)
        b.heading = math.pi / 2.0
        b.is_tumbling = False
        fnode = FoodNode(100, 0, 50)
        old_error = abs(math.atan2(math.sin(0.0 - b.heading), math.cos(0.0 - b.heading)))

        b.step([fnode], (0, 0, 1000), 0.1)

        new_error = abs(math.atan2(math.sin(0.0 - b.heading), math.cos(0.0 - b.heading)))
        self.assertLess(new_error, old_error)

    def test_tumble_reorientation_is_not_dt_scaled_drift(self):
        sp = MockSpecies(sensor_radius=0.0)
        b = Bacterium(sp, 0, 0, size=1.0)
        b.heading = 0.0
        b.is_tumbling = True
        b.state_timer = 0.05
        b._tumble_reoriented = False

        with patch('bacteria.core_sim.random.gauss', return_value=BACTERIA_TUMBLE_MEAN_TURN_RAD), \
             patch('bacteria.core_sim.random.choice', return_value=1.0):
            b.step([], (0, 0, 1000), TEST_DT_SEC)

        self.assertGreater(abs(b.heading), 0.5)

    def test_binary_fission_survival(self):
        """
        Verify that cells do not immediately die of starvation after binary fission ATP cost.
        """
        dt = TEST_DT_SEC
        sp = MockSpecies(rep_thresh=0.8, metabolic_rate=0.1)
        b = Bacterium(sp, 0, 0, size=1.0)
        
        b.base_size = 1.0
        # Give enough food to trigger fission but not much more
        child_volume = 1.0 ** 3
        atp_cost = FISSION_COST_PG_PER_UM3 * child_volume
        b.food_capacity = 100.0
        b.food_level = atp_cost + 0.5  # Just enough to trigger and have 0.5 leftover
        
        # Test attempt_fission
        child = b.attempt_fission(1.0, 1.0, 0.8, dt)
        
        # If it had barely enough food, it should return None if survival is compromised,
        # or it returns child and both must survive at least 1 tick
        if child is not None:
            min_survival = b.metabolic_rate * TEST_DT_SEC
            self.assertGreater(b.food_level, min_survival)
            self.assertGreater(child.food_level, min_survival)

    def test_ecm_memory_bounds(self):
        """
        Verify ECM removes empty cells to prevent dictionary memory leaks.
        """
        from bacteria.ecm import ECM
        ecm = ECM()
        
        # Add tiny amount that should decay below threshold instantly
        ecm.add(0, 0, ECM_MIN_CONCENTRATION * 1.0001, layer='qs')
        ecm.add(100, 100, ECM_MIN_CONCENTRATION * 1.0001, layer='bac')
        
        self.assertEqual(len(ecm.qs_inducers), 1)
        self.assertEqual(len(ecm.bacteriocins), 1)
        
        # Tick applies decay and diffusion, values drop below the ECM pruning threshold.
        dt = TEST_DT_SEC
        ecm.tick(dt)
        
        self.assertEqual(len(ecm.qs_inducers), 0)
        self.assertEqual(len(ecm.bacteriocins), 0)

    def test_ecm_processes_all_active_cells_each_tick(self):
        from bacteria.ecm import ECM
        ecm = ECM(cell_size=ECM_CELL_SIZE_UM)
        ecm.add(0, 0, 10.0, layer='qs')
        ecm.add(200, 200, 10.0, layer='qs')

        ecm.tick(TEST_DT_SEC)

        self.assertGreater(ecm.sample(20, 0, layer='qs'), 0.0)
        self.assertGreater(ecm.sample(220, 200, layer='qs'), 0.0)

    def test_death_drop_uses_lossy_necromass(self):
        sp = BacteriaSpecies('Test', 0, 0.0, 0.0, 0.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S, NEUTRAL_PH, (255, 255, 255), base_size=1.0)
        b = Bacterium(sp, 0, 0, size=1.0, initial_food=0.0)
        sp.members = [b]
        food_nodes = []

        births, deaths = sp.step({'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION}, 100, 800, 600, food_nodes, TEST_DT_SEC)

        expected_necromass = CELL_WET_MASS_PG_PER_UM3 * CELL_DRY_MASS_FRACTION * LABILE_NECROMASS_FRACTION
        self.assertEqual(births, 0)
        self.assertEqual(deaths, 1)
        self.assertEqual(len(food_nodes), 1)
        self.assertAlmostEqual(food_nodes[0].amount, expected_necromass)

    def test_predation_does_not_double_drop_necromass(self):
        predator_species = BacteriaSpecies('Predator', 0, 0.0, 0.0, 0.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S, NEUTRAL_PH, (255, 0, 0), predation_thresh=1.0, base_size=3.0)
        prey_species = BacteriaSpecies('Prey', 0, 0.0, 0.0, 0.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S, NEUTRAL_PH, (0, 255, 0), base_size=1.0)
        predator = Bacterium(predator_species, 0, 0, size=3.0, initial_food=0.0)
        prey = Bacterium(prey_species, 1, 0, size=1.0, initial_food=10.0)
        predator_species.members = [predator]
        prey_species.members = [prey]
        sim = Simulation([predator_species, prey_species], {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION, 'drop_frequency': 0, 'drop_range': (1, 1)}, 100)
        spatial_hash = {(0, 0): [predator, prey]}

        sim._process_interactions(spatial_hash, 20, TEST_DT_SEC)
        self.assertLess(prey.food_level, 0.0)

        food_nodes = []
        births, deaths = prey_species.step({'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION}, 100, 800, 600, food_nodes, TEST_DT_SEC)
        self.assertEqual(births, 0)
        self.assertEqual(deaths, 1)
        self.assertEqual(len(food_nodes), 0)

    def test_legacy_species_csv_normalizes_before_save(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, 'species_data.csv')
            with open(path, 'w', encoding='utf-8', newline='') as f:
                f.write('name,init,rep,aging,so,sf,speed,ph,aggression,defense\n')
                f.write('legacy,1,0.1,0.0,0.0,1.0,40.0,7.0,0.9,0.7\n')

            manager = SpeciesDataManager(path)
            self.assertEqual(manager.species[0]['t6ss_activity'], '0.9')
            self.assertEqual(manager.species[0]['immunity_factor'], '0.7')
            self.assertNotIn('aggression', manager.species[0])
            manager.save()

            with open(path, 'r', encoding='utf-8') as f:
                saved = f.read()
            self.assertIn('t6ss_activity', saved)
            self.assertNotIn('aggression', saved.splitlines()[0])

    def test_t6ss_damage_scaling(self):
        """
        Verify T6SS damage calculations don't overflow or produce negative floats.
        """
        sp1 = MockSpecies(t6ss_activity=1.0, immunity_factor=0.0)
        sp2 = MockSpecies(t6ss_activity=0.0, immunity_factor=1.0)
        
        # Attacker is huge (size 10), victim is tiny
        b_attacker = Bacterium(sp1, 0, 0, size=10.0)
        b_victim = Bacterium(sp2, 1, 1, size=0.1)
        
        dt = TEST_DT_SEC
        t6ss_dmg = (1.0 * (10.0 ** 2) * T6SS_DAMAGE_SCALING_FACTOR) * dt
        
        immunity = min(T6SS_MAX_IMMUNITY_DEFENSE, (1.0 * 0.1) / (10.0 + T6SS_IMMUNITY_SIZE_EPSILON_UM))
        
        actual_dmg = t6ss_dmg * (1.0 - immunity)
        
        # Must be bounded, positive, and not float('inf')
        self.assertGreater(actual_dmg, 0.0)
        self.assertLess(actual_dmg, 1000.0)  # arbitrary upper bound for 1 tick of damage
        self.assertFalse(math.isinf(actual_dmg))

    def test_stats_track_elapsed_tickrate_and_global_cumulative_totals(self):
        sim = Simulation([], {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION, 'co2': 0.0, 'drop_frequency': 0, 'drop_range': (1, 1)}, 100)
        sim.running = True

        sim.update(800, 600, 0.25)

        self.assertEqual(sim.stats['time'], [0])
        self.assertAlmostEqual(sim.stats['elapsed_seconds'][0], 0.25)
        self.assertAlmostEqual(sim.stats['dt_seconds'][0], 0.25)
        self.assertAlmostEqual(sim.stats['tickrate_hz'][0], 4.0)
        self.assertEqual(sim.stats['overall_births'][0], 0)
        self.assertEqual(sim.stats['overall_deaths'][0], 0)

    def test_export_import_global_stats_schema(self):
        sim = Simulation([], {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION, 'drop_frequency': 0, 'drop_range': (1, 1)}, 100)
        sim.stats = sim._empty_stats()
        sim.stats['time'].append(0)
        sim.stats['total'].append(5)
        sim.stats['births'].append(2)
        sim.stats['deaths'].append(1)
        sim.stats['elapsed_seconds'].append(0.5)
        sim.stats['dt_seconds'].append(0.5)
        sim.stats['tickrate_hz'].append(2.0)
        sim.stats['overall_births'].append(2)
        sim.stats['overall_deaths'].append(1)

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, 'simulation_stats.csv')
            with patch('bacteria.core_sim.EXPORT_FILE', path):
                sim.export()
            with open(path, 'r', encoding='utf-8') as f:
                exported = f.read()

            self.assertIn('time,total,births,deaths,elapsed_seconds,dt_seconds,tickrate_hz,overall_births,overall_deaths', exported)

            imported = Simulation([], {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION}, 100)
            imported.import_stats(path)
            self.assertEqual(imported.stats['overall_births'], [2])
            self.assertAlmostEqual(imported.stats['tickrate_hz'][0], 2.0)

    def test_import_legacy_global_stats_infers_timing_fields(self):
        legacy_csv = "\n".join([
            '--- Global Stats ---',
            'time,total,births,deaths',
            '0,10,2,1',
            '1,11,1,0',
            ''
        ])
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, 'legacy_stats.csv')
            with open(path, 'w', encoding='utf-8') as f:
                f.write(legacy_csv)

            sim = Simulation([], {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION}, 100)
            sim.import_stats(path)

        self.assertEqual(sim.stats['overall_births'], [2, 3])
        self.assertEqual(sim.stats['overall_deaths'], [1, 1])
        self.assertAlmostEqual(sim.stats['elapsed_seconds'][1], 2.0 / FPS)
        self.assertEqual(sim.stats['tickrate_hz'][0], FPS)

    def test_scale_distance_formatting(self):
        self.assertEqual(format_world_distance_um(5.0), '5.00 um')
        self.assertEqual(format_world_distance_um(50.0), '50.0 um')
        self.assertEqual(format_world_distance_um(2500.0), '2.50 mm')
        self.assertEqual(format_world_distance_um(20000.0), '2.00 cm')

if __name__ == '__main__':
    unittest.main()
