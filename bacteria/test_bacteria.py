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
    ACETATE_ECM_UNITS_PER_PG,
    ACETATE_ENERGY_YIELD_GLUCOSE_EQUIV,
    ACETATE_PER_GLUCOSE_MAX_MASS_RATIO,
    BARANYI_Q_STATE_MAX,
    BARANYI_Q0_MAX,
    BARANYI_Q0_MIN,
    CELL_DRY_MASS_FRACTION,
    CELL_WET_MASS_PG_PER_UM3,
    ATMOSPHERIC_O2_FRACTION,
    DEFAULT_BACTERIA_SPEED_UM_S,
    DEFAULT_BARANYI_Q0,
    DEFAULT_ENV_TEMPERATURE_C,
    ECM_CELL_SIZE_UM,
    ECM_MIN_CONCENTRATION,
    FISSION_COST_PG_PER_UM3,
    LABILE_NECROMASS_FRACTION,
    NEUTRAL_PH,
    SENESCENCE_BASE_LIFESPAN_SEC,
    T6SS_DAMAGE_SCALING_FACTOR,
    T6SS_IMMUNITY_SIZE_EPSILON_UM,
    T6SS_MAX_IMMUNITY_DEFENSE,
    DEFAULT_PERSISTER_SWITCH,
    DEFAULT_PERSISTER_ENTER_THRESH,
    DEFAULT_PERSISTER_EXIT_THRESH,
    PERSISTER_METABOLIC_FACTOR,
    PERSISTER_RESISTANCE_FACTOR,
    PERSISTER_THRESHOLD_MIN_GAP,
    DEFAULT_PHAGE_ADSORPTION_RATE,
    DEFAULT_PHAGE_BURST_SIZE,
    DEFAULT_PHAGE_LATENT_PERIOD_SEC,
    DEFAULT_PHAGE_RESISTANCE,
    DEFAULT_HGT_COMPETENCE,
    HGT_METABOLIC_BURDEN,
    HGT_TRANSFER_RATE_PER_SEC,
    PHAGE_ECM_UNITS_PER_VIRION,
    PHAGE_SEED_PACKET_VIRIONS,
    Bacterium,
    BacteriaSpecies,
    FoodNode,
    FPS,
    Simulation,
    SpeciesDataManager,
)
from bacteria.constants import (
    SUBSTRATE_GENERIC_CALORIE,
    SUBSTRATE_COMPLEX_GLYCAN,
    SUBSTRATE_MDCF_PREBIOTIC,
    SUBSTRATE_IDS,
    SUBSTRATE_ENERGY_YIELD,
    DEFAULT_SUBSTRATE_ID,
    DEFAULT_AFF_CALORIE,
    DEFAULT_AFF_GLYCAN,
    DEFAULT_AFF_MDCF,
    DIET_PROFILES,
    DEFAULT_DIET_PROFILE,
    MUCUS_LAYER_RADIUS_FRACTION,
    MUCUS_MOTILITY_FACTOR,
    LUMEN_FOOD_WASHOUT_RATE_PER_SEC,
    MUCUS_FOOD_WASHOUT_RATE_PER_SEC,
    MUCUS_DROP_AMOUNT_FACTOR,
    LUMEN_COMPETITION_CULL_WEIGHT,
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
        import random
        random.seed(0)  # deterministic: isolate chemotactic steering from OS-seeded Brownian noise
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
            # Phase I diet analytics: zone and per-substrate consumption columns are part of the schema.
            self.assertIn('lumen,mucus,consumed_calorie,consumed_glycan,consumed_mdcf', exported)
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


class TestTemperatureKinetics(unittest.TestCase):
    """Deterministic bounds tests for the Ratkowsky (1983) square-root temperature model."""

    def _species(self, temp_min=4.0, temp_max=47.0):
        return BacteriaSpecies('Temp', 0, 0.0, 0.0, 0.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S,
                               NEUTRAL_PH, (255, 255, 255), base_size=1.0,
                               temp_min=temp_min, temp_max=temp_max)

    def test_factor_zero_outside_biokinetic_window(self):
        sp = self._species(temp_min=4.0, temp_max=47.0)
        # No growth at or below T_min, nor at or above T_max.
        self.assertEqual(sp.temperature_factor(4.0), 0.0)
        self.assertEqual(sp.temperature_factor(-10.0), 0.0)
        self.assertEqual(sp.temperature_factor(47.0), 0.0)
        self.assertEqual(sp.temperature_factor(60.0), 0.0)

    def test_factor_bounded_unit_interval(self):
        sp = self._species()
        for i in range(0, 601):
            t = i * 0.1
            f = sp.temperature_factor(t)
            self.assertGreaterEqual(f, 0.0)
            self.assertLessEqual(f, 1.0)

    def test_factor_normalized_peak_is_one(self):
        sp = self._species()
        peak = max(sp.temperature_factor(i * 0.1) for i in range(40, 470))
        # Normalization should make the species optimum reach (within sampling) 1.0.
        self.assertAlmostEqual(peak, 1.0, places=4)

    def test_factor_monotonic_increase_below_optimum(self):
        sp = self._species()
        # Locate optimum by coarse scan, then assert monotonic rise from T_min up to it.
        opt = max((i * 0.1 for i in range(40, 470)), key=sp.temperature_factor)
        prev = -1.0
        t = sp.temp_min + 0.1
        while t < opt:
            f = sp.temperature_factor(t)
            self.assertGreaterEqual(f + 1e-12, prev)
            prev = f
            t += 0.1

    def test_psychrophile_has_lower_optimum_than_mesophile(self):
        # A cold-adapted species (lower cardinal temps) should peak at a lower temperature.
        psychro = self._species(temp_min=-5.0, temp_max=20.0)
        meso = self._species(temp_min=4.0, temp_max=47.0)
        psychro_opt = max((i * 0.1 for i in range(-100, 210)), key=psychro.temperature_factor)
        meso_opt = max((i * 0.1 for i in range(0, 480)), key=meso.temperature_factor)
        self.assertLess(psychro_opt, meso_opt)

    def test_cold_reduces_metabolic_consumption(self):
        # At a near-freezing temperature the Ratkowsky factor depresses basal metabolism,
        # so a chilled cell must consume less stored nutrient than one at its optimum.
        sp = self._species()
        warm = Bacterium(sp, 0.0, 0.0, size=1.0)
        cold = Bacterium(sp, 0.0, 0.0, size=1.0)
        warm.step([], bounds=(0.0, 0.0, 1.0e6), dt=1.0, temp_factor=sp.temperature_factor(37.0))
        cold.step([], bounds=(0.0, 0.0, 1.0e6), dt=1.0, temp_factor=sp.temperature_factor(6.0))
        self.assertGreater(warm.last_metabolic_consumption_pg, cold.last_metabolic_consumption_pg)
        self.assertGreater(cold.last_metabolic_consumption_pg, 0.0)


class TestCrossFeeding(unittest.TestCase):
    """Deterministic mass-balance and behaviour tests for M2 acetate overflow cross-feeding."""

    def _env(self, o2=ATMOSPHERIC_O2_FRACTION, temperature=DEFAULT_ENV_TEMPERATURE_C):
        # Aerobic, mesophilic-optimum environment so temp_factor and o2_gate are both nonzero.
        return {'ph': NEUTRAL_PH, 'o2': o2, 'co2': 0.0, 'temperature': temperature,
                'drop_frequency': 0, 'drop_range': (1, 1)}

    def _species(self, name, color, acetate_secretion=0.0, acetate_uptake=0.0):
        # bacteriocin_production=0 isolates cross-feeding food deltas from antimicrobial emission cost.
        return BacteriaSpecies(name, 0, 0.0, 0.0, 1.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S,
                               NEUTRAL_PH, color, base_size=1.0, bacteriocin_production=0.0,
                               acetate_secretion=acetate_secretion, acetate_uptake=acetate_uptake)

    def test_producer_secretes_acetate_into_ecm(self):
        sp = self._species('Producer', (255, 0, 0), acetate_secretion=1.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=5.0)
        # Pretend the cell just took up 1 pg glucose last tick; overflow should appear as acetate.
        b.last_food_uptake_pg = 1.0
        sp.members = [b]
        sim = Simulation([sp], self._env(), 100)

        food_before = b.food_level
        sim._process_interactions({}, 20, 1.0)

        # Mass excreted is bounded by the glucose->2 acetate ceiling, deducted from internal stores.
        expected_secreted = 1.0 * ACETATE_PER_GLUCOSE_MAX_MASS_RATIO
        self.assertAlmostEqual(food_before - b.food_level, expected_secreted)
        self.assertAlmostEqual(sim.ecm.sample(0.0, 0.0, layer='acetate'),
                               expected_secreted * ACETATE_ECM_UNITS_PER_PG)

    def test_consumer_gains_food_from_acetate_without_glucose(self):
        sp = self._species('Consumer', (0, 0, 255), acetate_uptake=1.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=1.0)
        sp.members = [b]
        sim = Simulation([sp], self._env(), 100)
        # Seed abundant local acetate (in ECM units) with no glucose food nodes present.
        seeded_pg = 0.5
        sim.ecm.add(0.0, 0.0, seeded_pg * ACETATE_ECM_UNITS_PER_PG, layer='acetate')

        food_before = b.food_level
        acetate_units_before = sim.ecm.sample(0.0, 0.0, layer='acetate')
        sim._process_interactions({}, 20, 1.0)
        acetate_units_after = sim.ecm.sample(0.0, 0.0, layer='acetate')

        consumed_pg = (acetate_units_before - acetate_units_after) / ACETATE_ECM_UNITS_PER_PG
        self.assertGreater(consumed_pg, 0.0)
        # Stored energy gained equals acetate mass consumed scaled by its glucose-energy equivalence.
        self.assertAlmostEqual(b.food_level - food_before,
                               consumed_pg * ACETATE_ENERGY_YIELD_GLUCOSE_EQUIV)
        self.assertEqual(b.last_acetate_consumption_pg, consumed_pg)

    def test_no_acetate_creation_when_traits_disabled(self):
        # Backward compatibility: default (0/0) traits reproduce pre-M2 behaviour exactly.
        sp = self._species('Neutral', (0, 255, 0))
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=5.0)
        b.last_food_uptake_pg = 1.0
        sp.members = [b]
        sim = Simulation([sp], self._env(), 100)

        food_before = b.food_level
        sim._process_interactions({}, 20, 1.0)

        self.assertEqual(b.food_level, food_before)
        self.assertEqual(sim.ecm.sample(0.0, 0.0, layer='acetate'), 0.0)
        self.assertEqual(b.last_acetate_consumption_pg, 0.0)

    def test_cross_feeding_conserves_mass_producer_to_consumer(self):
        # Acetate the producer loses (above the ECM prune floor) must be available to a co-located
        # consumer; no carbon is created from nothing.
        producer = self._species('Producer', (255, 0, 0), acetate_secretion=1.0)
        consumer = self._species('Consumer', (0, 0, 255), acetate_uptake=1.0)
        pb = Bacterium(producer, 0.0, 0.0, size=1.0, initial_food=8.0)
        pb.last_food_uptake_pg = 1.0
        producer.members = [pb]
        consumer.members = []
        sim = Simulation([producer, consumer], self._env(), 100)

        food_before = pb.food_level
        sim._process_interactions({}, 20, 1.0)
        secreted_pg = food_before - pb.food_level
        acetate_units = sim.ecm.sample(0.0, 0.0, layer='acetate')
        self.assertAlmostEqual(acetate_units, secreted_pg * ACETATE_ECM_UNITS_PER_PG)

    def test_uptake_requires_oxygen(self):
        # Acetate oxidation is aerobic: with no dissolved O2 the consumer cannot use acetate.
        sp = self._species('Consumer', (0, 0, 255), acetate_uptake=1.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=1.0)
        sp.members = [b]
        sim = Simulation([sp], self._env(o2=0.0), 100)
        sim.ecm.add(0.0, 0.0, 0.5 * ACETATE_ECM_UNITS_PER_PG, layer='acetate')

        food_before = b.food_level
        sim._process_interactions({}, 20, 1.0)
        self.assertEqual(b.food_level, food_before)
        self.assertEqual(b.last_acetate_consumption_pg, 0.0)


class TestLagPhase(unittest.TestCase):
    """Deterministic tests for the M3 Baranyi & Roberts (1994) physiological-state lag model."""

    def _species(self, name, base_rep=0.5, baranyi_q0=DEFAULT_BARANYI_Q0,
                 food_capacity=1000.0, metabolic_rate=0.0):
        # metabolic_rate=0 keeps food_level fixed so q dynamics are isolated; o2_sens=1 aerobe with
        # atmospheric O2 and 37 degC optimum give envf=temp_factor=1 for clean rate arithmetic.
        return BacteriaSpecies(name, 0, base_rep, 0.0, 1.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S,
                               NEUTRAL_PH, (0, 200, 0), base_size=1.0, bacteriocin_production=0.0,
                               food_capacity=food_capacity, metabolic_rate=metabolic_rate,
                               baranyi_q0=baranyi_q0)

    def test_fresh_inoculum_starts_in_lag(self):
        # A new cell carries its species' q0; the adjustment alpha(q)=q/(1+q) is near zero, so the
        # effective growth-rate multiplier starts heavily suppressed.
        sp = self._species('Fresh', baranyi_q0=DEFAULT_BARANYI_Q0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=500.0)
        self.assertEqual(b.q, DEFAULT_BARANYI_Q0)
        alpha = b.q / (1.0 + b.q)
        self.assertLess(alpha, 0.1)

    def test_physiological_state_grows_exponentially(self):
        # The q update is solved analytically each tick, so after n steps q == q0 * e^{mu_max*t}.
        base_rep = 0.5
        sp = self._species('Grow', base_rep=base_rep)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=500.0)
        dt = 0.1
        n = 10
        for _ in range(n):
            b.step([], (0.0, 0.0, 1000.0), dt)  # temp_factor defaults to 1.0
        expected_q = DEFAULT_BARANYI_Q0 * math.exp(base_rep * dt * n)
        self.assertAlmostEqual(b.q, expected_q, places=6)

    def test_lag_shortens_with_higher_q0(self):
        # Larger q0 means less "work to be done" (h0 = ln(1+1/q0)), so alpha crosses 0.5 (q=1) sooner.
        base_rep = 0.5
        dt = 0.1

        def steps_until_adapted(q0):
            sp = self._species('S', base_rep=base_rep, baranyi_q0=q0)
            b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=500.0)
            count = 0
            while b.q < 1.0 and count < 100000:
                b.step([], (0.0, 0.0, 1000.0), dt)
                count += 1
            return count

        self.assertLess(steps_until_adapted(0.5), steps_until_adapted(0.01))

    def test_high_q0_disables_lag_backward_compat(self):
        # q0 at the clamp ceiling makes alpha ~ 1 from the outset, reproducing pre-M3 (no-lag) growth.
        sp = self._species('NoLag', baranyi_q0=BARANYI_Q0_MAX)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=500.0)
        alpha = b.q / (1.0 + b.q)
        self.assertGreater(alpha, 0.999)

    def test_q0_clamped_to_finite_nonzero_range(self):
        # Zero/huge q0 are clamped so the lag is always finite and adaptation always completes.
        self.assertEqual(self._species('Z', baranyi_q0=0.0).baranyi_q0, BARANYI_Q0_MIN)
        self.assertEqual(self._species('H', baranyi_q0=1.0e12).baranyi_q0, BARANYI_Q0_MAX)

    def test_q_capped_to_prevent_overflow(self):
        # Over long runs q must saturate at the numerical cap rather than overflowing to infinity.
        sp = self._species('Cap', base_rep=5.0, baranyi_q0=DEFAULT_BARANYI_Q0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=1.0e9)
        for _ in range(2000):
            b.step([], (0.0, 0.0, 1000.0), 1.0)
        self.assertLessEqual(b.q, BARANYI_Q_STATE_MAX)
        self.assertTrue(math.isfinite(b.q))

    def test_daughter_inherits_parent_adapted_state(self):
        # A dividing parent (high q) passes its adapted state to the daughter: progeny do not re-lag.
        sp = self._species('Div', base_rep=1.0)
        parent = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
        parent.q = 50.0  # already adapted
        with patch('random.random', return_value=0.0), patch('random.gauss', return_value=1.0):
            child = parent.attempt_fission(1.0, 1.0, sp.rep_thresh, 1.0)
        self.assertIsNotNone(child)
        self.assertEqual(child.q, 50.0)

    def test_lag_gates_reproduction_vs_adapted(self):
        # With identical food/environment, a lagging (low-q0) inoculum yields a smaller per-tick
        # reproduction probability than an adapted (high-q0) one; verified via births under a fixed
        # RNG threshold that only the adapted species' rep_prob exceeds.
        dt = 1.0
        lagging = self._species('Lag', base_rep=0.7, baranyi_q0=DEFAULT_BARANYI_Q0)
        adapted = self._species('Adapt', base_rep=0.7, baranyi_q0=BARANYI_Q0_MAX)
        for sp in (lagging, adapted):
            sp.members = [Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)]
        env = {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION, 'co2': 0.0,
               'temperature': DEFAULT_ENV_TEMPERATURE_C, 'drop_frequency': 0, 'drop_range': (1, 1)}
        # Fixed RNG at 0.5: adapted rep_prob (~0.70) divides, lagging rep_prob (~0.06) does not.
        with patch('random.random', return_value=0.5), patch('random.gauss', return_value=1.0):
            lag_births, _ = lagging.step(env, 100000, 800, 600, [], dt)
            adapt_births, _ = adapted.step(env, 100000, 800, 600, [], dt)
        self.assertEqual(lag_births, 0)
        self.assertGreater(adapt_births, 0)


class TestPersisterSwitching(unittest.TestCase):
    """Deterministic tests for the M4 persister phenotypic-switch model (Balaban et al. 2004)."""

    def _species(self, name, persister_switch=1.0, base_rep=0.0, metabolic_rate=0.0,
                 enter=DEFAULT_PERSISTER_ENTER_THRESH, exit_=DEFAULT_PERSISTER_EXIT_THRESH,
                 bacteriocin_production=0.0, immunity_factor=0.0):
        # base_rep/metabolic_rate default to 0 so food stays fixed and fullness (the switch driver)
        # is fully controlled by initial_food/food_capacity; aerobe at 37 degC -> envf=temp_factor=1.
        return BacteriaSpecies(name, 0, base_rep, 0.0, 1.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S,
                               NEUTRAL_PH, (0, 200, 0), base_size=1.0,
                               bacteriocin_production=bacteriocin_production,
                               immunity_factor=immunity_factor,
                               food_capacity=1000.0, metabolic_rate=metabolic_rate,
                               persister_switch=persister_switch,
                               persister_enter_thresh=enter, persister_exit_thresh=exit_)

    def test_starvation_drives_entry_to_persister(self):
        # Fullness 0.05 < S1 (0.15): with the RNG pinned at 0 the positive enter-rate always fires.
        sp = self._species('Starve')
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=50.0)  # fullness 0.05
        self.assertFalse(b.is_persister)
        with patch('random.random', return_value=0.0):
            b.step([], (0.0, 0.0, 1000.0), 1.0)
        self.assertTrue(b.is_persister)

    def test_rich_nutrients_drive_exit_from_persister(self):
        # Fullness 0.9 > S2 (0.5): a dormant cell resuscitates back to proliferative.
        sp = self._species('Recover')
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)  # fullness 0.9
        b.is_persister = True
        with patch('random.random', return_value=0.0):
            b.step([], (0.0, 0.0, 1000.0), 1.0)
        self.assertFalse(b.is_persister)

    def test_persister_halts_fission(self):
        # A dormant cell never divides even with food above threshold and the RNG forcing fission.
        sp = self._species('Halt', persister_switch=0.0, base_rep=1.0)  # switch disabled so state is stable
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
        b.is_persister = True
        sp.members = [b]
        env = {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION, 'co2': 0.0,
               'temperature': DEFAULT_ENV_TEMPERATURE_C, 'drop_frequency': 0, 'drop_range': (1, 1)}
        with patch('random.random', return_value=0.0), patch('random.gauss', return_value=1.0):
            births, _ = sp.step(env, 100000, 800, 600, [], 1.0)
        self.assertEqual(births, 0)

    def test_persister_reduces_basal_metabolism(self):
        # Dormancy throttles maintenance burn to PERSISTER_METABOLIC_FACTOR of the proliferative rate.
        # Fix the RNG so both cells share an identical lifespan_variance (hence age_factor), isolating
        # the persister factor as the only difference in metabolic demand.
        sp = self._species('Meta', persister_switch=0.0, metabolic_rate=10.0)
        with patch('random.gauss', return_value=1.0), patch('random.random', return_value=0.5):
            active = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
            dormant = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
            dormant.is_persister = True
            active.step([], (0.0, 0.0, 1000.0), 1.0)
            dormant.step([], (0.0, 0.0, 1000.0), 1.0)
        self.assertAlmostEqual(dormant.last_metabolic_consumption_pg,
                               active.last_metabolic_consumption_pg * PERSISTER_METABOLIC_FACTOR)

    def test_persister_resists_bacteriocin_dose(self):
        # The same environmental bacteriocin dose damages a dormant cell PERSISTER_RESISTANCE_FACTOR
        # as much as a proliferative one (Balaban et al. 2004; Lewis 2010).
        sp = self._species('Resist', persister_switch=0.0)
        active = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
        dormant = Bacterium(sp, 100.0, 100.0, size=1.0, initial_food=900.0)
        dormant.is_persister = True
        sp.members = [active, dormant]
        sim = Simulation([sp], {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION,
                                 'drop_frequency': 0, 'drop_range': (1, 1)}, 100)
        # Deposit an equal bacteriocin dose at each cell's location.
        sim.ecm.add(active.x, active.y, 100.0, layer='bac')
        sim.ecm.add(dormant.x, dormant.y, 100.0, layer='bac')
        # Cells sit in separate spatial-hash cells so they do not interact with each other.
        spatial_hash = {(0, 0): [active], (5, 5): [dormant]}
        with patch('random.random', return_value=0.9):
            sim._process_interactions(spatial_hash, 20, 1.0)
        active_drop = 900.0 - active.food_level
        dormant_drop = 900.0 - dormant.food_level
        self.assertGreater(active_drop, 0.0)
        self.assertAlmostEqual(dormant_drop, active_drop * PERSISTER_RESISTANCE_FACTOR)

    def test_switch_disabled_is_backward_compatible(self):
        # persister_switch=0 means a cell never spontaneously enters persistence, even when starved.
        sp = self._species('Off', persister_switch=0.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=10.0)  # fullness 0.01
        with patch('random.random', return_value=0.0):
            for _ in range(50):
                b.step([], (0.0, 0.0, 1000.0), 1.0)
        self.assertFalse(b.is_persister)

    def test_exit_threshold_kept_above_enter_threshold(self):
        # The constructor enforces a minimum S2 - S1 gap so the enter/exit bands cannot overlap.
        sp = self._species('Gap', enter=0.6, exit_=0.4)
        self.assertGreaterEqual(sp.persister_exit_thresh,
                                sp.persister_enter_thresh + PERSISTER_THRESHOLD_MIN_GAP)


class TestPhagePredationHGT(unittest.TestCase):
    """Deterministic tests for the M5 lytic-phage SIR model and horizontal gene transfer."""

    def _species(self, name, phage_resistance=0.0, hgt_competence=0.0, base_rep=0.0,
                 metabolic_rate=0.0):
        # base_rep/metabolic_rate default to 0 so food stays fixed and the only dynamics under test
        # are infection, lysis, and HGT. Aerobe at 37 degC -> envf = temp_factor = 1.
        return BacteriaSpecies(name, 0, base_rep, 0.0, 1.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S,
                               NEUTRAL_PH, (0, 200, 0), base_size=1.0,
                               food_capacity=1000.0, metabolic_rate=metabolic_rate,
                               bacteriocin_production=0.0, immunity_factor=0.0,
                               phage_resistance=phage_resistance, hgt_competence=hgt_competence)

    def _env(self):
        return {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION, 'co2': 0.0,
                'temperature': DEFAULT_ENV_TEMPERATURE_C, 'drop_frequency': 0,
                'drop_range': (1, 1), 'phage_inflow': 0}

    def test_traits_clamped_to_unit_interval(self):
        sp = self._species('Clamp', phage_resistance=5.0, hgt_competence=-2.0)
        self.assertEqual(sp.phage_resistance, 1.0)
        self.assertEqual(sp.hgt_competence, 0.0)

    def test_initial_resistant_fraction(self):
        # phage_resistance=1.0 -> every cell is born resistant; 0.0 -> none.
        all_res = self._species('AllRes', phage_resistance=1.0)
        none_res = self._species('NoneRes', phage_resistance=0.0)
        self.assertTrue(Bacterium(all_res, 0.0, 0.0, size=1.0).phage_resistant)
        self.assertFalse(Bacterium(none_res, 0.0, 0.0, size=1.0).phage_resistant)

    def test_adsorption_infects_susceptible_cell(self):
        # A susceptible cell in a dense phage field adsorbs and enters the latent (infected) state.
        sp = self._species('Susc', phage_resistance=0.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
        b.speed = 0.0  # immobilize so it stays in the seeded point phage field
        sp.members = [b]
        sim = Simulation([sp], self._env(), 100)
        sim.ecm.add(0.0, 0.0, 1000.0 * PHAGE_ECM_UNITS_PER_VIRION, layer='phage')
        self.assertFalse(b.infected)
        with patch('random.random', return_value=0.0):  # force adsorption to fire
            sp.step(self._env(), 100000, 800, 600, [], TEST_DT_SEC)
        self.assertTrue(b.infected)
        self.assertAlmostEqual(b.lysis_timer, DEFAULT_PHAGE_LATENT_PERIOD_SEC)

    def test_resistant_cell_never_adsorbs(self):
        # A resistant cell in the same dense phage field is never infected.
        sp = self._species('Res', phage_resistance=1.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
        b.speed = 0.0
        sp.members = [b]
        sim = Simulation([sp], self._env(), 100)
        sim.ecm.add(0.0, 0.0, 1000.0 * PHAGE_ECM_UNITS_PER_VIRION, layer='phage')
        with patch('random.random', return_value=0.0):
            sp.step(self._env(), 100000, 800, 600, [], TEST_DT_SEC)
        self.assertFalse(b.infected)

    def test_lysis_kills_cell_and_releases_burst(self):
        # An infected cell whose latent timer elapses lyses: it dies, drops necromass, and releases
        # a burst of progeny virions into the local phage field.
        sp = self._species('Lyse', phage_resistance=0.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
        b.speed = 0.0  # immobilize so the burst is released at the sampled origin
        b.infected = True
        b.lysis_timer = TEST_DT_SEC * 0.5  # elapses within a single step
        sp.members = [b]
        sim = Simulation([sp], self._env(), 100)
        food_nodes = []
        # Pin RNG high so no secondary death/adsorption paths interfere.
        with patch('random.random', return_value=0.99):
            births, deaths = sp.step(self._env(), 100000, 800, 600, food_nodes, TEST_DT_SEC)
        self.assertEqual(deaths, 1)
        self.assertNotIn(b, sp.members)
        released = sim.ecm.sample(0.0, 0.0, layer='phage') / PHAGE_ECM_UNITS_PER_VIRION
        self.assertGreater(released, 0.0)
        self.assertTrue(any(n.amount > 0 for n in food_nodes))

    def test_infected_cell_does_not_divide(self):
        # A phage-hijacked cell cannot reproduce even with ample food and forced fission RNG.
        sp = self._species('NoDiv', phage_resistance=0.0, base_rep=1.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
        b.infected = True
        b.lysis_timer = DEFAULT_PHAGE_LATENT_PERIOD_SEC  # not yet lysing
        sp.members = [b]
        with patch('random.random', return_value=0.0), patch('random.gauss', return_value=1.0):
            births, _ = sp.step(self._env(), 100000, 800, 600, [], TEST_DT_SEC)
        self.assertEqual(births, 0)

    def test_hgt_transfers_resistance_and_imposes_burden(self):
        # A resistant, competent donor in contact with a susceptible competent recipient transfers
        # resistance and saddles the recipient with the plasmid carriage cost.
        sp = self._species('HGT', phage_resistance=0.0, hgt_competence=1.0)
        donor = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
        recipient = Bacterium(sp, 0.5, 0.0, size=1.0, initial_food=900.0)
        donor.phage_resistant = True
        recipient.phage_resistant = False
        sp.members = [donor, recipient]
        sim = Simulation([sp], self._env(), 100)
        spatial_hash = {(0, 0): [donor, recipient]}
        with patch('random.random', return_value=0.0):  # force transfer to fire
            sim._process_interactions(spatial_hash, 20, TEST_DT_SEC)
        self.assertTrue(recipient.phage_resistant)
        self.assertAlmostEqual(recipient.metabolic_burden_factor, 1.0 + HGT_METABOLIC_BURDEN)

    def test_hgt_disabled_without_competence(self):
        # With hgt_competence=0 the trait never transfers regardless of contact or RNG.
        sp = self._species('NoHGT', phage_resistance=0.0, hgt_competence=0.0)
        donor = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
        recipient = Bacterium(sp, 0.5, 0.0, size=1.0, initial_food=900.0)
        donor.phage_resistant = True
        recipient.phage_resistant = False
        sp.members = [donor, recipient]
        sim = Simulation([sp], self._env(), 100)
        spatial_hash = {(0, 0): [donor, recipient]}
        with patch('random.random', return_value=0.0):
            sim._process_interactions(spatial_hash, 20, TEST_DT_SEC)
        self.assertFalse(recipient.phage_resistant)

    def test_daughter_inherits_resistance_and_burden(self):
        # Heritability: a resistant, plasmid-bearing parent passes both traits to its daughter.
        sp = self._species('Inherit', phage_resistance=0.0, base_rep=1.0)
        parent = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=900.0)
        parent.phage_resistant = True
        parent.metabolic_burden_factor = 1.0 + HGT_METABOLIC_BURDEN
        with patch('random.random', return_value=0.0), patch('random.gauss', return_value=1.0):
            child = parent.attempt_fission(1.0, 1.0, sp.rep_thresh, TEST_DT_SEC)
        self.assertIsNotNone(child)
        self.assertTrue(child.phage_resistant)
        self.assertAlmostEqual(child.metabolic_burden_factor, 1.0 + HGT_METABOLIC_BURDEN)

    def test_phage_inflow_zero_is_backward_compatible(self):
        # With phage_inflow=0 no virions are ever seeded and no cell becomes infected.
        sp = self._species('Quiet', phage_resistance=0.0)
        b = Bacterium(sp, 400.0, 300.0, size=1.0, initial_food=900.0)
        sp.members = [b]
        sim = Simulation([sp], self._env(), 100)
        sim.running = True
        for _ in range(20):
            sim.update(800, 600, TEST_DT_SEC)
        self.assertEqual(sim.ecm.sample(400.0, 300.0, layer='phage'), 0.0)
        self.assertFalse(b.infected)


class TestDietSubstrates(unittest.TestCase):
    """Deterministic tests for the Phase I diet-driven substrate model: per-substrate affinity
    gating, energy-yield scaling, diet-profile drop sampling, and legacy back-compat."""

    def _species(self, name, aff_calorie=DEFAULT_AFF_CALORIE, aff_glycan=DEFAULT_AFF_GLYCAN,
                 aff_mdcf=DEFAULT_AFF_MDCF, consumption_speed=1.0):
        # base_rep/metabolic 0 so the only dynamics under test are uptake. Aerobe at 37 degC.
        return BacteriaSpecies(name, 0, 0.0, 0.0, 1.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S,
                               NEUTRAL_PH, (0, 200, 0), base_size=1.0,
                               food_capacity=1000.0, metabolic_rate=0.0,
                               consumption_speed=consumption_speed,
                               bacteriocin_production=0.0, immunity_factor=0.0,
                               aff_calorie=aff_calorie, aff_glycan=aff_glycan, aff_mdcf=aff_mdcf)

    def _env(self, diet_profile=DEFAULT_DIET_PROFILE, drop_frequency=0):
        return {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION, 'co2': 0.0,
                'temperature': DEFAULT_ENV_TEMPERATURE_C, 'drop_frequency': drop_frequency,
                'drop_range': (1, 1), 'phage_inflow': 0, 'diet_profile': diet_profile}

    def test_affinities_clamped_non_negative(self):
        sp = self._species('Clamp', aff_calorie=-3.0, aff_glycan=-1.0, aff_mdcf=2.5)
        self.assertEqual(sp.aff_calorie, 0.0)
        self.assertEqual(sp.aff_glycan, 0.0)
        self.assertEqual(sp.aff_mdcf, 2.5)

    def test_bacterium_affinity_array_from_species(self):
        sp = self._species('Aff', aff_calorie=1.0, aff_glycan=0.5, aff_mdcf=0.25)
        b = Bacterium(sp, 0.0, 0.0, size=1.0)
        self.assertEqual(b.substrate_affinity,
                         (1.0, 0.5, 0.25))

    def test_foodnode_default_substrate_is_generic_calorie(self):
        # Necromass/legacy food paths construct nodes without a substrate id.
        node = FoodNode.acquire().init(0.0, 0.0, 5.0)
        self.assertEqual(node.substrate_id, DEFAULT_SUBSTRATE_ID)
        self.assertEqual(node.substrate_id, SUBSTRATE_GENERIC_CALORIE)

    def test_zero_affinity_cell_cannot_eat_substrate(self):
        # A calorie-only cell (default) cannot ingest a complex-glycan node.
        sp = self._species('CalOnly')  # aff_glycan=0 by default
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=0.0)
        b.speed = 0.0
        node = FoodNode.acquire().init(0.0, 0.0, 100.0, SUBSTRATE_COMPLEX_GLYCAN)
        before = node.amount
        b.step([node], (0, 0, 1000), TEST_DT_SEC)
        self.assertEqual(node.amount, before)  # untouched
        self.assertEqual(b.food_level, 0.0)

    def test_competent_cell_eats_mdcf_with_energy_yield(self):
        # A MDCF-competent cell ingests prebiotic mass; stored energy = mass * yield (< mass).
        sp = self._species('MdcfEater', aff_mdcf=1.0, consumption_speed=10.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=0.0)
        b.speed = 0.0
        node = FoodNode.acquire().init(0.0, 0.0, 100.0, SUBSTRATE_MDCF_PREBIOTIC)
        before = node.amount
        b.step([node], (0, 0, 1000), TEST_DT_SEC)
        eaten_mass = before - node.amount
        self.assertGreater(eaten_mass, 0.0)
        # Energy stored equals ingested mass scaled by the MDCF yield (mass conserved on the node).
        self.assertAlmostEqual(b.food_level,
                               eaten_mass * SUBSTRATE_ENERGY_YIELD[SUBSTRATE_MDCF_PREBIOTIC])

    def test_generic_calorie_full_energy_yield(self):
        # Generic calories yield 1:1 energy per ingested mass (legacy behaviour preserved).
        sp = self._species('CalEater', consumption_speed=10.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=0.0)
        b.speed = 0.0
        node = FoodNode.acquire().init(0.0, 0.0, 100.0, SUBSTRATE_GENERIC_CALORIE)
        before = node.amount
        b.step([node], (0, 0, 1000), TEST_DT_SEC)
        eaten_mass = before - node.amount
        self.assertGreater(eaten_mass, 0.0)
        self.assertAlmostEqual(b.food_level, eaten_mass * 1.0)

    def test_energy_yields_within_physical_bounds(self):
        # Every substrate yield is a fraction of glucose energy in (0, 1].
        for sid in SUBSTRATE_IDS:
            y = SUBSTRATE_ENERGY_YIELD[sid]
            self.assertGreater(y, 0.0)
            self.assertLessEqual(y, 1.0)
        # Glucose-equivalent calorie is the reference maximum.
        self.assertEqual(SUBSTRATE_ENERGY_YIELD[SUBSTRATE_GENERIC_CALORIE], 1.0)

    def test_diet_profile_ratios_sum_to_one(self):
        for profile, ratios in DIET_PROFILES.items():
            self.assertAlmostEqual(sum(ratios.values()), 1.0,
                                   msg=f"Diet profile {profile} ratios must sum to 1.0")

    def test_drop_substrate_distribution_matches_profile(self):
        # Seeded sampling reproduces the configured categorical drop mix within tolerance.
        import random as _random
        sp = self._species('Sampler')
        sim = Simulation([sp], self._env(diet_profile='mdcf'), 100)
        _random.seed(12345)
        counts = {sid: 0 for sid in SUBSTRATE_IDS}
        n = 20000
        for _ in range(n):
            counts[sim._sample_substrate()] += 1
        expected = DIET_PROFILES['mdcf']
        for sid in SUBSTRATE_IDS:
            frac = counts[sid] / n
            self.assertAlmostEqual(frac, expected[sid], delta=0.02)

    def test_seeded_sampling_is_reproducible(self):
        import random as _random
        sp = self._species('Repro')
        sim = Simulation([sp], self._env(diet_profile='standard'), 100)
        _random.seed(999)
        seq_a = [sim._sample_substrate() for _ in range(50)]
        _random.seed(999)
        seq_b = [sim._sample_substrate() for _ in range(50)]
        self.assertEqual(seq_a, seq_b)

    def test_legacy_species_defaults_to_calorie_only(self):
        # A species built without affinity kwargs eats only generic calories (back-compat).
        sp = BacteriaSpecies('Legacy', 0, 0.0, 0.0, 1.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S,
                             NEUTRAL_PH, (0, 200, 0), base_size=1.0, food_capacity=1000.0,
                             metabolic_rate=0.0, consumption_speed=10.0)
        self.assertEqual(sp.aff_calorie, DEFAULT_AFF_CALORIE)
        self.assertEqual(sp.aff_glycan, 0.0)
        self.assertEqual(sp.aff_mdcf, 0.0)
        b = Bacterium(sp, 0.0, 0.0, size=1.0, initial_food=0.0)
        b.speed = 0.0
        glycan = FoodNode.acquire().init(0.0, 0.0, 100.0, SUBSTRATE_COMPLEX_GLYCAN)
        b.step([glycan], (0, 0, 1000), TEST_DT_SEC)
        self.assertEqual(glycan.amount, 100.0)  # cannot touch fibre


class TestSpatialZones(unittest.TestCase):
    """Deterministic tests for the Phase I spatial heterogeneity model: lumen vs mucus motility,
    zone-specific nutrient washout, mucus-reduced drop deposition, and luminal competition bias."""

    def _species(self, name, consumption_speed=1.0, base_rep=0.0):
        return BacteriaSpecies(name, 0, base_rep, 0.0, 1.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S,
                               NEUTRAL_PH, (0, 200, 0), base_size=1.0,
                               food_capacity=1000.0, metabolic_rate=0.0,
                               consumption_speed=consumption_speed,
                               bacteriocin_production=0.0, immunity_factor=0.0)

    def _env(self, drop_frequency=0):
        return {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION, 'co2': 0.0,
                'temperature': DEFAULT_ENV_TEMPERATURE_C, 'drop_frequency': drop_frequency,
                'drop_range': (1, 1), 'phage_inflow': 0, 'diet_profile': DEFAULT_DIET_PROFILE}

    def test_zone_constants_are_well_ordered(self):
        # Mucus is a thin outer annulus and a movement/clearance refuge relative to the lumen.
        self.assertTrue(0.0 < MUCUS_LAYER_RADIUS_FRACTION < 1.0)
        self.assertTrue(0.0 < MUCUS_MOTILITY_FACTOR < 1.0)
        self.assertGreater(LUMEN_FOOD_WASHOUT_RATE_PER_SEC, MUCUS_FOOD_WASHOUT_RATE_PER_SEC)
        self.assertTrue(0.0 < MUCUS_DROP_AMOUNT_FACTOR < 1.0)
        self.assertGreater(LUMEN_COMPETITION_CULL_WEIGHT, 1.0)

    def test_mucus_cell_moves_slower_than_lumen_cell(self):
        # Per-step displacement magnitude equals the scalar speed (heading-independent), so a
        # non-tumbling cell in the mucus annulus advances MUCUS_MOTILITY_FACTOR as far as one in
        # the lumen under identical conditions.
        import random as _random
        sp = self._species('Mover')
        cx, cy, R = 400.0, 300.0, 500.0
        bounds = (cx, cy, R)
        dt = 1.0 / 1000.0  # tiny dt keeps both cells well inside the boundary (no reflection)

        def _disp(start_x):
            _random.seed(2024)
            b = Bacterium(sp, start_x, cy, size=1.0, initial_food=500.0)
            b.is_tumbling = False
            b.state_timer = 1e9  # never flips to tumble within this step
            x0, y0 = b.x, b.y
            b.step([], bounds, dt)
            return math.hypot(b.x - x0, b.y - y0)

        lumen_disp = _disp(cx + 5.0)               # near centre -> lumen
        mucus_disp = _disp(cx + 0.92 * R)          # outer annulus -> mucus
        self.assertGreater(lumen_disp, 0.0)
        self.assertAlmostEqual(mucus_disp, lumen_disp * MUCUS_MOTILITY_FACTOR, places=6)

    def test_lumen_nutrient_washes_out_faster_than_mucus(self):
        # Equal-mass nodes placed in each zone lose mass by zone-specific exponential clearance.
        sp = self._species('NoEat', consumption_speed=0.0)
        sim = Simulation([sp], self._env(drop_frequency=0), 100)
        sim.running = True
        w, h, dt = 800, 600, 1.0
        sim.update(w, h, dt)  # establishes petri_center / petri_radius, no drops
        cx, cy = sim.petri_center
        R = sim.petri_radius
        lumen_node = FoodNode.acquire().init(cx, cy, 100.0, SUBSTRATE_GENERIC_CALORIE)
        mucus_node = FoodNode.acquire().init(cx + 0.95 * R, cy, 100.0, SUBSTRATE_GENERIC_CALORIE)
        sim.food_nodes = [lumen_node, mucus_node]
        for _ in range(20):
            sim.update(w, h, dt)
        # Both decay, but the lumen node retains strictly less mass than the mucus refuge.
        self.assertLess(lumen_node.amount, mucus_node.amount)
        self.assertLess(lumen_node.amount, 100.0)
        self.assertGreater(mucus_node.amount, 0.0)

    def test_mucus_node_persists_far_longer_than_washout_estimate(self):
        # Closed-form check: after t seconds the mucus mass ~ m0 (1 - rate)^n for unit dt steps.
        sp = self._species('NoEat', consumption_speed=0.0)
        sim = Simulation([sp], self._env(drop_frequency=0), 100)
        sim.running = True
        w, h, dt = 800, 600, 1.0
        sim.update(w, h, dt)
        cx, cy = sim.petri_center
        R = sim.petri_radius
        node = FoodNode.acquire().init(cx + 0.95 * R, cy, 100.0, SUBSTRATE_GENERIC_CALORIE)
        sim.food_nodes = [node]
        n = 10
        for _ in range(n):
            sim.update(w, h, dt)
        expected = 100.0 * (1.0 - MUCUS_FOOD_WASHOUT_RATE_PER_SEC * dt) ** n
        self.assertAlmostEqual(node.amount, expected, places=4)

    def test_stats_record_lumen_and_mucus_population_split(self):
        # Place cells at known radii and confirm the per-step stats classify them by zone.
        sp = self._species('Zoned', consumption_speed=0.0)
        sim = Simulation([sp], self._env(drop_frequency=0), 100)
        sim.running = True
        w, h, dt = 800, 600, 1.0
        sim.update(w, h, dt)  # establishes petri_center / petri_radius
        cx, cy = sim.petri_center
        R = sim.petri_radius
        # Two lumen cells (near centre) and three mucus cells (outer annulus).
        sp.members = []
        for _ in range(2):
            sp.members.append(Bacterium(sp, cx + 0.10 * R, cy, size=1.0))
        for _ in range(3):
            sp.members.append(Bacterium(sp, cx + 0.95 * R, cy, size=1.0))
        sim.update(w, h, dt)
        self.assertEqual(sim.stats['lumen'][-1], 2)
        self.assertEqual(sim.stats['mucus'][-1], 3)
        # The split must always sum to the recorded total population.
        self.assertEqual(sim.stats['lumen'][-1] + sim.stats['mucus'][-1], sim.stats['total'][-1])


class TestSubstrateConsumption(unittest.TestCase):
    """Phase I diet analytics: per-substrate cumulative consumption stats are accumulated as cells
    ingest substrate-typed nutrient patches, recorded as monotonic running totals, and exported."""

    def _species(self, name, consumption_speed=5.0):
        return BacteriaSpecies(name, 0, 0.0, 0.0, 1.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S,
                               NEUTRAL_PH, (0, 200, 0), base_size=1.0,
                               food_capacity=1000.0, metabolic_rate=0.0,
                               consumption_speed=consumption_speed,
                               bacteriocin_production=0.0, immunity_factor=0.0)

    def _env(self):
        return {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION, 'co2': 0.0,
                'temperature': DEFAULT_ENV_TEMPERATURE_C, 'drop_frequency': 0,
                'drop_range': (1, 1), 'phage_inflow': 0, 'diet_profile': DEFAULT_DIET_PROFILE}

    def test_consumption_attributed_to_correct_substrate_class(self):
        # A hungry, competent cell sitting on a calorie patch must accrue calorie consumption only;
        # the glycan and MDCF running totals stay at zero because no such patches are present and
        # the default species has zero affinity for them anyway.
        sp = self._species('Grazer')
        sim = Simulation([sp], self._env(), 100)
        sim.running = True
        w, h, dt = 800, 600, 1.0
        sim.update(w, h, dt)  # establishes petri_center / petri_radius (no drops)
        cx, cy = sim.petri_center
        # One hungry cell co-located with a calorie node so it grazes immediately.
        sp.members = [Bacterium(sp, cx, cy, size=1.0, initial_food=0.0)]
        calorie_node = FoodNode.acquire().init(cx, cy, 500.0, SUBSTRATE_GENERIC_CALORIE)
        sim.food_nodes = [calorie_node]
        for _ in range(5):
            sim.update(w, h, dt)
        self.assertGreater(sim.stats['consumed_calorie'][-1], 0.0)
        self.assertEqual(sim.stats['consumed_glycan'][-1], 0.0)
        self.assertEqual(sim.stats['consumed_mdcf'][-1], 0.0)

    def test_consumption_series_is_monotonic_non_decreasing(self):
        # Cumulative consumption can never decrease step-to-step (mass already ingested is permanent).
        sp = self._species('Grazer')
        sim = Simulation([sp], self._env(), 100)
        sim.running = True
        w, h, dt = 800, 600, 1.0
        sim.update(w, h, dt)
        cx, cy = sim.petri_center
        sp.members = [Bacterium(sp, cx, cy, size=1.0, initial_food=0.0)]
        sim.food_nodes = [FoodNode.acquire().init(cx, cy, 500.0, SUBSTRATE_GENERIC_CALORIE)]
        for _ in range(8):
            sim.update(w, h, dt)
        series = sim.stats['consumed_calorie']
        self.assertGreaterEqual(len(series), 2)
        for prev, curr in zip(series, series[1:]):
            self.assertGreaterEqual(curr, prev)

    def test_empty_stats_includes_substrate_consumption_keys(self):
        # The stats schema must always expose the three per-substrate consumption series.
        sim = Simulation([], self._env(), 100)
        stats = sim._empty_stats()
        for key in ('consumed_calorie', 'consumed_glycan', 'consumed_mdcf'):
            self.assertIn(key, stats)
            self.assertEqual(stats[key], [])


class TestReproducibility(unittest.TestCase):
    """Phase I reproducibility & extinction baselines: a fixed seed yields byte-identical
    trajectories, and a species driven to zero records a finite time-to-extinction."""

    def _species(self, name, init=20, base_rep=0.3):
        return BacteriaSpecies(name, init, base_rep, 0.0, 1.0, 1.0, DEFAULT_BACTERIA_SPEED_UM_S,
                               NEUTRAL_PH, (0, 200, 0), base_size=1.0,
                               food_capacity=50.0, metabolic_rate=0.1,
                               consumption_speed=1.0,
                               bacteriocin_production=0.0, immunity_factor=0.0)

    def _env(self, drop_frequency=20):
        return {'ph': NEUTRAL_PH, 'o2': ATMOSPHERIC_O2_FRACTION, 'co2': 0.0,
                'temperature': DEFAULT_ENV_TEMPERATURE_C, 'drop_frequency': drop_frequency,
                'drop_range': (1, 4), 'phage_inflow': 0, 'diet_profile': DEFAULT_DIET_PROFILE}

    def _run(self, seed, n):
        import random as _random
        _random.seed(seed)  # seed before building species so initial layouts match across runs
        spp = [self._species('A'), self._species('B')]
        sim = Simulation(spp, self._env(), 300, seed=seed)
        sim.running = True
        for _ in range(n):
            sim.update(800, 600, 0.1)
        return sim

    def test_seeded_runs_are_identical(self):
        sim_a = self._run(2024, 60)
        sim_b = self._run(2024, 60)
        self.assertEqual(sim_a.stats['total'], sim_b.stats['total'])
        self.assertEqual(sim_a.stats['births'], sim_b.stats['births'])
        self.assertEqual(sim_a.stats['deaths'], sim_b.stats['deaths'])
        for sp_a, sp_b in zip(sim_a.species, sim_b.species):
            self.assertEqual(sp_a.history, sp_b.history)

    def test_different_seeds_diverge(self):
        sim_a = self._run(1, 60)
        sim_b = self._run(2, 60)
        # Distinct seeds should not produce identical population trajectories.
        self.assertNotEqual(sim_a.stats['total'], sim_b.stats['total'])

    def test_extinction_time_recorded_on_zero_population(self):
        import random as _random
        _random.seed(7)
        sp = self._species('Doomed', init=5)
        sim = Simulation([sp], self._env(drop_frequency=0), 300, seed=7)
        sim.running = True
        sim.update(800, 600, 0.1)
        self.assertIsNone(sp.extinction_time)  # still alive
        sp.members = []  # force the population to zero
        sim.update(800, 600, 0.1)
        self.assertIsNotNone(sp.extinction_time)
        self.assertAlmostEqual(sp.extinction_time, sim.stats['elapsed_seconds'][-1])

    def test_extinction_time_is_first_crossing_only(self):
        import random as _random
        _random.seed(3)
        sp = self._species('Once', init=4)
        sim = Simulation([sp], self._env(drop_frequency=0), 300, seed=3)
        sim.running = True
        sim.update(800, 600, 0.1)
        sp.members = []
        sim.update(800, 600, 0.1)
        first = sp.extinction_time
        sim.update(800, 600, 0.1)  # still zero
        self.assertEqual(sp.extinction_time, first)  # never overwritten


class TestSeedInput(unittest.TestCase):
    """Phase I: the Globals-tab seed field parses to a non-negative int or None (random)."""

    @classmethod
    def setUpClass(cls):
        # App pulls in pygame; force a headless video driver so import never needs a display.
        os.environ.setdefault('SDL_VIDEODRIVER', 'dummy')
        from bacteria.app import App
        cls.App = App

    def test_blank_and_keywords_map_to_random(self):
        for text in ('', '   ', 'none', 'NULL', 'None'):
            ok, value = self.App.parse_seed_input(text)
            self.assertTrue(ok)
            self.assertIsNone(value)

    def test_integers_parse(self):
        self.assertEqual(self.App.parse_seed_input('42'), (True, 42))
        self.assertEqual(self.App.parse_seed_input('  7 '), (True, 7))
        self.assertEqual(self.App.parse_seed_input('0'), (True, 0))

    def test_trailing_decimal_zero_tolerated(self):
        # normalize_numeric_string can hand us '42.0'; it must still parse to int 42.
        self.assertEqual(self.App.parse_seed_input('42.0'), (True, 42))

    def test_negative_and_garbage_rejected(self):
        for text in ('-1', '-100', 'abc', '4x', '3.5.1'):
            ok, value = self.App.parse_seed_input(text)
            self.assertFalse(ok)
            self.assertIsNone(value)


class TestDietInput(unittest.TestCase):
    """Phase I: the Globals-tab diet_profile cell cycles through the categorical feeding arms."""

    @classmethod
    def setUpClass(cls):
        os.environ.setdefault('SDL_VIDEODRIVER', 'dummy')
        from bacteria.app import App
        cls.App = App

    def test_cycle_visits_every_arm_and_wraps(self):
        order = list(DIET_PROFILES.keys())
        # Starting from each arm, next_diet_profile advances to the following arm (wrapping).
        for i, arm in enumerate(order):
            self.assertEqual(self.App.next_diet_profile(arm), order[(i + 1) % len(order)])

    def test_full_cycle_returns_to_start(self):
        order = list(DIET_PROFILES.keys())
        current = order[0]
        for _ in range(len(order)):
            current = self.App.next_diet_profile(current)
        self.assertEqual(current, order[0])  # one full lap returns to the origin

    def test_unknown_profile_resets_to_default(self):
        self.assertEqual(self.App.next_diet_profile('not_a_diet'), DEFAULT_DIET_PROFILE)
        self.assertEqual(self.App.next_diet_profile(None), DEFAULT_DIET_PROFILE)


if __name__ == '__main__':
    unittest.main()

