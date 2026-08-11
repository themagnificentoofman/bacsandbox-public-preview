import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from bacteria.constants import SUBSTRATE_NAMES
from bacteria.constants import ATMOSPHERIC_O2_FRACTION, NEUTRAL_PH
from bacteria.core_sim import Bacterium, BacteriaSpecies, Simulation
from bacteria.diet_translator import DietTranslationEngine
from bacteria.ecm import ECM, SubstrateLayerSpec
from bacteria.metabolic_integration import (
    DietBolus,
    FbaAgentState,
    FbaConservationError,
    FbaIntegrationAdapter,
    FbaPhysiology,
    FbaSpatialConfig,
    FbaTickDisposition,
    GrowthOnlyFbaLifecyclePolicy,
    build_coordinate_environment_batch,
    deposit_diet_bolus,
)
from bacteria.metabolism import (
    MetabolicOrchestrator,
    build_species_phenotype_config,
    load_metabolic_model,
    load_metabolic_priors,
)


def synthetic_substrate_specs(diffusion_um2_s=0.0):
    return tuple(
        SubstrateLayerSpec(
            name=name,
            diffusion_um2_s=diffusion_um2_s,
            decay_rate_s_inv=0.0,
            min_concentration_m=1.0e-15,
            provenance='Own (BacSandbox synthetic test)',
        )
        for name in SUBSTRATE_NAMES
    )


class SyntheticAgent:
    def __init__(self, species, x, y):
        self.species = species
        self.x = x
        self.y = y


class SyntheticSpecies:
    def __init__(self, aff_calorie=1.0, aff_glycan=1.0, aff_mdcf=1.0):
        self.aff_calorie = aff_calorie
        self.aff_glycan = aff_glycan
        self.aff_mdcf = aff_mdcf
        self.members = []


class FixedBiomassLifecycle:
    def __init__(self, biomass_gdw=1.0e-9):
        self.biomass_gdw = biomass_gdw

    def initialize_agent(self, bacterium):
        return FbaAgentState(biomass_gdw=self.biomass_gdw)

    def evaluate(
        self,
        bacterium,
        state,
        *,
        dt_seconds,
        environmental_fitness,
        carrying_capacity_factor,
    ):
        return FbaPhysiology(
            render_scale=1.0,
            motility_scale=0.0,
            chemotactic_drive=0.0,
            division_probability=0.0,
            starved=False,
        )

    def partition_for_fission(self, parent, child, state):
        partitioned = FbaAgentState(
            biomass_gdw=state.biomass_gdw * 0.5,
            growth_rate_s=state.growth_rate_s,
            biomass_delta_gdw=state.biomass_delta_gdw * 0.5,
            substrate_consumed_mol=tuple(
                value * 0.5 for value in state.substrate_consumed_mol
            ),
            acetate_secreted_mol=state.acetate_secreted_mol * 0.5,
        )
        return partitioned, partitioned


class TestGrowthOnlyFbaLifecyclePolicy(unittest.TestCase):
    def test_divides_on_biomass_doubling_without_inventing_starvation(self):
        policy = GrowthOnlyFbaLifecyclePolicy(initial_biomass_gdw=1.0e-9)
        agent = object()
        initial = policy.initialize_agent(agent)

        quiescent = policy.evaluate(
            agent,
            initial,
            dt_seconds=1.0,
            environmental_fitness=0.0,
            carrying_capacity_factor=1.0,
        )
        self.assertEqual(quiescent.division_probability, 0.0)
        self.assertFalse(quiescent.starved)

        doubled = FbaAgentState(biomass_gdw=2.0e-9)
        dividing = policy.evaluate(
            agent,
            doubled,
            dt_seconds=1.0,
            environmental_fitness=0.0,
            carrying_capacity_factor=1.0,
        )
        self.assertEqual(dividing.division_probability, 1.0)
        self.assertFalse(dividing.starved)

        parent_state, child_state = policy.partition_for_fission(
            agent,
            object(),
            doubled,
        )
        self.assertEqual(parent_state.biomass_gdw, 1.0e-9)
        self.assertEqual(child_state.biomass_gdw, 1.0e-9)
        self.assertEqual(
            parent_state.biomass_gdw + child_state.biomass_gdw,
            doubled.biomass_gdw,
        )


def make_python_adapter(species, *, event_sink=None, conservation_atol_m=1.0e-15):
    model = load_metabolic_model()
    priors = load_metabolic_priors()
    phenotype_config = build_species_phenotype_config(species, model=model, priors=priors)
    orchestrator = MetabolicOrchestrator(
        mode='fba',
        backend='python',
        model=model,
        priors=priors,
        phenotype_config=phenotype_config,
    )
    config = FbaSpatialConfig(
        environment_volume_l=1.0e-3,
        substrate_layers=synthetic_substrate_specs(),
        conservation_atol_m=conservation_atol_m,
    )
    return FbaIntegrationAdapter(
        orchestrator,
        DietTranslationEngine(),
        config,
        {profile: index for index, profile in enumerate(species)},
        FixedBiomassLifecycle(),
        event_sink=event_sink,
    )


class TestContinuousSubstrateLayers(unittest.TestCase):
    def test_registration_deposit_and_batch_sampling_are_strict(self):
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(synthetic_substrate_specs())
        self.assertEqual(ecm.substrate_layer_names, SUBSTRATE_NAMES)
        self.assertEqual(tuple(ecm._layers), ('qs', 'bac', 'acetate', 'phage'))

        coordinate = ecm.world_to_cell(-0.1, 10.0)
        self.assertEqual(coordinate, (-1, 1))
        increment = np.asarray((0.2, 0.3, 0.4), dtype=np.float64)
        ecm.deposit_substrates_m(coordinate, increment)
        sampled = ecm.sample_substrate_batch_m((coordinate, (4, 5)))

        np.testing.assert_array_equal(sampled[0], increment)
        np.testing.assert_array_equal(sampled[1], np.zeros(3, dtype=np.float64))
        self.assertEqual(sampled.dtype, np.float64)
        self.assertTrue(sampled.flags.c_contiguous)
        self.assertTrue(sampled.flags.owndata)
        self.assertTrue(sampled.flags.writeable)

        with self.assertRaisesRegex(ValueError, 'unique'):
            ecm.sample_substrate_batch_m((coordinate, coordinate))
        with self.assertRaisesRegex(TypeError, 'float64'):
            ecm.deposit_substrates_m(coordinate, increment.astype(np.float32))

    def test_positive_subfloor_deposit_is_rejected_without_mass_loss(self):
        ecm = ECM(cell_size=10.0, backend='python')
        specs = tuple(
            SubstrateLayerSpec(
                name=name,
                diffusion_um2_s=0.0,
                decay_rate_s_inv=0.0,
                min_concentration_m=1.0e-3,
                provenance='Own (BacSandbox synthetic test)',
            )
            for name in SUBSTRATE_NAMES
        )
        ecm.register_substrate_layers(specs)
        before = ecm.snapshot_substrate_layers()

        with self.assertRaisesRegex(ValueError, 'retention floor'):
            ecm.deposit_substrates_m(
                (0, 0),
                np.asarray((1.0e-4, 0.0, 0.0), dtype=np.float64),
            )

        self.assertEqual(ecm.snapshot_substrate_layers(), before)

    def test_replacement_is_atomic_and_snapshot_restores_complete_maps(self):
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(synthetic_substrate_specs())
        coordinates = ((0, 0), (1, 0))
        initial = np.asarray(((0.2, 0.3, 0.4), (0.5, 0.6, 0.7)), dtype=np.float64)
        for coordinate, increment in zip(coordinates, initial):
            ecm.deposit_substrates_m(coordinate, increment)
        snapshot = ecm.snapshot_substrate_layers()

        invalid = initial.copy()
        invalid[1, 2] = np.nan
        with self.assertRaisesRegex(ValueError, 'finite nonnegative'):
            ecm.replace_substrate_batch_m(coordinates, invalid)
        np.testing.assert_array_equal(ecm.sample_substrate_batch_m(coordinates), initial)

        replacement = np.asarray(((0.1, 0.0, 0.2), (0.0, 0.4, 0.0)), dtype=np.float64)
        ecm.replace_substrate_batch_m(coordinates, replacement)
        np.testing.assert_array_equal(ecm.sample_substrate_batch_m(coordinates), replacement)
        ecm.restore_substrate_layers(snapshot)
        np.testing.assert_array_equal(ecm.sample_substrate_batch_m(coordinates), initial)

    def test_python_diffusion_conserves_concentration_sum_without_decay(self):
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(synthetic_substrate_specs(diffusion_um2_s=1.0))
        origin = (0, 0)
        ecm.deposit_substrates_m(
            origin,
            np.asarray((1.0, 0.0, 0.0), dtype=np.float64),
        )

        ecm.tick_substrates(1.0)
        generic_layer = ecm._substrate_layers[SUBSTRATE_NAMES[0]]
        self.assertEqual(set(generic_layer), {(0, 0), (-1, 0), (1, 0), (0, -1), (0, 1)})
        self.assertAlmostEqual(sum(generic_layer.values()), 1.0, delta=1.0e-15)
        self.assertEqual(ecm._substrate_layers[SUBSTRATE_NAMES[1]], {})
        self.assertEqual(ecm._substrate_layers[SUBSTRATE_NAMES[2]], {})


class TestDietBolusDeposition(unittest.TestCase):
    def test_weighted_footprint_conserves_translated_amount(self):
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(synthetic_substrate_specs())
        engine = DietTranslationEngine()
        bolus = DietBolus(
            event_id='mdcf-test',
            diet_name='mdcf_2',
            feed_mass_g=50.0,
            footprint=(((0, 0), 0.25), ((1, 0), 0.75)),
        )
        deposited = deposit_diet_bolus(
            engine,
            ecm,
            bolus,
            environment_volume_l=1.0,
            dt_seconds=0.1,
        )
        concentrations = ecm.sample_substrate_batch_m(((0, 0), (1, 0)))
        represented_mol = np.sum(concentrations, axis=0)

        np.testing.assert_allclose(
            represented_mol,
            np.asarray(deposited.total_substrate_mol),
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            concentrations[1],
            concentrations[0] * 3.0,
            rtol=0.0,
            atol=1.0e-15,
        )

    def test_tiny_bolus_below_retention_floor_is_rejected_atomically(self):
        ecm = ECM(cell_size=10.0, backend='python')
        specs = tuple(
            SubstrateLayerSpec(
                name=name,
                diffusion_um2_s=0.0,
                decay_rate_s_inv=0.0,
                min_concentration_m=1.0,
                provenance='Own (BacSandbox synthetic test)',
            )
            for name in SUBSTRATE_NAMES
        )
        ecm.register_substrate_layers(specs)
        before = ecm.snapshot_substrate_layers()
        with self.assertRaisesRegex(ValueError, 'retention floor'):
            deposit_diet_bolus(
                DietTranslationEngine(),
                ecm,
                DietBolus(
                    event_id='subfloor',
                    diet_name='mdcf_2',
                    feed_mass_g=1.0,
                    footprint=(((0, 0), 0.999999), ((1, 0), 0.000001)),
                ),
                environment_volume_l=0.1,
                dt_seconds=0.1,
            )
        self.assertEqual(ecm.snapshot_substrate_layers(), before)


class TestFbaIntegrationAdapter(unittest.TestCase):
    def test_batch_order_is_deterministic_and_shared_cells_share_environment(self):
        first = SyntheticSpecies()
        second = SyntheticSpecies()
        first.members = (SyntheticAgent(first, 15.0, 0.0), SyntheticAgent(first, 1.0, 1.0))
        second.members = (SyntheticAgent(second, 2.0, 2.0),)
        species = (first, second)
        adapter = make_python_adapter(species)
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(synthetic_substrate_specs())

        batch = build_coordinate_environment_batch(
            species,
            ecm,
            adapter.state_store,
            adapter.phenotype_index,
            adapter.config.environment_volume_l,
        )

        self.assertEqual(batch.coordinates, ((0, 0), (1, 0)))
        np.testing.assert_array_equal(batch.environment_id, np.asarray((1, 0, 0), dtype=np.int32))
        np.testing.assert_array_equal(batch.phenotype_id, np.asarray((0, 0, 1), dtype=np.int32))

    def test_successful_tick_commits_grid_biomass_and_event(self):
        profile = SyntheticSpecies()
        agent = SyntheticAgent(profile, 1.0, 1.0)
        profile.members = (agent,)
        events = []
        adapter = make_python_adapter((profile,), event_sink=events.append)
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(synthetic_substrate_specs())
        coordinate = (0, 0)
        initial = np.asarray((0.1, 0.1, 0.1), dtype=np.float64)
        ecm.deposit_substrates_m(coordinate, initial)

        result = adapter.advance_tick(
            species=(profile,),
            ecm=ecm,
            tick_index=0,
            dt_seconds=0.01,
        )

        self.assertEqual(result.disposition, FbaTickDisposition.COMMITTED)
        self.assertEqual(result.backend_used, 'python')
        self.assertLess(result.max_conservation_error_mol, 1.0e-18)
        self.assertGreater(result.biomass_delta_gdw, 0.0)
        self.assertGreater(result.acetate_secreted_mol, 0.0)
        self.assertTrue(np.all(ecm.sample_substrate_batch_m((coordinate,))[0] < initial))
        self.assertGreater(adapter.state_store.state_for(agent).biomass_gdw, 1.0e-9)
        self.assertEqual(events[-1].disposition, FbaTickDisposition.COMMITTED)

    def test_solver_failure_restores_grid_and_sidecar_state(self):
        profile = SyntheticSpecies()
        agent = SyntheticAgent(profile, 1.0, 1.0)
        profile.members = (agent,)
        events = []
        adapter = make_python_adapter((profile,), event_sink=events.append)
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(synthetic_substrate_specs())
        initial = np.asarray((0.1, 0.2, 0.3), dtype=np.float64)
        ecm.deposit_substrates_m((0, 0), initial)
        adapter.state_store.reconcile((agent,))
        original_state = adapter.state_store.state_for(agent)

        def fail_after_write(**batch):
            batch['concentrations_m'].fill(0.0)
            batch['biomass_gdw'].fill(99.0)
            raise RuntimeError('synthetic solver failure')

        adapter.orchestrator.advance_batch = fail_after_write
        with self.assertRaisesRegex(RuntimeError, 'synthetic solver failure'):
            adapter.advance_tick(
                species=(profile,),
                ecm=ecm,
                tick_index=1,
                dt_seconds=0.01,
                boluses=(DietBolus(
                    event_id='rollback-bolus',
                    diet_name='mdcf_2',
                    feed_mass_g=1.0,
                    footprint=(((0, 0), 1.0),),
                ),),
            )

        np.testing.assert_array_equal(ecm.sample_substrate_batch_m(((0, 0),))[0], initial)
        self.assertEqual(adapter.state_store.state_for(agent), original_state)
        self.assertEqual(events[-1].disposition, FbaTickDisposition.ROLLED_BACK)
        self.assertEqual(events[-1].phase, 'solve')

    def test_first_use_failure_restores_unregistered_ecm_state(self):
        profile = SyntheticSpecies()
        profile.members = (SyntheticAgent(profile, 1.0, 1.0),)
        adapter = make_python_adapter((profile,))
        ecm = ECM(cell_size=10.0, backend='python')
        before = ecm.snapshot_substrate_layers()

        def fail(**batch):
            raise RuntimeError('first-use failure')

        adapter.orchestrator.advance_batch = fail
        with self.assertRaisesRegex(RuntimeError, 'first-use failure'):
            adapter.advance_tick(
                species=(profile,),
                ecm=ecm,
                tick_index=0,
                dt_seconds=0.01,
            )

        self.assertEqual(ecm.snapshot_substrate_layers(), before)
        self.assertEqual(ecm.substrate_layer_names, ())


class TestLiveFbaSimulation(unittest.TestCase):
    @staticmethod
    def make_species():
        profile = BacteriaSpecies(
            'Synthetic live FBA species',
            0,
            0.1,
            0.01,
            0.0,
            1.0,
            20.0,
            NEUTRAL_PH,
            (100, 150, 200),
            aff_calorie=1.0,
            aff_glycan=1.0,
            aff_mdcf=1.0,
        )
        profile.members = [Bacterium(profile, 0.0, 0.0)]
        return profile

    @staticmethod
    def environment():
        return {
            'ph': NEUTRAL_PH,
            'o2': ATMOSPHERIC_O2_FRACTION,
            'co2': 0.0,
            'temperature': 37.0,
            'drop_frequency': 100.0,
            'drop_range': (1.0, 1.0),
        }

    def test_live_fba_tick_bypasses_food_nodes_and_pg_reserve(self):
        profile = self.make_species()
        agent = profile.members[0]
        adapter = make_python_adapter((profile,))
        simulation = Simulation(
            [profile],
            self.environment(),
            cap=100,
            seed=11,
            ecm_backend='python',
            metabolic_adapter=adapter,
        )
        simulation._populated = True
        simulation.petri_center = (0.0, 0.0)
        simulation.petri_radius = 100.0
        simulation.running = True
        simulation.queue_diet_bolus(DietBolus(
            event_id='live-bolus',
            diet_name='mdcf_2',
            feed_mass_g=1.0,
            footprint=(((0, 0), 1.0),),
        ))
        initial_food_pg = agent.food_level

        result = simulation.update(100, 100, 0.01)

        self.assertEqual(result.disposition, FbaTickDisposition.COMMITTED)
        self.assertEqual(simulation.step_idx, 1)
        self.assertEqual(simulation.food_nodes, [])
        self.assertEqual(agent.food_level, initial_food_pg)
        self.assertGreater(adapter.state_store.state_for(agent).biomass_gdw, 1.0e-9)
        self.assertGreater(simulation.stats['fba_biomass_gdw'][-1], 1.0e-9)
        self.assertGreater(simulation.stats['fba_acetate_mol'][-1], 0.0)
        self.assertEqual(simulation._pending_diet_boluses, [])

    def test_failed_live_fba_tick_retains_pending_bolus_and_tick_index(self):
        profile = self.make_species()
        adapter = make_python_adapter((profile,))
        simulation = Simulation(
            [profile],
            self.environment(),
            cap=100,
            seed=12,
            ecm_backend='python',
            metabolic_adapter=adapter,
        )
        simulation._populated = True
        simulation.petri_center = (0.0, 0.0)
        simulation.petri_radius = 100.0
        simulation.running = True
        bolus = DietBolus(
            event_id='retry-bolus',
            diet_name='mdcf_2',
            feed_mass_g=1.0,
            footprint=(((0, 0), 1.0),),
        )
        simulation.queue_diet_bolus(bolus)

        def fail(**batch):
            raise RuntimeError('live synthetic failure')

        adapter.orchestrator.advance_batch = fail
        with self.assertRaisesRegex(RuntimeError, 'live synthetic failure'):
            simulation.update(100, 100, 0.01)

        self.assertEqual(simulation.step_idx, 0)
        self.assertEqual(simulation._pending_diet_boluses, [bolus])
        self.assertEqual(simulation.stats['time'], [])

    def test_reset_clears_live_fba_grid_state_queue_result_and_cache(self):
        profile = self.make_species()
        adapter = make_python_adapter((profile,))
        simulation = Simulation(
            [profile],
            self.environment(),
            cap=100,
            seed=13,
            ecm_backend='python',
            metabolic_adapter=adapter,
        )
        simulation._populated = True
        simulation.petri_center = (0.0, 0.0)
        simulation.petri_radius = 100.0
        simulation.running = True
        simulation.queue_diet_bolus(DietBolus(
            event_id='first-run',
            diet_name='mdcf_2',
            feed_mass_g=1.0,
            footprint=(((0, 0), 1.0),),
        ))
        simulation.update(100, 100, 0.01)
        simulation.queue_diet_bolus(DietBolus(
            event_id='must-clear',
            diet_name='mdcf_2',
            feed_mass_g=1.0,
            footprint=(((0, 0), 1.0),),
        ))
        self.assertTrue(any(simulation.ecm._substrate_layers.values()))
        self.assertTrue(adapter.state_store.snapshot())
        self.assertIsNotNone(simulation.last_fba_result)

        simulation.reset()

        self.assertEqual(simulation._pending_diet_boluses, [])
        self.assertIsNone(simulation.last_fba_result)
        self.assertFalse(any(simulation.ecm._substrate_layers.values()))
        self.assertEqual(adapter.state_store.snapshot(), {})
        self.assertEqual(adapter.orchestrator._kernel.cache_stats(), (0, 0, 0, 0, 0))
        self.assertEqual(simulation.stats['fba_biomass_gdw'], [])

    def test_export_includes_live_fba_metrics(self):
        profile = self.make_species()
        adapter = make_python_adapter((profile,))
        simulation = Simulation(
            [profile],
            self.environment(),
            cap=100,
            seed=14,
            ecm_backend='python',
            metabolic_adapter=adapter,
        )
        simulation._populated = True
        simulation.petri_center = (0.0, 0.0)
        simulation.petri_radius = 100.0
        simulation.running = True
        simulation.queue_diet_bolus(DietBolus(
            event_id='export',
            diet_name='mdcf_2',
            feed_mass_g=1.0,
            footprint=(((0, 0), 1.0),),
        ))
        simulation.update(100, 100, 0.01)

        with tempfile.TemporaryDirectory() as directory:
            export_path = Path(directory) / 'stats.csv'
            with patch('bacteria.core_sim.EXPORT_FILE', str(export_path)):
                simulation.export()
            with export_path.open(newline='', encoding='utf-8') as handle:
                rows = list(csv.reader(handle))

        header = next(row for row in rows if row and row[0] == 'time')
        data = rows[rows.index(header) + 1]
        for field_name in (
            'fba_biomass_gdw',
            'fba_acetate_mol',
            'fba_consumed_calorie_mol',
            'fba_consumed_glycan_mol',
            'fba_consumed_mdcf_mol',
        ):
            index = header.index(field_name)
            self.assertNotEqual(data[index], '')

    def test_conservation_failure_rolls_back_before_grid_commit(self):
        profile = SyntheticSpecies()
        agent = SyntheticAgent(profile, 1.0, 1.0)
        profile.members = (agent,)
        adapter = make_python_adapter((profile,), conservation_atol_m=0.0)
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(synthetic_substrate_specs())
        initial = np.asarray((0.1, 0.1, 0.1), dtype=np.float64)
        ecm.deposit_substrates_m((0, 0), initial)

        def violate_conservation(**batch):
            batch['concentrations_m'] -= 0.01
            return (1, 0, 0, 0, 1, 0.01)

        adapter.orchestrator.advance_batch = violate_conservation
        with self.assertRaises(FbaConservationError):
            adapter.advance_tick(
                species=(profile,),
                ecm=ecm,
                tick_index=2,
                dt_seconds=0.01,
            )
        np.testing.assert_array_equal(ecm.sample_substrate_batch_m(((0, 0),))[0], initial)


if __name__ == '__main__':
    unittest.main()