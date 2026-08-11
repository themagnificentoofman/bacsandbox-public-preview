import copy
import gc
import importlib
import logging
import math
import random
import tracemalloc
import unittest
import weakref
from unittest import mock

import numpy as np

from bacteria.constants import ATMOSPHERIC_O2_FRACTION, NEUTRAL_PH, SUBSTRATE_NAMES
from bacteria.core_sim import Bacterium, BacteriaSpecies, Simulation
from bacteria.diet_translator import DietTranslationEngine
from bacteria.ecm import ECM, SubstrateLayerSpec
from bacteria.metabolic_integration import (
    DietBolus,
    FbaConservationError,
    FbaIntegrationAdapter,
    FbaSpatialConfig,
    FbaTickDisposition,
    GrowthOnlyFbaLifecyclePolicy,
    build_coordinate_environment_batch,
)
from bacteria.metabolism import (
    MetabolicOrchestrator,
    build_species_phenotype_config,
    load_metabolic_model,
    load_metabolic_priors,
)


try:
    bacteria_c_kernel = importlib.import_module('bacteria_c_kernel')
except ModuleNotFoundError as error:
    if error.name != 'bacteria_c_kernel':
        raise
    bacteria_c_kernel = None


NATIVE_METABOLIC_AVAILABLE = (
    bacteria_c_kernel is not None and hasattr(bacteria_c_kernel, 'MetabolicKernel')
)
BIOMASS_ASSERT_ATOL_GDW = 1.0e-18
RETENTION_LIMIT_BYTES = 1_048_576
QUIET_TEST_LOGGER = logging.getLogger('bacsandbox.integration-stress')
QUIET_TEST_LOGGER.disabled = True


def synthetic_substrate_specs(diffusion_um2_s=0.0):
    return tuple(
        SubstrateLayerSpec(
            name=name,
            diffusion_um2_s=diffusion_um2_s,
            decay_rate_s_inv=0.0,
            min_concentration_m=1.0e-15,
            provenance='Own (BacSandbox synthetic stress fixture)',
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


def make_adapter(
    species,
    *,
    backend='python',
    environment_volume_l=1.0e-3,
    diffusion_um2_s=0.0,
    initial_biomass_gdw=1.0e-9,
    event_sink=None,
    logger=None,
):
    species = tuple(species)
    if logger is None:
        logger = QUIET_TEST_LOGGER
    model = load_metabolic_model()
    priors = load_metabolic_priors()
    orchestrator = MetabolicOrchestrator(
        mode='fba',
        backend=backend,
        model=model,
        priors=priors,
        phenotype_config=build_species_phenotype_config(
            species,
            model=model,
            priors=priors,
        ),
    )
    return FbaIntegrationAdapter(
        orchestrator,
        DietTranslationEngine(),
        FbaSpatialConfig(
            environment_volume_l=environment_volume_l,
            substrate_layers=synthetic_substrate_specs(diffusion_um2_s),
        ),
        {profile: index for index, profile in enumerate(species)},
        GrowthOnlyFbaLifecyclePolicy(initial_biomass_gdw=initial_biomass_gdw),
        event_sink=event_sink,
        logger=logger,
    )


def make_live_fixture():
    species = BacteriaSpecies(
        'Synthetic integration stress species',
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
    agent = Bacterium(species, 0.0, 0.0)
    species.members = [agent]
    adapter = make_adapter((species,))
    simulation = Simulation(
        [species],
        {
            'ph': NEUTRAL_PH,
            'o2': ATMOSPHERIC_O2_FRACTION,
            'co2': 0.0,
            'temperature': 37.0,
            'drop_frequency': 0.0,
            'drop_range': (0.0, 0.0),
        },
        cap=100,
        seed=4101,
        ecm_backend='python',
        metabolic_adapter=adapter,
    )
    simulation._populated = True
    simulation.petri_center = (0.0, 0.0)
    simulation.petri_radius = 100.0
    simulation.running = True
    simulation.ecm.register_substrate_layers(adapter.config.substrate_layers)
    simulation.ecm.deposit_substrates_m(
        (-3, 4),
        np.asarray((0.125, 0.25, 0.375), dtype=np.float64),
    )
    simulation.ecm.deposit_substrates_m(
        (7, -2),
        np.asarray((0.5, 0.625, 0.75), dtype=np.float64),
    )
    adapter.state_store.reconcile((agent,))
    bolus = DietBolus(
        event_id='transaction-retry',
        diet_name='mdcf_2',
        feed_mass_g=1.0,
        footprint=(((0, 0), 0.4), ((1, 0), 0.6)),
    )
    simulation.queue_diet_bolus(bolus)
    return simulation, species, agent, bolus


def make_synthetic_community(agent_count, *, species_count=1, shared_cell=True):
    profiles = tuple(SyntheticSpecies() for _ in range(species_count))
    agents = []
    remaining = agent_count
    for profile_index, profile in enumerate(profiles):
        profile_count = remaining // (species_count - profile_index)
        remaining -= profile_count
        profile.members = [
            SyntheticAgent(
                profile,
                1.0 if shared_cell else 10.0 * (len(agents) + index),
                1.0,
            )
            for index in range(profile_count)
        ]
        agents.extend(profile.members)
    return profiles, tuple(agents)


def snapshot_live_simulation(simulation, species, agent):
    memo = {
        id(simulation): simulation,
        id(species): species,
        id(agent): agent,
    }
    return {
        'substrates': simulation.ecm.snapshot_substrate_layers(),
        'signals': tuple(
            (name, tuple(sorted(layer.items())))
            for name, layer in simulation.ecm._layers.items()
        ),
        'states': simulation.metabolic_adapter.state_store.snapshot(),
        'species': copy.deepcopy(species.__dict__, memo),
        'agent': copy.deepcopy(agent.__dict__, memo),
        'stats': copy.deepcopy(simulation.stats),
        'species_step_stats': copy.deepcopy(simulation.species_step_stats),
        'environment': copy.deepcopy(simulation.env),
        'pending_boluses': tuple(simulation._pending_diet_boluses),
        'last_fba_result': simulation.last_fba_result,
        'step_idx': simulation.step_idx,
        'random_state': random.getstate(),
    }


class InjectedPhaseError(RuntimeError):
    pass


class AppendThenRaiseList(list):
    def append(self, value):
        super().append(value)
        raise InjectedPhaseError('injected statistics commit failure')


class TestIntegrationOwnerConservation(unittest.TestCase):
    def test_each_environment_conserves_substrate_and_separates_output_units(self):
        first = SyntheticSpecies()
        second = SyntheticSpecies(aff_calorie=0.8, aff_glycan=1.2, aff_mdcf=1.5)
        first.members = [
            SyntheticAgent(first, 1.0, 1.0),
            SyntheticAgent(first, 11.0, 1.0),
        ]
        second.members = [
            SyntheticAgent(second, 2.0, 2.0),
            SyntheticAgent(second, 21.0, 1.0),
        ]
        species = (first, second)
        agents = tuple(agent for profile in species for agent in profile.members)
        adapter = make_adapter(species, environment_volume_l=1.0e-3)
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(adapter.config.substrate_layers)
        coordinates = ((0, 0), (1, 0), (2, 0))
        initial_concentrations = np.asarray(
            ((0.05, 0.04, 0.03), (0.06, 0.05, 0.04), (0.07, 0.06, 0.05)),
            dtype=np.float64,
        )
        for coordinate, concentrations in zip(coordinates, initial_concentrations):
            ecm.deposit_substrates_m(coordinate, concentrations.copy())
        adapter.state_store.reconcile(agents)
        states_before = adapter.state_store.snapshot()
        captured = {}
        original_advance = adapter.orchestrator.advance_batch

        def capture_solver_outputs(**arguments):
            summary = original_advance(**arguments)
            for name in (
                'environment_id',
                'environment_volume_l',
                'biomass_gdw',
                'biomass_delta_gdw',
                'substrate_consumed_mol',
                'acetate_secreted_mol',
            ):
                captured[name] = arguments[name].copy()
            return summary

        adapter.orchestrator.advance_batch = capture_solver_outputs
        try:
            result = adapter.advance_tick(
                species=species,
                ecm=ecm,
                tick_index=0,
                dt_seconds=0.01,
            )
        finally:
            del adapter.orchestrator.advance_batch

        self.assertEqual(result.disposition, FbaTickDisposition.COMMITTED)
        final_concentrations = ecm.sample_substrate_batch_m(coordinates)
        for environment in range(len(coordinates)):
            members = captured['environment_id'] == environment
            for substrate in range(len(SUBSTRATE_NAMES)):
                depleted_mol = (
                    initial_concentrations[environment, substrate]
                    - final_concentrations[environment, substrate]
                ) * captured['environment_volume_l'][environment]
                consumed_mol = float(np.sum(
                    captured['substrate_consumed_mol'][members, substrate]
                ))
                allowed_mol = (
                    adapter.config.conservation_atol_m
                    * captured['environment_volume_l'][environment]
                )
                self.assertLessEqual(
                    abs(depleted_mol - consumed_mol),
                    allowed_mol,
                    msg=f'environment={environment}, substrate={substrate}',
                )

        states_after = adapter.state_store.snapshot()
        for agent_index, agent in enumerate(agents):
            before = states_before[agent]
            after = states_after[agent]
            self.assertAlmostEqual(
                after.biomass_gdw - before.biomass_gdw,
                after.biomass_delta_gdw,
                delta=BIOMASS_ASSERT_ATOL_GDW,
            )
            self.assertAlmostEqual(
                after.biomass_gdw,
                captured['biomass_gdw'][agent_index],
                delta=BIOMASS_ASSERT_ATOL_GDW,
            )
            self.assertEqual(
                after.substrate_consumed_mol,
                tuple(captured['substrate_consumed_mol'][agent_index]),
            )
            self.assertEqual(
                after.acetate_secreted_mol,
                captured['acetate_secreted_mol'][agent_index],
            )
            self.assertTrue(math.isfinite(after.acetate_secreted_mol))
            self.assertGreaterEqual(after.acetate_secreted_mol, 0.0)

        self.assertAlmostEqual(
            result.biomass_delta_gdw,
            float(np.sum(captured['biomass_delta_gdw'])),
            delta=BIOMASS_ASSERT_ATOL_GDW,
        )
        self.assertAlmostEqual(
            result.acetate_secreted_mol,
            float(np.sum(captured['acetate_secreted_mol'])),
            delta=BIOMASS_ASSERT_ATOL_GDW,
        )
        self.assertEqual(ecm.acetate, {})

    def test_inconsistent_biomass_delta_is_rejected_before_commit(self):
        species, agents = make_synthetic_community(1)
        adapter = make_adapter(species)
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(adapter.config.substrate_layers)
        ecm.deposit_substrates_m(
            (0, 0),
            np.asarray((0.1, 0.1, 0.1), dtype=np.float64),
        )
        adapter.state_store.reconcile(agents)
        substrate_before = ecm.snapshot_substrate_layers()
        states_before = adapter.state_store.snapshot()

        def inconsistent_biomass(**arguments):
            arguments['biomass_gdw'] += 1.0e-6
            return (1, 0, 0, 0, 0, 0.01)

        adapter.orchestrator.advance_batch = inconsistent_biomass
        with self.assertRaisesRegex(FbaConservationError, 'biomass'):
            adapter.advance_tick(
                species=species,
                ecm=ecm,
                tick_index=1,
                dt_seconds=0.01,
            )
        self.assertEqual(ecm.snapshot_substrate_layers(), substrate_before)
        self.assertEqual(adapter.state_store.snapshot(), states_before)

    def test_default_tolerance_rejects_ordinary_scale_unreported_depletion(self):
        species, _ = make_synthetic_community(1)
        adapter = make_adapter(species, environment_volume_l=1.0e-3)
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(adapter.config.substrate_layers)
        ecm.deposit_substrates_m(
            (0, 0),
            np.asarray((0.1, 0.1, 0.1), dtype=np.float64),
        )
        before = ecm.snapshot_substrate_layers()

        def omit_consumption(**arguments):
            # 5e-11 M in a 1e-3 L environment is 5e-14 mol of unreported loss.
            arguments['concentrations_m'][0, 0] -= 5.0e-11
            return (1, 0, 0, 0, 0, 0.01)

        adapter.orchestrator.advance_batch = omit_consumption
        with self.assertRaises(FbaConservationError):
            adapter.advance_tick(
                species=species,
                ecm=ecm,
                tick_index=0,
                dt_seconds=0.01,
            )
        self.assertEqual(ecm.snapshot_substrate_layers(), before)


class TestSimulationTransactionStress(unittest.TestCase):
    def assert_transaction_restored(self, simulation, species, agent, expected):
        self.assertEqual(simulation.ecm.snapshot_substrate_layers(), expected['substrates'])
        self.assertEqual(
            tuple(
                (name, tuple(sorted(layer.items())))
                for name, layer in simulation.ecm._layers.items()
            ),
            expected['signals'],
        )
        self.assertEqual(
            simulation.metabolic_adapter.state_store.snapshot(),
            expected['states'],
        )
        self.assertEqual(species.__dict__, expected['species'])
        self.assertEqual(agent.__dict__, expected['agent'])
        self.assertEqual(simulation.stats, expected['stats'])
        self.assertEqual(simulation.species_step_stats, expected['species_step_stats'])
        self.assertEqual(simulation.env, expected['environment'])
        self.assertEqual(
            tuple(simulation._pending_diet_boluses),
            expected['pending_boluses'],
        )
        self.assertIs(simulation.last_fba_result, expected['last_fba_result'])
        self.assertEqual(simulation.step_idx, expected['step_idx'])
        self.assertEqual(random.getstate(), expected['random_state'])

    def test_each_injected_phase_propagates_and_restores_owner_state(self):
        for phase in (
            'translator',
            'bolus deposit',
            'diffusion',
            'solver',
            'validate',
            'grid commit',
            'state-store commit',
        ):
            with self.subTest(phase=phase):
                simulation, species, agent, _ = make_live_fixture()
                expected = snapshot_live_simulation(simulation, species, agent)
                message = f'injected {phase} failure'
                adapter = simulation.metabolic_adapter

                if phase == 'translator':
                    patcher = mock.patch.object(
                        adapter.diet_engine,
                        'get_molar_substrate_bounds',
                        side_effect=InjectedPhaseError(message),
                    )
                elif phase == 'bolus deposit':
                    original_replace = simulation.ecm.replace_substrate_batch_m

                    def fail_after_bolus_commit(coordinates, concentrations_m):
                        original_replace(coordinates, concentrations_m)
                        raise InjectedPhaseError(message)

                    patcher = mock.patch.object(
                        simulation.ecm,
                        'replace_substrate_batch_m',
                        side_effect=fail_after_bolus_commit,
                    )
                elif phase == 'diffusion':
                    def fail_diffusion(dt_seconds):
                        for layer in simulation.ecm._substrate_layers.values():
                            layer[(91, -17)] = 9.25
                        raise InjectedPhaseError(message)

                    patcher = mock.patch.object(
                        simulation.ecm,
                        'tick_substrates',
                        side_effect=fail_diffusion,
                    )
                elif phase == 'solver':
                    def fail_solver(**arguments):
                        arguments['concentrations_m'].fill(0.0)
                        arguments['biomass_gdw'].fill(99.0)
                        raise InjectedPhaseError(message)

                    patcher = mock.patch.object(
                        adapter.orchestrator,
                        'advance_batch',
                        side_effect=fail_solver,
                    )
                elif phase == 'validate':
                    patcher = mock.patch.object(
                        adapter,
                        '_validate_batch',
                        side_effect=InjectedPhaseError(message),
                    )
                elif phase == 'grid commit':
                    original_replace = simulation.ecm.replace_substrate_batch_m

                    def fail_after_grid_commit(coordinates, concentrations_m):
                        original_replace(coordinates, concentrations_m)
                        raise InjectedPhaseError(message)

                    patcher = mock.patch.object(
                        simulation.ecm,
                        'replace_substrate_batch_m',
                        side_effect=fail_after_grid_commit,
                    )
                else:
                    original_commit = adapter.state_store.commit_batch

                    def fail_after_state_commit(batch):
                        original_commit(batch)
                        raise InjectedPhaseError(message)

                    patcher = mock.patch.object(
                        adapter.state_store,
                        'commit_batch',
                        side_effect=fail_after_state_commit,
                    )

                with patcher:
                    with self.assertRaises(InjectedPhaseError) as raised:
                        simulation.update(100, 100, 0.01)
                self.assertIs(type(raised.exception), InjectedPhaseError)
                self.assertEqual(str(raised.exception), message)
                self.assert_transaction_restored(
                    simulation,
                    species,
                    agent,
                    expected,
                )

    def test_post_commit_state_store_failure_restores_complete_simulation_state(self):
        simulation, species, agent, bolus = make_live_fixture()
        expected = snapshot_live_simulation(simulation, species, agent)
        original_state_for = simulation.metabolic_adapter.state_store.state_for
        call_count = 0

        def fail_after_adapter_commit(target):
            nonlocal call_count
            call_count += 1
            if call_count == 2:
                raise InjectedPhaseError('injected post-commit state-store failure')
            return original_state_for(target)

        with mock.patch.object(
            simulation.metabolic_adapter.state_store,
            'state_for',
            side_effect=fail_after_adapter_commit,
        ):
            with self.assertRaises(InjectedPhaseError) as raised:
                simulation.update(100, 100, 0.01)
        self.assertIs(type(raised.exception), InjectedPhaseError)
        self.assertEqual(
            str(raised.exception),
            'injected post-commit state-store failure',
        )
        self.assert_transaction_restored(simulation, species, agent, expected)
        self.assertEqual(simulation._pending_diet_boluses, [bolus])
        self.assertIs(species.members[0], agent)

    def test_late_statistics_failure_restores_lifecycle_and_prior_appends(self):
        simulation, species, agent, _ = make_live_fixture()
        events = []
        simulation.metabolic_adapter.event_sink = events.append
        simulation.stats['fba_biomass_gdw'] = AppendThenRaiseList()
        expected = snapshot_live_simulation(simulation, species, agent)

        with self.assertRaises(InjectedPhaseError) as raised:
            simulation.update(100, 100, 0.01)

        self.assertIs(type(raised.exception), InjectedPhaseError)
        self.assertEqual(str(raised.exception), 'injected statistics commit failure')
        self.assert_transaction_restored(simulation, species, agent, expected)
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].disposition, FbaTickDisposition.ROLLED_BACK)
        self.assertEqual(events[0].phase, 'simulation')
        self.assertEqual(
            simulation.metabolic_adapter.orchestrator._kernel.cache_stats(),
            (0, 0, 0, 0, 0),
        )

    def test_late_failure_restores_preexisting_tolerant_cache_exactly(self):
        simulation, species, agent, _ = make_live_fixture()
        adapter = simulation.metabolic_adapter
        events = []
        adapter.event_sink = events.append
        concentrations = np.asarray((0.05, 0.04, 0.03), dtype=np.float64)
        warm_flux, warm_hit = adapter.orchestrator.solve_flux(
            concentrations_m=concentrations,
            cache_enabled=True,
        )
        self.assertFalse(warm_hit)
        starting_stats = adapter.orchestrator._kernel.cache_stats()
        simulation.stats['fba_biomass_gdw'] = AppendThenRaiseList()

        with self.assertRaises(InjectedPhaseError):
            simulation.update(100, 100, 0.01)

        self.assertEqual(adapter.orchestrator._kernel.cache_stats(), starting_stats)
        restored_flux, restored_hit = adapter.orchestrator.solve_flux(
            concentrations_m=concentrations,
            cache_enabled=True,
        )
        self.assertTrue(restored_hit)
        np.testing.assert_array_equal(restored_flux, warm_flux)
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].disposition, FbaTickDisposition.ROLLED_BACK)

    def test_extreme_glycan_depletion_and_nonfinite_outputs_roll_back_exactly(self):
        for case in ('extreme glycan depletion', 'nan', 'positive infinity'):
            with self.subTest(case=case):
                simulation, species, agent, _ = make_live_fixture()
                if case == 'extreme glycan depletion':
                    simulation.ecm.deposit_substrates_m(
                        (0, 0),
                        np.asarray((0.0, 1.0e200, 0.0), dtype=np.float64),
                    )
                expected = snapshot_live_simulation(simulation, species, agent)

                def invalid_solver_output(**arguments):
                    if case == 'extreme glycan depletion':
                        arguments['concentrations_m'][:, 1] = 0.0
                    elif case == 'nan':
                        arguments['concentrations_m'][0, 0] = np.nan
                    else:
                        arguments['concentrations_m'][0, 0] = np.inf
                    return (1, 0, 0, 0, 0, 0.01)

                with mock.patch.object(
                    simulation.metabolic_adapter.orchestrator,
                    'advance_batch',
                    side_effect=invalid_solver_output,
                ):
                    expected_exception = (
                        FbaConservationError
                        if case == 'extreme glycan depletion'
                        else ValueError
                    )
                    with self.assertRaises(expected_exception):
                        simulation.update(100, 100, 0.01)
                self.assert_transaction_restored(
                    simulation,
                    species,
                    agent,
                    expected,
                )


class TestLargeDeterministicBatch(unittest.TestCase):
    def test_5001_shared_cell_agents_keep_order_dtypes_and_rollback_exactly(self):
        agent_count = 5_001
        species, agents = make_synthetic_community(
            agent_count,
            species_count=2,
            shared_cell=True,
        )
        adapter = make_adapter(
            species,
            environment_volume_l=1.0e-3,
            initial_biomass_gdw=1.0e-12,
        )
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(adapter.config.substrate_layers)
        ecm.deposit_substrates_m(
            (0, 0),
            np.asarray((0.1, 0.1, 0.1), dtype=np.float64),
        )

        first = build_coordinate_environment_batch(
            species,
            ecm,
            adapter.state_store,
            adapter.phenotype_index,
            adapter.config.environment_volume_l,
        )
        second = build_coordinate_environment_batch(
            species,
            ecm,
            adapter.state_store,
            adapter.phenotype_index,
            adapter.config.environment_volume_l,
        )
        self.assertEqual(first.agents, agents)
        self.assertEqual(second.agents, first.agents)
        self.assertEqual(first.coordinates, ((0, 0),))
        self.assertEqual(first.phenotype_id.shape, (agent_count,))
        self.assertEqual(first.environment_id.shape, (agent_count,))
        self.assertEqual(first.concentrations_m.shape, (1, 3))
        self.assertEqual(first.substrate_consumed_mol.shape, (agent_count, 3))
        self.assertEqual(first.phenotype_id.dtype, np.dtype(np.int32))
        self.assertEqual(first.environment_id.dtype, np.dtype(np.int32))
        for array in (
            first.concentrations_m,
            first.environment_volume_l,
            first.biomass_gdw,
            first.growth_rate_s,
            first.biomass_delta_gdw,
            first.substrate_consumed_mol,
            first.acetate_secreted_mol,
        ):
            self.assertEqual(array.dtype, np.dtype(np.float64))
            self.assertTrue(array.flags.c_contiguous)
        np.testing.assert_array_equal(second.phenotype_id, first.phenotype_id)
        np.testing.assert_array_equal(second.environment_id, first.environment_id)

        observed = {}
        original_advance = adapter.orchestrator.advance_batch

        def inspect_integer_inputs(**arguments):
            phenotype_before = arguments['phenotype_id'].copy()
            environment_before = arguments['environment_id'].copy()
            summary = original_advance(**arguments)
            observed['phenotype'] = arguments['phenotype_id'].copy()
            observed['environment'] = arguments['environment_id'].copy()
            observed['phenotype_dtype'] = arguments['phenotype_id'].dtype
            observed['environment_dtype'] = arguments['environment_id'].dtype
            np.testing.assert_array_equal(arguments['phenotype_id'], phenotype_before)
            np.testing.assert_array_equal(arguments['environment_id'], environment_before)
            return summary

        adapter.orchestrator.advance_batch = inspect_integer_inputs
        try:
            result = adapter.advance_tick(
                species=species,
                ecm=ecm,
                tick_index=0,
                dt_seconds=1.0e-4,
            )
        finally:
            del adapter.orchestrator.advance_batch
        self.assertEqual(result.disposition, FbaTickDisposition.COMMITTED)
        self.assertEqual(result.agent_count, agent_count)
        self.assertEqual(result.environment_count, 1)
        self.assertEqual(observed['phenotype_dtype'], np.dtype(np.int32))
        self.assertEqual(observed['environment_dtype'], np.dtype(np.int32))
        self.assertEqual(observed['phenotype'].shape, (agent_count,))
        self.assertEqual(observed['environment'].shape, (agent_count,))
        last_state = adapter.state_store.state_for(agents[-1])
        self.assertGreater(last_state.biomass_delta_gdw, 0.0)
        self.assertGreater(sum(last_state.substrate_consumed_mol), 0.0)
        self.assertGreater(last_state.acetate_secreted_mol, 0.0)
        self.assertGreater(result.biomass_delta_gdw, 0.0)
        self.assertTrue(all(value > 0.0 for value in result.substrate_consumed_mol))

        substrate_before_failure = ecm.snapshot_substrate_layers()
        states_before_failure = adapter.state_store.snapshot()

        def fail_large_solver(**arguments):
            arguments['concentrations_m'].fill(0.0)
            arguments['biomass_gdw'].fill(1.0)
            arguments['substrate_consumed_mol'].fill(2.0)
            raise InjectedPhaseError('injected 5001-agent solver failure')

        with mock.patch.object(
            adapter.orchestrator,
            'advance_batch',
            side_effect=fail_large_solver,
        ):
            with self.assertRaisesRegex(
                InjectedPhaseError,
                '^injected 5001-agent solver failure$',
            ):
                adapter.advance_tick(
                    species=species,
                    ecm=ecm,
                    tick_index=1,
                    dt_seconds=1.0e-4,
                )
        self.assertEqual(ecm.snapshot_substrate_layers(), substrate_before_failure)
        self.assertEqual(adapter.state_store.snapshot(), states_before_failure)
        self.assertEqual(
            tuple(agent for profile in species for agent in profile.members),
            agents,
        )


class TestIntegrationObjectOwnership(unittest.TestCase):
    """Python reachability checks; native leak evidence remains the Valgrind job."""

    def test_repeated_batches_results_and_events_are_collectible(self):
        species, _ = make_synthetic_community(3)
        result_references = []
        event_references = []
        array_references = []

        def observe_event(event):
            event_references.append(weakref.ref(event))

        adapter = make_adapter(species, event_sink=observe_event)
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(adapter.config.substrate_layers)
        original_advance = adapter.orchestrator.advance_batch

        def observe_arrays(**arguments):
            array_references.extend(
                weakref.ref(arguments[name])
                for name in (
                    'phenotype_id',
                    'environment_id',
                    'concentrations_m',
                    'environment_volume_l',
                    'biomass_gdw',
                    'growth_rate_s',
                    'biomass_delta_gdw',
                    'substrate_consumed_mol',
                    'acetate_secreted_mol',
                )
            )
            return original_advance(**arguments)

        adapter.orchestrator.advance_batch = observe_arrays
        try:
            for tick_index in range(40):
                result = adapter.advance_tick(
                    species=species,
                    ecm=ecm,
                    tick_index=tick_index,
                    dt_seconds=0.0,
                )
                result_references.append(weakref.ref(result))
                del result
        finally:
            del adapter.orchestrator.advance_batch

        gc.collect()
        self.assertTrue(result_references)
        self.assertTrue(event_references)
        self.assertTrue(array_references)
        self.assertTrue(all(reference() is None for reference in result_references))
        self.assertTrue(all(reference() is None for reference in event_references))
        self.assertTrue(all(reference() is None for reference in array_references))

    def test_python_transaction_retention_is_bounded_after_warmup(self):
        species, _ = make_synthetic_community(8)
        quiet_logger = logging.getLogger('bacsandbox.integration-retention')
        quiet_logger.disabled = True
        adapter = make_adapter(species, logger=quiet_logger)
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(adapter.config.substrate_layers)

        def run_tick(tick_index):
            adapter.advance_tick(
                species=species,
                ecm=ecm,
                tick_index=tick_index,
                dt_seconds=0.0,
            )

        for tick_index in range(30):
            run_tick(tick_index)
        gc.collect()
        tracemalloc.start()
        try:
            for tick_index in range(30, 60):
                run_tick(tick_index)
            gc.collect()
            baseline_bytes = tracemalloc.get_traced_memory()[0]
            for tick_index in range(60, 310):
                run_tick(tick_index)
            gc.collect()
            retained_growth = tracemalloc.get_traced_memory()[0] - baseline_bytes
            self.assertLessEqual(
                retained_growth,
                RETENTION_LIMIT_BYTES,
                msg=f'Python integration retention grew by {retained_growth} bytes',
            )
        finally:
            tracemalloc.stop()


class TestBackendCacheSnapshots(unittest.TestCase):
    def test_python_and_native_cache_snapshots_restore_payload_and_counters(self):
        backends = ['python']
        if NATIVE_METABOLIC_AVAILABLE:
            backends.append('cpp')
        for backend in backends:
            with self.subTest(backend=backend):
                species, _ = make_synthetic_community(1)
                adapter = make_adapter(species, backend=backend)
                concentrations = np.asarray((0.05, 0.04, 0.03), dtype=np.float64)
                first_flux, first_hit = adapter.orchestrator.solve_flux(
                    concentrations_m=concentrations,
                    cache_enabled=True,
                )
                self.assertFalse(first_hit)
                snapshot = adapter.orchestrator.snapshot_cache()
                stats = adapter.orchestrator._kernel.cache_stats()

                adapter.orchestrator.solve_flux(
                    concentrations_m=concentrations,
                    cache_enabled=True,
                )
                self.assertNotEqual(adapter.orchestrator._kernel.cache_stats(), stats)
                adapter.orchestrator.restore_cache(snapshot)
                self.assertEqual(adapter.orchestrator._kernel.cache_stats(), stats)

                restored_flux, restored_hit = adapter.orchestrator.solve_flux(
                    concentrations_m=concentrations,
                    cache_enabled=True,
                )
                self.assertTrue(restored_hit)
                np.testing.assert_array_equal(restored_flux, first_flux)


@unittest.skipUnless(NATIVE_METABOLIC_AVAILABLE, 'native MetabolicKernel is not built')
class TestNativeIntegrationParity(unittest.TestCase):
    def run_backend(self, backend):
        species, agents = make_synthetic_community(
            24,
            species_count=2,
            shared_cell=False,
        )
        for index, agent in enumerate(agents):
            agent.x = float((index % 4) * 10 + 1)
            agent.y = 1.0
        adapter = make_adapter(
            species,
            backend=backend,
            environment_volume_l=1.0e-3,
            initial_biomass_gdw=1.0e-10,
        )
        ecm = ECM(cell_size=10.0, backend='python')
        ecm.register_substrate_layers(adapter.config.substrate_layers)
        coordinates = ((0, 0), (1, 0), (2, 0), (3, 0))
        for index, coordinate in enumerate(coordinates):
            ecm.deposit_substrates_m(
                coordinate,
                np.asarray(
                    (0.05 + index * 0.01, 0.04 + index * 0.01, 0.03 + index * 0.01),
                    dtype=np.float64,
                ),
            )
        result = adapter.advance_tick(
            species=species,
            ecm=ecm,
            tick_index=0,
            dt_seconds=0.005,
        )
        states = adapter.state_store.snapshot()
        state_matrix = np.asarray(
            [
                (
                    states[agent].biomass_gdw,
                    states[agent].growth_rate_s,
                    states[agent].biomass_delta_gdw,
                    *states[agent].substrate_consumed_mol,
                    states[agent].acetate_secreted_mol,
                )
                for agent in agents
            ],
            dtype=np.float64,
        )
        return result, ecm.sample_substrate_batch_m(coordinates), state_matrix

    def test_python_and_native_adapter_outputs_match(self):
        python_result, python_concentrations, python_states = self.run_backend('python')
        native_result, native_concentrations, native_states = self.run_backend('cpp')

        self.assertEqual(python_result.backend_used, 'python')
        self.assertEqual(native_result.backend_used, 'cpp')
        self.assertEqual(native_result.kernel_summary[:5], python_result.kernel_summary[:5])
        self.assertAlmostEqual(
            native_result.kernel_summary[5],
            python_result.kernel_summary[5],
            delta=1.0e-12,
        )
        np.testing.assert_allclose(
            native_concentrations,
            python_concentrations,
            rtol=0.0,
            atol=1.0e-15,
        )
        self.assertTrue(np.all(python_states[:, 2:] > 0.0))
        self.assertTrue(np.all(native_states[:, 2:] > 0.0))
        np.testing.assert_allclose(
            native_states[:, :2],
            python_states[:, :2],
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            native_states[:, 2:],
            python_states[:, 2:],
            rtol=0.0,
            atol=1.0e-24,
        )


if __name__ == '__main__':
    unittest.main()