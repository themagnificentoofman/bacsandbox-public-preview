import copy
import importlib
from importlib import resources
import json
import math
from pathlib import Path
import tempfile
import time
import unittest

import numpy as np

from bacteria.diet_translator import (
    DEFAULT_DIET_PROFILES_PATH,
    DietProfileValidationError,
    DietTranslationEngine,
)
from bacteria.metabolism import (
    MetabolicOrchestrator,
    SPECIES_BASE_UPTAKE_VMAX_MOL_GDW_S,
    SPECIES_UPTAKE_KM_MOLAR,
    build_species_phenotype_config,
    load_metabolic_model,
    load_metabolic_priors,
)


try:
    bacteria_c_kernel = importlib.import_module('bacteria_c_kernel')
except ImportError:
    bacteria_c_kernel = None


NATIVE_METABOLIC_AVAILABLE = (
    bacteria_c_kernel is not None and hasattr(bacteria_c_kernel, 'MetabolicKernel')
)
EXPECTED_DIET_ORDER = ('malnutrition', 'rusf', 'mdcf_2')
EXPECTED_SUBSTRATE_ORDER = ('generic_calorie', 'complex_glycan', 'mdcf_prebiotic')
EXPECTED_COMPONENT_ORDER = (
    'simple_hexose',
    'digestible_starch',
    'resistant_starch',
    'cellulose',
    'arabinan',
    'xylan',
    'beta_1_4_galactan',
    'beta_1_4_mannan',
)
# These independently frozen 50 g/1 L concentrations lock the documented
# glucose-equivalent bolus calibration; see DOCUMENTATION.md section 10.4.
EXPECTED_MDCF_50G_COMPONENTS_M = np.asarray((
    0.08275575261007312,
    0.03478716108750197,
    0.0,
    0.001064944097366153,
    0.00010138550371321273,
    0.00006039452350327761,
    0.0002780913974742124,
    0.000020245970462459088,
), dtype=np.float64)
EXPECTED_MDCF_50G_SUBSTRATES_M = np.asarray((
    0.11754291369757508,
    0.0012267241245826435,
    0.00029833736793667145,
), dtype=np.float64)


class DietTranslatorAssertions:
    def assert_owned_float64_array(self, value, shape):
        self.assertIsInstance(value, np.ndarray)
        self.assertEqual(value.shape, shape)
        self.assertEqual(value.dtype, np.dtype(np.float64))
        self.assertTrue(value.flags.c_contiguous)
        self.assertTrue(value.flags.writeable)
        self.assertTrue(value.flags.owndata)
        self.assertTrue(np.all(np.isfinite(value)))
        self.assertTrue(np.all(value >= 0.0))


class TestDietTranslationContracts(DietTranslatorAssertions, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = DietTranslationEngine()

    def test_production_asset_loads_with_exact_public_ordering(self):
        self.assertTrue(DEFAULT_DIET_PROFILES_PATH.is_file())
        self.assertEqual(self.engine.diet_names, EXPECTED_DIET_ORDER)
        self.assertEqual(self.engine.substrate_order, EXPECTED_SUBSTRATE_ORDER)
        self.assertEqual(self.engine.component_order, EXPECTED_COMPONENT_ORDER)

        packaged_asset = resources.files('bacteria').joinpath('models', 'diet_profiles.json')
        self.assertTrue(packaged_asset.is_file())

    def test_mdcf_50g_formula_outputs_are_independently_frozen(self):
        components = self.engine.get_molar_component_concentrations(
            'mdcf_2', 50.0, 1.0, 1.0
        )
        bounds = self.engine.get_molar_substrate_bounds('mdcf_2', 50.0, 1.0, 1.0)

        galactan_index = self.engine.component_order.index('beta_1_4_galactan')
        mannan_index = self.engine.component_order.index('beta_1_4_mannan')
        resistant_starch_index = self.engine.component_order.index('resistant_starch')
        self.assertEqual(components[resistant_starch_index], 0.0)
        self.assertAlmostEqual(
            components[galactan_index],
            0.0002780913974742124,
            delta=1.0e-15,
        )
        self.assertAlmostEqual(
            components[mannan_index],
            0.000020245970462459088,
            delta=1.0e-15,
        )
        combined_prebiotic = components[galactan_index] + components[mannan_index]
        self.assertLessEqual(
            abs(combined_prebiotic - 0.00029833736793667145),
            1.0e-6,
        )
        self.assertAlmostEqual(
            combined_prebiotic,
            0.00029833736793667145,
            delta=1.0e-15,
        )

        for component, expected in zip(components, EXPECTED_MDCF_50G_COMPONENTS_M):
            self.assertAlmostEqual(component, expected, delta=1.0e-15)
        for substrate, expected in zip(bounds, EXPECTED_MDCF_50G_SUBSTRATES_M):
            self.assertAlmostEqual(substrate, expected, delta=1.0e-15)

    def test_mass_volume_zero_and_timestep_bolus_semantics(self):
        baseline = self.engine.get_molar_substrate_bounds('mdcf_2', 50.0, 1.0, 1.0)
        double_mass = self.engine.get_molar_substrate_bounds('mdcf_2', 100.0, 1.0, 1.0)
        double_volume = self.engine.get_molar_substrate_bounds('mdcf_2', 50.0, 2.0, 1.0)
        # Linear feed-mass and inverse-volume scaling conserve delivered component
        # amount under bolus semantics; see DOCUMENTATION.md section 10.4.
        np.testing.assert_allclose(double_mass, 2.0 * baseline, rtol=0.0, atol=1.0e-15)
        np.testing.assert_allclose(double_volume, 0.5 * baseline, rtol=0.0, atol=1.0e-15)

        zero_bounds = self.engine.get_molar_substrate_bounds('mdcf_2', 0.0, 1.0, 1.0)
        zero_components = self.engine.get_molar_component_concentrations(
            'mdcf_2', 0.0, 1.0, 1.0
        )
        np.testing.assert_array_equal(zero_bounds, np.zeros(3, dtype=np.float64))
        np.testing.assert_array_equal(
            zero_components,
            np.zeros(len(EXPECTED_COMPONENT_ORDER), dtype=np.float64),
        )

        for dt_s in (0.01, 0.1, 10.0, 3600.0):
            with self.subTest(dt_s=dt_s):
                actual = self.engine.get_molar_substrate_bounds(
                    'mdcf_2', 50.0, 1.0, dt_s
                )
                np.testing.assert_array_equal(actual, baseline)

    def test_return_arrays_are_fresh_owned_writable_and_finite(self):
        first_bounds = self.engine.get_molar_substrate_bounds('mdcf_2', 50.0, 1.0, 1.0)
        second_bounds = self.engine.get_molar_substrate_bounds('mdcf_2', 50.0, 1.0, 1.0)
        first_components = self.engine.get_molar_component_concentrations(
            'mdcf_2', 50.0, 1.0, 1.0
        )
        second_components = self.engine.get_molar_component_concentrations(
            'mdcf_2', 50.0, 1.0, 1.0
        )

        self.assert_owned_float64_array(first_bounds, (3,))
        self.assert_owned_float64_array(second_bounds, (3,))
        self.assert_owned_float64_array(first_components, (len(EXPECTED_COMPONENT_ORDER),))
        self.assert_owned_float64_array(second_components, (len(EXPECTED_COMPONENT_ORDER),))
        self.assertFalse(np.shares_memory(first_bounds, second_bounds))
        self.assertFalse(np.shares_memory(first_components, second_components))
        self.assertFalse(np.shares_memory(first_bounds, first_components))

        expected_second = second_bounds.copy()
        first_bounds.fill(0.0)
        np.testing.assert_array_equal(second_bounds, expected_second)

    def test_request_validation_rejects_bad_names_types_and_ranges(self):
        for diet_name, exception_type in (
            (None, TypeError),
            (True, TypeError),
            (1, TypeError),
            (b'mdcf_2', TypeError),
            ('', ValueError),
            ('MDCF-2', ValueError),
            ('unknown', ValueError),
        ):
            with self.subTest(field='diet_name', value=diet_name):
                with self.assertRaises(exception_type):
                    self.engine.get_molar_substrate_bounds(diet_name, 50.0, 1.0, 1.0)

        for field in ('feed_mass_g', 'volume_L', 'dt_s'):
            for value in (True, False, '1', None, 1.0 + 0.0j, object()):
                arguments = {
                    'diet_name': 'mdcf_2',
                    'feed_mass_g': 50.0,
                    'volume_L': 1.0,
                    'dt_s': 1.0,
                }
                arguments[field] = value
                with self.subTest(field=field, value=repr(value)):
                    with self.assertRaises(TypeError):
                        self.engine.get_molar_substrate_bounds(**arguments)

        for value in (-1.0, float('nan'), float('inf'), float('-inf')):
            with self.subTest(field='feed_mass_g', value=value):
                with self.assertRaises(ValueError):
                    self.engine.get_molar_substrate_bounds('mdcf_2', value, 1.0, 1.0)

        for field in ('volume_L', 'dt_s'):
            for value in (0.0, -1.0, float('nan'), float('inf'), float('-inf')):
                arguments = {
                    'diet_name': 'mdcf_2',
                    'feed_mass_g': 50.0,
                    'volume_L': 1.0,
                    'dt_s': 1.0,
                }
                arguments[field] = value
                with self.subTest(field=field, value=value):
                    with self.assertRaises(ValueError):
                        self.engine.get_molar_substrate_bounds(**arguments)

    def test_extreme_finite_request_fails_cleanly_on_overflow(self):
        # The largest finite mass divided by the smallest positive volume probes
        # representability without accepting NaN/Inf output; see DOCUMENTATION.md section 10.4.
        with self.assertRaises(OverflowError):
            self.engine.get_molar_substrate_bounds(
                'mdcf_2',
                np.finfo(np.float64).max,
                np.finfo(np.float64).tiny,
                1.0,
            )


class TestDietProfileAssetValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.valid_payload = json.loads(DEFAULT_DIET_PROFILES_PATH.read_text(encoding='utf-8'))

    @staticmethod
    def set_path(payload, path, value):
        target = payload
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value

    @staticmethod
    def delete_path(payload, path):
        target = payload
        for key in path[:-1]:
            target = target[key]
        del target[path[-1]]

    def assert_payload_rejected(self, mutation, message_pattern=None):
        payload = copy.deepcopy(self.valid_payload)
        mutation(payload)
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / 'diet_profiles.json'
            path.write_text(json.dumps(payload, allow_nan=True), encoding='utf-8')
            if message_pattern is None:
                context = self.assertRaises(DietProfileValidationError)
            else:
                context = self.assertRaisesRegex(
                    DietProfileValidationError,
                    message_pattern,
                )
            with context:
                DietTranslationEngine(path)

    def assert_raw_asset_rejected(self, raw_text, message_pattern=None):
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / 'diet_profiles.json'
            path.write_text(raw_text, encoding='utf-8')
            if message_pattern is None:
                context = self.assertRaises(DietProfileValidationError)
            else:
                context = self.assertRaisesRegex(
                    DietProfileValidationError,
                    message_pattern,
                )
            with context:
                DietTranslationEngine(path)

    def test_invalid_json_and_duplicate_keys_at_every_depth_are_rejected(self):
        self.assert_raw_asset_rejected('{ invalid JSON', 'cannot load')
        duplicate_assets = (
            '{"schema_version": 1, "schema_version": 1}',
            '{"outer": {"source_hook": "first", "source_hook": "second"}}',
        )
        for raw_text in duplicate_assets:
            with self.subTest(raw_text=raw_text):
                self.assert_raw_asset_rejected(raw_text, 'duplicate JSON key')

    def test_missing_and_unknown_fields_are_rejected(self):
        cases = (
            (
                'missing root field',
                lambda payload: self.delete_path(payload, ('scope',)),
            ),
            (
                'unknown root field',
                lambda payload: self.set_path(payload, ('unexpected',), True),
            ),
            (
                'missing nested field',
                lambda payload: self.delete_path(
                    payload,
                    ('components', 'simple_hexose', 'display_name'),
                ),
            ),
            (
                'unknown nested field',
                lambda payload: self.set_path(
                    payload,
                    ('components', 'simple_hexose', 'unexpected'),
                    'value',
                ),
            ),
        )
        for name, mutation in cases:
            with self.subTest(case=name):
                self.assert_payload_rejected(mutation, 'missing keys|unknown keys')

    def test_wrong_schema_identity_version_and_orders_are_rejected(self):
        cases = (
            (
                'schema version',
                lambda payload: self.set_path(payload, ('schema_version',), 2),
            ),
            (
                'profile revision',
                lambda payload: self.set_path(payload, ('profile_revision',), 2),
            ),
            (
                'profile id',
                lambda payload: self.set_path(payload, ('profile_id',), 'other_profile'),
            ),
            (
                'substrate order',
                lambda payload: payload['substrate_order'].reverse(),
            ),
            (
                'component order',
                lambda payload: payload['component_order'].reverse(),
            ),
            (
                'ingredient order',
                lambda payload: payload['ingredient_order'].reverse(),
            ),
            (
                'diet order',
                lambda payload: payload['diet_order'].reverse(),
            ),
        )
        for name, mutation in cases:
            with self.subTest(case=name):
                self.assert_payload_rejected(mutation)

    def test_unresolved_source_hooks_and_bad_cross_references_are_rejected(self):
        cases = (
            (
                'component source hook',
                lambda payload: self.set_path(
                    payload,
                    ('components', 'simple_hexose', 'source_hook'),
                    'missing-hook',
                ),
            ),
            (
                'ingredient source hook',
                lambda payload: self.set_path(
                    payload,
                    ('ingredients', 'chickpea', 'source_hook'),
                    'missing-hook',
                ),
            ),
            (
                'diet source hook',
                lambda payload: self.set_path(
                    payload,
                    ('diets', 'mdcf_2', 'source_hook'),
                    'missing-hook',
                ),
            ),
            (
                'component substrate',
                lambda payload: self.set_path(
                    payload,
                    ('components', 'simple_hexose', 'substrate'),
                    'missing-substrate',
                ),
            ),
        )
        for name, mutation in cases:
            with self.subTest(case=name):
                self.assert_payload_rejected(mutation, 'unknown hook|unknown substrate')

    def test_fraction_ranges_and_component_sum_are_rejected(self):
        fraction_path = (
            'ingredients',
            'chickpea',
            'component_mass_fractions',
            'simple_hexose',
        )
        for value in (-0.01, 1.01):
            with self.subTest(case='fraction range', value=value):
                self.assert_payload_rejected(
                    lambda payload, value=value: self.set_path(payload, fraction_path, value),
                    r'\[0.0, 1.0\]',
                )

        def component_sum_above_one(payload):
            fractions = payload['ingredients']['added_sugar']['component_mass_fractions']
            fractions['simple_hexose'] = 0.75
            fractions['digestible_starch'] = 0.5

        self.assert_payload_rejected(component_sum_above_one, 'sum above 1.0')

    def test_invalid_recipe_totals_and_rounding_are_rejected(self):
        def mismatched_reported_total(payload):
            payload['diets']['mdcf_2']['reported_total_g_per_100g'] = 99.0

        def outside_rounding_tolerance(payload):
            diet = payload['diets']['mdcf_2']
            diet['reported_total_g_per_100g'] = 101.0
            diet['ingredient_mass_g_per_100g']['added_sugar'] = 30.86

        def zero_dry_mass(payload):
            diet = payload['diets']['mdcf_2']
            diet['reported_total_g_per_100g'] = 0.0
            for ingredient in diet['ingredient_mass_g_per_100g']:
                diet['ingredient_mass_g_per_100g'][ingredient] = 0.0

        cases = (
            ('reported mismatch', mismatched_reported_total, 'does not match'),
            ('rounding tolerance', outside_rounding_tolerance, 'rounding tolerance'),
            ('zero dry mass', zero_dry_mass, 'positive dry mass'),
        )
        for name, mutation, message_pattern in cases:
            with self.subTest(case=name):
                self.assert_payload_rejected(mutation, message_pattern)

    def test_nonfinite_json_numbers_are_rejected_when_parser_represents_them(self):
        value_path = (
            'ingredients',
            'chickpea',
            'component_mass_fractions',
            'simple_hexose',
        )
        for value in (float('nan'), float('inf'), float('-inf')):
            with self.subTest(value=value):
                self.assert_payload_rejected(
                    lambda payload, value=value: self.set_path(payload, value_path, value),
                    'finite',
                )


class TestDietTranslatorMetabolicIntegration(DietTranslatorAssertions, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = DietTranslationEngine()
        cls.model = load_metabolic_model()
        cls.priors = load_metabolic_priors()

    def phenotype_config(self):
        # Capacity multipliers isolate one high-MDCF phenotype and one exact
        # non-utilizer without changing Km; see DOCUMENTATION.md section 10.4.
        profiles = (
            {'aff_calorie': 1.0, 'aff_glycan': 2.0, 'aff_mdcf': 4.0},
            {'aff_calorie': 1.0, 'aff_glycan': 1.0, 'aff_mdcf': 0.0},
        )
        return build_species_phenotype_config(
            profiles,
            self.model,
            self.priors,
            8,
        )

    def orchestrator(self, backend):
        return MetabolicOrchestrator(
            mode='fba',
            backend=backend,
            model=self.model,
            priors=self.priors,
            phenotype_config=self.phenotype_config(),
        )

    def translated_batch(self):
        increment = self.engine.get_molar_substrate_bounds('mdcf_2', 50.0, 1.0, 1.0)
        concentrations = np.zeros((2, 3), dtype=np.float64, order='C')
        concentrations[0] += increment
        concentrations[1] += increment
        # One gram dry biomass in each 1 L environment makes uptake measurable
        # without approaching depletion; see DOCUMENTATION.md section 10.4.
        return {
            'phenotype_id': np.asarray((0, 1), dtype=np.int32),
            'environment_id': np.asarray((0, 1), dtype=np.int32),
            'concentrations_m': concentrations,
            'environment_volume_l': np.ones(2, dtype=np.float64),
            'biomass_gdw': np.ones(2, dtype=np.float64),
            'growth_rate_s': np.zeros(2, dtype=np.float64),
            'biomass_delta_gdw': np.zeros(2, dtype=np.float64),
            'substrate_consumed_mol': np.zeros((2, 3), dtype=np.float64),
            'acetate_secreted_mol': np.zeros(2, dtype=np.float64),
        }

    def advance(self, backend):
        orchestrator = self.orchestrator(backend)
        batch = self.translated_batch()
        initial_concentrations = batch['concentrations_m'].copy()
        summary = orchestrator.advance_batch(
            **batch,
            dt_seconds=0.1,
            depletion_safety=0.9,
            max_substeps=8,
            cache_enabled=False,
        )
        return orchestrator, batch, initial_concentrations, summary

    def assert_valid_advance(self, orchestrator, batch, initial_concentrations, summary):
        self.assertEqual(batch['concentrations_m'].shape, (2, 3))
        self.assertTrue(batch['concentrations_m'].flags.c_contiguous)
        self.assertTrue(batch['concentrations_m'].flags.writeable)
        self.assertTrue(np.all(np.isfinite(batch['concentrations_m'])))
        self.assertTrue(np.all(batch['concentrations_m'] >= 0.0))
        for name in (
            'biomass_gdw',
            'growth_rate_s',
            'biomass_delta_gdw',
            'substrate_consumed_mol',
            'acetate_secreted_mol',
        ):
            self.assertTrue(np.all(np.isfinite(batch[name])), name)
        self.assertTrue(all(math.isfinite(float(value)) for value in summary))
        self.assertGreater(batch['substrate_consumed_mol'][0, 2], 0.0)
        self.assertEqual(batch['substrate_consumed_mol'][1, 2], 0.0)
        self.assertLess(
            batch['concentrations_m'][0, 2],
            initial_concentrations[0, 2],
        )
        self.assertEqual(
            batch['concentrations_m'][1, 2],
            initial_concentrations[1, 2],
        )
        self.assertIn(orchestrator.backend_used, ('python', 'cpp'))

    def test_species_affinities_map_to_exact_vmax_and_shared_km(self):
        config = self.phenotype_config()
        # The hard-coded base capacities and Km values are the Phase 3 ABI
        # calibration under test; see DOCUMENTATION.md section 10.4.
        base_vmax = np.asarray((2.8e-6, 1.4e-6, 1.8e-6), dtype=np.float64)
        expected_km_row = np.asarray((1.0e-4, 2.0e-4, 1.5e-4), dtype=np.float64)
        expected_vmax = np.asarray((
            (2.8e-6, 2.8e-6, 7.2e-6),
            (2.8e-6, 1.4e-6, 0.0),
        ), dtype=np.float64)
        np.testing.assert_array_equal(
            np.asarray(SPECIES_BASE_UPTAKE_VMAX_MOL_GDW_S),
            base_vmax,
        )
        np.testing.assert_array_equal(
            np.asarray(SPECIES_UPTAKE_KM_MOLAR),
            expected_km_row,
        )
        np.testing.assert_array_equal(config['uptake_vmax'], expected_vmax)
        np.testing.assert_array_equal(
            config['uptake_km'],
            np.tile(expected_km_row, (2, 1)),
        )
        self.assertEqual(config['uptake_vmax'][1, 2], 0.0)
        for name in (
            'uptake_vmax',
            'uptake_km',
            'inhibition_alpha',
            'inhibition_beta',
            'max_inhibition',
        ):
            self.assertEqual(config[name].shape, (2, 3))
            self.assertEqual(config[name].dtype, np.dtype(np.float64))
            self.assertTrue(config[name].flags.c_contiguous)
            self.assertTrue(np.all(np.isfinite(config[name])))

    def test_python_backend_accepts_row_injection_and_preserves_non_utilizer(self):
        result = self.advance('python')
        self.assertEqual(result[0].backend, 'python')
        self.assert_valid_advance(*result)

    @unittest.skipUnless(
        NATIVE_METABOLIC_AVAILABLE,
        'native metabolic kernel is not built',
    )
    def test_native_backend_matches_python_for_translated_batch(self):
        python_result = self.advance('python')
        native_result = self.advance('cpp')
        self.assertEqual(native_result[0].backend_used, 'cpp')
        self.assert_valid_advance(*native_result)

        python_batch = python_result[1]
        native_batch = native_result[1]
        for name in (
            'concentrations_m',
            'biomass_gdw',
            'growth_rate_s',
            'biomass_delta_gdw',
            'substrate_consumed_mol',
            'acetate_secreted_mol',
        ):
            np.testing.assert_allclose(
                native_batch[name],
                python_batch[name],
                rtol=1.0e-12,
                atol=1.0e-15,
            )
        np.testing.assert_allclose(
            np.asarray(native_result[3], dtype=np.float64),
            np.asarray(python_result[3], dtype=np.float64),
            rtol=1.0e-12,
            atol=1.0e-15,
        )

    def test_broadcast_view_is_rejected_but_allocated_row_storage_succeeds(self):
        orchestrator = self.orchestrator('python')
        increment = self.engine.get_molar_substrate_bounds('mdcf_2', 50.0, 1.0, 1.0)
        broadcast_view = np.broadcast_to(increment, (2, 3))
        self.assertEqual(broadcast_view.strides[0], 0)
        self.assertFalse(broadcast_view.flags.writeable)
        self.assertFalse(broadcast_view.flags.owndata)

        rejected_batch = self.translated_batch()
        rejected_batch['concentrations_m'] = broadcast_view
        with self.assertRaisesRegex(ValueError, 'C-contiguous|writable'):
            orchestrator.advance_batch(
                **rejected_batch,
                dt_seconds=0.1,
                max_substeps=8,
            )

        direct_vector_batch = self.translated_batch()
        direct_vector_batch['concentrations_m'] = increment
        with self.assertRaisesRegex(ValueError, 'two-dimensional|shape'):
            orchestrator.advance_batch(
                **direct_vector_batch,
                dt_seconds=0.1,
                max_substeps=8,
            )

        allocated_batch = self.translated_batch()
        summary = orchestrator.advance_batch(
            **allocated_batch,
            dt_seconds=0.1,
            max_substeps=8,
        )
        self.assertIsNotNone(summary)
        self.assertTrue(np.all(np.isfinite(allocated_batch['concentrations_m'])))

    def test_repeated_translation_allocations_are_bounded_and_unaliased(self):
        # Retaining 256 three-double outputs detects aliasing while a five-second
        # budget tolerates instrumented CI; see DOCUMENTATION.md section 10.4.
        call_count = 256
        start = time.perf_counter()
        outputs = [
            self.engine.get_molar_substrate_bounds('mdcf_2', 50.0, 1.0, 1.0)
            for _ in range(call_count)
        ]
        elapsed = time.perf_counter() - start

        self.assertTrue(math.isfinite(elapsed))
        self.assertLess(elapsed, 5.0)
        self.assertEqual(
            sum(output.nbytes for output in outputs),
            call_count * 3 * np.dtype(np.float64).itemsize,
        )
        self.assertEqual(len({output.ctypes.data for output in outputs}), call_count)
        for previous, current in zip(outputs, outputs[1:]):
            self.assertFalse(np.shares_memory(previous, current))
        print(
            f'Diet translation allocation smoke: {call_count} calls in {elapsed:.6f}s; '
            f'{sum(output.nbytes for output in outputs)} payload bytes'
        )

    def test_repeated_python_metabolic_batches_complete_within_smoke_budget(self):
        # Sixty-four cached two-agent batches expose accidental scaling regressions;
        # a ten-second ceiling remains non-microbenchmark-grade; see DOCUMENTATION.md section 10.4.
        iteration_count = 64
        orchestrator = self.orchestrator('python')
        batch = self.translated_batch()
        initial_concentrations = batch['concentrations_m'].copy()
        initial_biomass = batch['biomass_gdw'].copy()

        orchestrator.advance_batch(
            **batch,
            dt_seconds=0.1,
            max_substeps=8,
            cache_enabled=True,
        )
        start = time.perf_counter()
        for _ in range(iteration_count):
            batch['concentrations_m'][...] = initial_concentrations
            batch['biomass_gdw'][...] = initial_biomass
            summary = orchestrator.advance_batch(
                **batch,
                dt_seconds=0.1,
                max_substeps=8,
                cache_enabled=True,
            )
        elapsed = time.perf_counter() - start

        self.assertTrue(math.isfinite(elapsed))
        self.assertLess(elapsed, 10.0)
        self.assertTrue(all(math.isfinite(float(value)) for value in summary))
        self.assertTrue(np.all(np.isfinite(batch['concentrations_m'])))
        self.assertTrue(np.all(batch['concentrations_m'] >= 0.0))
        print(
            f'Python translated metabolic smoke: {iteration_count} batches in '
            f'{elapsed:.6f}s'
        )


if __name__ == '__main__':
    unittest.main()