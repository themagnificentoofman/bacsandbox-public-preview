import importlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from bacteria.metabolism import (
    CACHE_TOLERANCE_MOLAR,
    DEFAULT_MODEL_PATH,
    DEFAULT_PRIOR_PATH,
    MetabolicAssetError,
    MetabolicOrchestrator,
    MetabolicSolverUnavailable,
    PythonMetabolicKernel,
    SPECIES_BASE_UPTAKE_VMAX_MOL_GDW_S,
    SPECIES_UPTAKE_KM_MOLAR,
    SUBSTRATE_ORDER,
    build_neutral_phenotype_config,
    build_species_phenotype_config,
    concentration_dependent_uptake_bound,
    enumerate_neighbor_cache_keys,
    load_metabolic_model,
    load_metabolic_priors,
    posterior_mean_inhibition_multiplier,
    quantize_cache_key,
)


try:
    bacteria_c_kernel = importlib.import_module('bacteria_c_kernel')
except ImportError:
    bacteria_c_kernel = None


NATIVE_METABOLIC_AVAILABLE = (
    bacteria_c_kernel is not None and hasattr(bacteria_c_kernel, 'MetabolicKernel')
)


def read_asset(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def duplicated_model_arguments(cache_capacity=8):
    model = load_metabolic_model()
    reaction_count = model.lower_bounds.size
    metabolite_count = model.model_metabolite_offsets[-1]
    nonzero_count = model.stoich_values.size
    return {
        'model_reaction_offsets': np.asarray((0, reaction_count, 2 * reaction_count), dtype=np.int64),
        'model_metabolite_offsets': np.asarray((0, metabolite_count, 2 * metabolite_count), dtype=np.int64),
        'stoich_row_offsets': np.concatenate((
            model.stoich_row_offsets,
            model.stoich_row_offsets[1:] + nonzero_count,
        )).astype(np.int64),
        'stoich_column_indices': np.concatenate((
            model.stoich_column_indices,
            model.stoich_column_indices + reaction_count,
        )).astype(np.int32),
        'stoich_values': np.tile(model.stoich_values, 2).astype(np.float64),
        'lower_bounds': np.tile(model.lower_bounds, 2).astype(np.float64),
        'upper_bounds': np.tile(model.upper_bounds, 2).astype(np.float64),
        'objective_coefficients': np.tile(model.objective_coefficients, 2).astype(np.float64),
        'exchange_reaction_indices': np.asarray(
            ((0, 1, 2), (8, 9, 10)), dtype=np.int32
        ),
        'biomass_reaction_indices': np.asarray((6, 14), dtype=np.int32),
        'acetate_reaction_indices': np.asarray((7, 15), dtype=np.int32),
        'cache_tolerance_molar': model.cache_tolerance_molar,
        'cache_capacity': cache_capacity,
    }


class TemporaryAssetMixin:
    def load_modified_asset(self, source_path, mutation, loader):
        payload = read_asset(source_path)
        mutation(payload)
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / 'asset.json'
            path.write_text(json.dumps(payload), encoding='utf-8')
            return loader(path)


class TestMetabolicModelLoading(TemporaryAssetMixin, unittest.TestCase):
    def test_canonical_model_metadata_and_topology(self):
        model = load_metabolic_model()

        self.assertEqual(model.schema_version, 1)
        self.assertEqual(model.model_id, 'reduced_metabolism_v1')
        self.assertEqual(model.model_revision, 1)
        self.assertEqual(model.provenance, 'Own (BacSandbox model)')
        self.assertEqual(model.substrate_order, SUBSTRATE_ORDER)
        self.assertIn('not a genome-scale model', model.scope_statement.lower())
        self.assertIn('never literal polymer molarity', model.scope_statement.lower())
        self.assertIn('not empirical', model.design_rationale['acetate_coupling'])
        self.assertIn('not a biological threshold', model.design_rationale['cache_tolerance'])
        self.assertEqual(model.units['substrate_concentration'], 'M glucose-equivalent monomer')
        self.assertIn('mol glucose-equivalent monomer / gDW / s', model.units['reaction_flux'])
        self.assertIn('specific growth rate in 1/s', model.units['reaction_flux'])
        self.assertIn('not a claim of dimensional GEM rigor', model.units['reaction_flux'])
        self.assertTrue(model.requires_explicit_environment_volume_liters)
        self.assertEqual(model.cache_tolerance_molar, CACHE_TOLERANCE_MOLAR)
        self.assertEqual(model.analytical_full_uptake_biomass_flux, 3.0)
        self.assertEqual(model.analytical_full_uptake_acetate_flux, 0.75)
        self.assertEqual(model.metabolite_order[-2:], ('biomass_precursor', 'acetate_internal'))
        self.assertEqual(model.reaction_order[-2:], ('biomass', 'acetate_export'))

    def test_flattened_abi_arrays_have_exact_shapes_dtypes_and_values(self):
        model = load_metabolic_model()
        expected = {
            'model_reaction_offsets': (np.int64, (2,), (0, 8)),
            'model_metabolite_offsets': (np.int64, (2,), (0, 5)),
            'stoich_row_offsets': (np.int64, (6,), (0, 2, 4, 6, 10, 14)),
            'stoich_column_indices': (
                np.int32,
                (14,),
                (0, 3, 1, 4, 2, 5, 3, 4, 5, 6, 3, 4, 5, 7),
            ),
            'stoich_values': (
                np.float64,
                (14,),
                (1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, -1.0,
                 0.25, 0.25, 0.25, -1.0),
            ),
            'lower_bounds': (np.float64, (8,), (0.0,) * 8),
            'upper_bounds': (np.float64, (8,), (1.0,) * 6 + (3.0, 3.0)),
            'objective_coefficients': (
                np.float64,
                (8,),
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
            ),
            'exchange_reaction_indices': (np.int32, (1, 3), ((0, 1, 2),)),
            'biomass_reaction_indices': (np.int32, (1,), (6,)),
            'acetate_reaction_indices': (np.int32, (1,), (7,)),
        }

        for name, (dtype, shape, values) in expected.items():
            with self.subTest(name=name):
                array = getattr(model, name)
                self.assertEqual(array.dtype, np.dtype(dtype))
                self.assertEqual(array.shape, shape)
                self.assertTrue(array.flags.c_contiguous)
                self.assertFalse(array.flags.writeable)
                np.testing.assert_array_equal(array, np.asarray(values, dtype=dtype))

    def test_wrong_substrate_order_rejected(self):
        with self.assertRaisesRegex(MetabolicAssetError, 'substrate_order.*exactly'):
            self.load_modified_asset(
                DEFAULT_MODEL_PATH,
                lambda value: value['substrate_order'].reverse(),
                load_metabolic_model,
            )

    def test_duplicate_reaction_name_rejected(self):
        def duplicate_name(value):
            value['reaction_order'][1] = value['reaction_order'][0]

        with self.assertRaisesRegex(MetabolicAssetError, 'reaction_order.*unique'):
            self.load_modified_asset(DEFAULT_MODEL_PATH, duplicate_name, load_metabolic_model)

    def test_stoichiometry_index_and_duplicate_coordinates_rejected(self):
        def invalid_index(value):
            value['stoichiometry'][0]['reaction_index'] = len(value['reaction_order'])

        with self.assertRaisesRegex(MetabolicAssetError, 'reaction_index.*\[0, 8\)'):
            self.load_modified_asset(DEFAULT_MODEL_PATH, invalid_index, load_metabolic_model)

        def duplicate_coordinate(value):
            value['stoichiometry'].append(dict(value['stoichiometry'][0]))

        with self.assertRaisesRegex(MetabolicAssetError, 'duplicate stoichiometry'):
            self.load_modified_asset(DEFAULT_MODEL_PATH, duplicate_coordinate, load_metabolic_model)

    def test_bound_order_and_marker_indices_rejected(self):
        def invalid_bounds(value):
            value['lower_bounds'][0] = value['upper_bounds'][0] + 1.0

        with self.assertRaisesRegex(MetabolicAssetError, 'bound ordering'):
            self.load_modified_asset(DEFAULT_MODEL_PATH, invalid_bounds, load_metabolic_model)

        def invalid_marker(value):
            value['biomass_reaction_index'] = value['exchange_reaction_indices'][0]

        with self.assertRaisesRegex(MetabolicAssetError, 'must be distinct'):
            self.load_modified_asset(DEFAULT_MODEL_PATH, invalid_marker, load_metabolic_model)

    def test_unknown_keys_wrong_revision_nonfinite_and_units_rejected(self):
        cases = (
            ('unexpected', lambda value: value.__setitem__('unexpected', True), 'keys differ'),
            ('revision', lambda value: value.__setitem__('model_revision', 2), 'unsupported model'),
            ('unit', lambda value: value['units'].__setitem__('substrate_concentration', 'M'), 'glucose-equivalent'),
            ('nonfinite', lambda value: value['lower_bounds'].__setitem__(0, float('nan')), 'finite'),
        )
        for name, mutation, message in cases:
            with self.subTest(name=name):
                with self.assertRaisesRegex(MetabolicAssetError, message):
                    self.load_modified_asset(
                        DEFAULT_MODEL_PATH,
                        mutation,
                        load_metabolic_model,
                    )


class TestMetabolicPriorLoading(TemporaryAssetMixin, unittest.TestCase):
    def test_production_priors_are_ordered_neutral_and_nonempirical(self):
        priors = load_metabolic_priors()

        self.assertEqual(priors.prior_id, 'metabolic_priors_v1')
        self.assertEqual(priors.prior_revision, 1)
        self.assertEqual(priors.substrate_order, SUBSTRATE_ORDER)
        self.assertIn('no empirical interaction', priors.scope.lower())
        self.assertIn('no empirical interaction', priors.parameter_rationale['alpha_beta'])
        self.assertIn('exactly neutral', priors.parameter_rationale['max_inhibition'])
        self.assertEqual(tuple(row.substrate for row in priors.rows), SUBSTRATE_ORDER)
        for row in priors.rows:
            self.assertGreater(row.alpha, 0.0)
            self.assertGreater(row.beta, 0.0)
            self.assertEqual(row.max_inhibition, 0.0)
            self.assertIn('no empirical interaction', row.evidence_scope)

    def test_nonpositive_parameters_and_nonneutral_production_prior_rejected(self):
        cases = (
            ('alpha', -1.0, 'positive'),
            ('beta', 0.0, 'positive'),
            # Nonzero inhibition is synthetic in this test and is forbidden in production assets.
            ('max_inhibition', 0.25, 'must be zero'),
        )
        for field, replacement, message in cases:
            with self.subTest(field=field):
                def mutate(value, field=field, replacement=replacement):
                    value['rows'][0][field] = replacement

                with self.assertRaisesRegex(MetabolicAssetError, message):
                    self.load_modified_asset(
                        DEFAULT_PRIOR_PATH,
                        mutate,
                        load_metabolic_priors,
                    )

    def test_prior_order_duplicate_rows_and_revision_rejected(self):
        def wrong_order(value):
            value['rows'][0], value['rows'][1] = value['rows'][1], value['rows'][0]

        with self.assertRaisesRegex(MetabolicAssetError, 'exact substrate order'):
            self.load_modified_asset(DEFAULT_PRIOR_PATH, wrong_order, load_metabolic_priors)

        def duplicate_row(value):
            value['rows'].append(dict(value['rows'][0]))

        with self.assertRaisesRegex(MetabolicAssetError, 'duplicate prior row'):
            self.load_modified_asset(DEFAULT_PRIOR_PATH, duplicate_row, load_metabolic_priors)

        with self.assertRaisesRegex(MetabolicAssetError, 'unsupported prior'):
            self.load_modified_asset(
                DEFAULT_PRIOR_PATH,
                lambda value: value.__setitem__('prior_revision', 2),
                load_metabolic_priors,
            )


class TestMetabolicArchitectureHelpers(unittest.TestCase):
    def test_cache_quantization_uses_floor_at_boundaries(self):
        tolerance = CACHE_TOLERANCE_MOLAR
        below_two_bins = np.nextafter(2.0 * tolerance, 0.0)

        self.assertEqual(quantize_cache_key((0.0, tolerance, below_two_bins)), (0, 1, 1))
        self.assertEqual(quantize_cache_key((tolerance, 2.0 * tolerance, 3.0 * tolerance)), (1, 2, 3))

    def test_neighbor_enumeration_is_stable_and_complete(self):
        neighbors = enumerate_neighbor_cache_keys((0, 4, 9))

        self.assertEqual(len(neighbors), 27)
        self.assertEqual(len(set(neighbors)), 27)
        self.assertEqual(neighbors[0], (-1, 3, 8))
        self.assertEqual(neighbors[13], (0, 4, 9))
        self.assertEqual(neighbors[-1], (1, 5, 10))

    def test_cache_helper_rejects_invalid_shapes_ranges_and_numbers(self):
        invalid_calls = (
            lambda: quantize_cache_key((1.0, 2.0)),
            lambda: quantize_cache_key((-1.0, 0.0, 0.0)),
            lambda: quantize_cache_key((0.0, 0.0, float('nan'))),
            lambda: quantize_cache_key((0.0, 0.0, 0.0), 0.0),
            lambda: quantize_cache_key((float(np.finfo(np.float64).max), 0.0, 0.0)),
            lambda: enumerate_neighbor_cache_keys((0, 1)),
            lambda: enumerate_neighbor_cache_keys((0, -1, 2)),
            lambda: enumerate_neighbor_cache_keys((0, 1.5, 2)),
            lambda: enumerate_neighbor_cache_keys((0, 1, np.iinfo(np.int64).max)),
        )
        for call in invalid_calls:
            with self.subTest(call=call):
                with self.assertRaises((TypeError, ValueError)):
                    call()


class TestPythonMetabolicSolver(unittest.TestCase):
    def make_kernel(self, cache_capacity=8, phenotype_count=1):
        model = load_metabolic_model()
        kernel = PythonMetabolicKernel(
            model.model_reaction_offsets,
            model.model_metabolite_offsets,
            model.stoich_row_offsets,
            model.stoich_column_indices,
            model.stoich_values,
            model.lower_bounds,
            model.upper_bounds,
            model.objective_coefficients,
            model.exchange_reaction_indices,
            model.biomass_reaction_indices,
            model.acetate_reaction_indices,
            model.cache_tolerance_molar,
            cache_capacity,
        )
        kernel.configure_phenotypes(
            np.zeros(phenotype_count, dtype=np.int32),
            np.full((phenotype_count, 3), 2.0, dtype=np.float64),
            np.ones((phenotype_count, 3), dtype=np.float64),
            np.ones((phenotype_count, 3), dtype=np.float64),
            np.ones((phenotype_count, 3), dtype=np.float64),
            np.zeros((phenotype_count, 3), dtype=np.float64),
            model.model_revision,
            1,
        )
        return model, kernel

    def make_batch(self, agent_count=2, environment_count=1, concentration=2.0e-4):
        return {
            'phenotype_id': np.zeros(agent_count, dtype=np.int32),
            'environment_id': np.arange(agent_count, dtype=np.int32) % environment_count,
            'concentrations_m': np.full(
                (environment_count, 3), concentration, dtype=np.float64
            ),
            'environment_volume_l': np.full(environment_count, 1.0e-6, dtype=np.float64),
            'biomass_gdw': np.full(agent_count, 5.0e-7, dtype=np.float64),
            'growth_rate_s': np.full(agent_count, -1.0, dtype=np.float64),
            'biomass_delta_gdw': np.full(agent_count, -2.0, dtype=np.float64),
            'substrate_consumed_mol': np.full((agent_count, 3), -3.0, dtype=np.float64),
            'acetate_secreted_mol': np.full(agent_count, -4.0, dtype=np.float64),
        }

    def test_reduced_lp_full_partial_zero_and_mass_balance(self):
        model, kernel = self.make_kernel()
        cases = (
            (np.ones(3, dtype=np.float64), np.asarray((1.0, 1.0, 1.0))),
            (np.asarray((1.0, 1.0 / 3.0, 0.0)), np.asarray((1.0, 0.5, 0.0))),
            (np.zeros(3, dtype=np.float64), np.zeros(3, dtype=np.float64)),
        )
        dense = np.zeros((5, 8), dtype=np.float64)
        for row in range(5):
            start = model.stoich_row_offsets[row]
            stop = model.stoich_row_offsets[row + 1]
            dense[row, model.stoich_column_indices[start:stop]] = model.stoich_values[start:stop]

        for concentrations, uptake in cases:
            with self.subTest(concentrations=concentrations):
                flux, cache_hit = kernel.solve_flux(0, concentrations, cache_enabled=False)
                expected = np.concatenate((
                    uptake,
                    uptake,
                    (np.sum(uptake), 0.25 * np.sum(uptake)),
                ))
                self.assertFalse(cache_hit)
                np.testing.assert_allclose(flux, expected, rtol=0.0, atol=1.0e-10)
                self.assertLessEqual(np.max(np.abs(dense @ flux)), 1.0e-10)

    def test_synthetic_posterior_mean_constrains_exchange(self):
        model, kernel = self.make_kernel()
        kernel.configure_phenotypes(
            np.asarray((0,), dtype=np.int32),
            np.ones((1, 3), dtype=np.float64),
            np.full((1, 3), 3.0, dtype=np.float64),
            np.full((1, 3), 2.0, dtype=np.float64),
            np.full((1, 3), 3.0, dtype=np.float64),
            np.full((1, 3), 0.5, dtype=np.float64),
            model.model_revision,
            7,
        )

        flux, _ = kernel.solve_flux(0, np.full(3, 2.0, dtype=np.float64), False)

        np.testing.assert_allclose(flux[:3], 0.32, rtol=0.0, atol=1.0e-12)
        self.assertAlmostEqual(flux[6], 0.96)
        self.assertAlmostEqual(flux[7], 0.24)

    def test_concatenated_models_and_explicit_infeasible_status(self):
        arguments = duplicated_model_arguments()
        kernel = PythonMetabolicKernel(**arguments)
        kernel.configure_phenotypes(
            np.asarray((0, 1), dtype=np.int32),
            np.full((2, 3), 2.0, dtype=np.float64),
            np.ones((2, 3), dtype=np.float64),
            np.ones((2, 3), dtype=np.float64),
            np.ones((2, 3), dtype=np.float64),
            np.zeros((2, 3), dtype=np.float64),
            1,
            1,
        )
        first, _ = kernel.solve_flux(0, np.ones(3, dtype=np.float64), False)
        second, _ = kernel.solve_flux(1, np.ones(3, dtype=np.float64), False)
        np.testing.assert_array_equal(first[8:], 0.0)
        np.testing.assert_array_equal(second[:8], 0.0)
        np.testing.assert_allclose(first[:8], second[8:], rtol=0.0, atol=1.0e-10)

        model = load_metabolic_model()
        lower = np.array(model.lower_bounds, copy=True)
        lower[0] = 0.5
        infeasible = PythonMetabolicKernel(
            model.model_reaction_offsets,
            model.model_metabolite_offsets,
            model.stoich_row_offsets,
            model.stoich_column_indices,
            model.stoich_values,
            lower,
            model.upper_bounds,
            model.objective_coefficients,
            model.exchange_reaction_indices,
            model.biomass_reaction_indices,
            model.acetate_reaction_indices,
            model.cache_tolerance_molar,
            8,
        )
        infeasible.configure_phenotypes(
            np.asarray((0,), dtype=np.int32),
            np.ones((1, 3), dtype=np.float64),
            np.ones((1, 3), dtype=np.float64),
            np.ones((1, 3), dtype=np.float64),
            np.ones((1, 3), dtype=np.float64),
            np.zeros((1, 3), dtype=np.float64),
            1,
            1,
        )
        with self.assertRaisesRegex(RuntimeError, 'infeasible'):
            infeasible.solve_flux(0, np.zeros(3, dtype=np.float64), False)

        infinite_upper = np.array(model.upper_bounds, copy=True)
        infinite_upper[0] = np.inf
        with self.assertRaisesRegex(ValueError, 'finite bounds'):
            PythonMetabolicKernel(
                model.model_reaction_offsets,
                model.model_metabolite_offsets,
                model.stoich_row_offsets,
                model.stoich_column_indices,
                model.stoich_values,
                model.lower_bounds,
                infinite_upper,
                model.objective_coefficients,
                model.exchange_reaction_indices,
                model.biomass_reaction_indices,
                model.acetate_reaction_indices,
                model.cache_tolerance_molar,
                8,
            )

    def test_cache_tolerance_phenotype_lru_and_disabled_counters(self):
        tolerance = CACHE_TOLERANCE_MOLAR
        _, kernel = self.make_kernel(cache_capacity=2, phenotype_count=2)
        capped_a = np.full(3, 100000.75 * tolerance, dtype=np.float64)
        capped_b = np.full(3, 100001.25 * tolerance, dtype=np.float64)

        _, first_hit = kernel.solve_flux(0, capped_a)
        _, adjacent_hit = kernel.solve_flux(0, capped_b)
        _, phenotype_hit = kernel.solve_flux(1, capped_b)

        self.assertFalse(first_hit)
        self.assertTrue(adjacent_hit)
        self.assertFalse(phenotype_hit)
        self.assertEqual(kernel.cache_stats(), (2, 1, 2, 0, 2))

        kernel.clear_cache()
        unsaturated_a = np.full(3, 2.25 * tolerance, dtype=np.float64)
        unsaturated_b = np.full(3, 2.75 * tolerance, dtype=np.float64)
        first_flux, first_hit = kernel.solve_flux(0, unsaturated_a)
        cached_flux, within_tolerance_hit = kernel.solve_flux(0, unsaturated_b)
        self.assertFalse(first_hit)
        self.assertTrue(within_tolerance_hit)
        np.testing.assert_array_equal(cached_flux, first_flux)
        self.assertEqual(kernel.cache_stats(), (1, 1, 1, 0, 1))

        kernel.clear_cache()
        concentrations = tuple(np.full(3, value, dtype=np.float64) for value in (0.1, 0.2, 0.3))
        kernel.solve_flux(0, concentrations[0])
        kernel.solve_flux(0, concentrations[1])
        kernel.solve_flux(0, concentrations[0])
        kernel.solve_flux(0, concentrations[2])
        self.assertEqual(kernel.cache_stats(), (2, 1, 3, 1, 3))
        _, evicted_hit = kernel.solve_flux(0, concentrations[1])
        self.assertFalse(evicted_hit)
        self.assertEqual(kernel.cache_stats(), (2, 1, 4, 2, 4))

        kernel.clear_cache()
        kernel.solve_flux(0, concentrations[0], cache_enabled=False)
        kernel.solve_flux(0, concentrations[0], cache_enabled=False)
        self.assertEqual(kernel.cache_stats(), (0, 0, 0, 0, 2))

        model = load_metabolic_model()
        kernel.configure_phenotypes(
            np.asarray((0,), dtype=np.int32),
            np.full((1, 3), 2.0, dtype=np.float64),
            np.ones((1, 3), dtype=np.float64),
            np.ones((1, 3), dtype=np.float64),
            np.ones((1, 3), dtype=np.float64),
            np.zeros((1, 3), dtype=np.float64),
            model.model_revision + 1,
            9,
        )
        self.assertEqual(kernel.cache_stats(), (0, 0, 0, 0, 0))
        _, reconfigured_hit = kernel.solve_flux(0, concentrations[0])
        self.assertFalse(reconfigured_hit)
        self.assertEqual(kernel.cache_stats(), (1, 0, 1, 0, 1))

    def test_adaptive_shared_environment_mass_accounting_and_weighted_growth(self):
        _, kernel = self.make_kernel()
        batch = self.make_batch(agent_count=3, environment_count=2, concentration=8.0e-5)
        batch['environment_id'][:] = (0, 0, 1)
        initial_concentrations = batch['concentrations_m'].copy()
        initial_biomass = batch['biomass_gdw'].copy()

        summary = kernel.advance_batch(
            **batch,
            dt_seconds=0.75,
            depletion_safety=0.8,
            max_substeps=32,
            cache_enabled=True,
        )

        self.assertGreater(summary[0], 1)
        self.assertGreater(summary[3], summary[1])
        self.assertLess(summary[5], 0.75)
        self.assertTrue(np.all(batch['concentrations_m'] >= 0.0))
        for environment in range(2):
            members = batch['environment_id'] == environment
            depleted = (
                initial_concentrations[environment] - batch['concentrations_m'][environment]
            ) * batch['environment_volume_l'][environment]
            np.testing.assert_allclose(
                np.sum(batch['substrate_consumed_mol'][members], axis=0),
                depleted,
                rtol=0.0,
                atol=1.0e-18,
            )
        np.testing.assert_allclose(
            batch['biomass_gdw'] - initial_biomass,
            batch['biomass_delta_gdw'],
            rtol=0.0,
            atol=1.0e-18,
        )
        self.assertTrue(np.all(batch['growth_rate_s'] >= 0.0))
        self.assertTrue(np.all(batch['acetate_secreted_mol'] >= 0.0))

    def test_low_substrate_threshold_tightens_only_below_boundary(self):
        _, above_kernel = self.make_kernel()
        _, below_kernel = self.make_kernel()
        above = self.make_batch(agent_count=1, environment_count=1, concentration=1.1e-4)
        below = self.make_batch(agent_count=1, environment_count=1, concentration=0.9e-4)

        above_summary = above_kernel.advance_batch(
            **above,
            dt_seconds=0.9,
            depletion_safety=0.8,
            max_substeps=8,
            cache_enabled=False,
        )
        below_summary = below_kernel.advance_batch(
            **below,
            dt_seconds=0.9,
            depletion_safety=0.8,
            max_substeps=8,
            cache_enabled=False,
        )

        self.assertEqual(above_summary[0], 1)
        self.assertGreater(below_summary[0], 1)
        self.assertEqual(above_summary[5], 0.9)
        self.assertLess(below_summary[5], 0.9)

    def test_dt_zero_and_max_substep_failure_are_transactional(self):
        _, kernel = self.make_kernel()
        zero_batch = self.make_batch(concentration=5.0e-5)
        original_concentrations = zero_batch['concentrations_m'].copy()
        original_biomass = zero_batch['biomass_gdw'].copy()

        summary = kernel.advance_batch(
            **zero_batch,
            dt_seconds=0.0,
            depletion_safety=0.8,
            max_substeps=2,
            cache_enabled=True,
        )

        self.assertEqual(summary, (0, 0, 0, 0, 0, 0.0))
        np.testing.assert_array_equal(zero_batch['concentrations_m'], original_concentrations)
        np.testing.assert_array_equal(zero_batch['biomass_gdw'], original_biomass)
        for name in ('growth_rate_s', 'biomass_delta_gdw', 'substrate_consumed_mol', 'acetate_secreted_mol'):
            np.testing.assert_array_equal(zero_batch[name], np.zeros_like(zero_batch[name]))
        self.assertEqual(kernel.cache_stats(), (0, 0, 0, 0, 0))

        failing_batch = self.make_batch(concentration=5.0e-5)
        snapshots = {name: value.copy() for name, value in failing_batch.items() if value.flags.writeable}
        with self.assertRaisesRegex(RuntimeError, 'max_substeps'):
            kernel.advance_batch(
                **failing_batch,
                dt_seconds=1.0,
                depletion_safety=0.5,
                max_substeps=1,
                cache_enabled=False,
            )
        for name, snapshot in snapshots.items():
            np.testing.assert_array_equal(failing_batch[name], snapshot)

    def test_posterior_multiplier_matches_synthetic_analytical_case(self):
        # Synthetic alpha=2, beta=3, and cap=0.5 give 1-(2/5)*0.5=0.8; see DOCUMENTATION.md §10.3.
        self.assertAlmostEqual(posterior_mean_inhibition_multiplier(2.0, 3.0, 0.5), 0.8)
        self.assertEqual(posterior_mean_inhibition_multiplier(1.0, 1.0, 0.0), 1.0)

    def test_uptake_bound_matches_synthetic_analytical_case(self):
        # Synthetic vmax=10, C=2, km=3, and multiplier=0.8 give 10*2/(3+2)*0.8=3.2; see DOCUMENTATION.md §10.3.
        self.assertAlmostEqual(
            concentration_dependent_uptake_bound(10.0, 2.0, 3.0, 0.8),
            3.2,
        )
        self.assertEqual(concentration_dependent_uptake_bound(10.0, 0.0, 3.0), 0.0)

    def test_formula_helpers_reject_invalid_ranges_and_nonfinite_values(self):
        invalid_calls = (
            lambda: posterior_mean_inhibition_multiplier(0.0, 1.0, 0.0),
            lambda: posterior_mean_inhibition_multiplier(1.0, 1.0, 1.1),
            lambda: posterior_mean_inhibition_multiplier(float('inf'), 1.0, 0.0),
            lambda: posterior_mean_inhibition_multiplier(1.0e308, 1.0e308, 0.5),
            lambda: concentration_dependent_uptake_bound(-1.0, 1.0, 1.0),
            lambda: concentration_dependent_uptake_bound(1.0, -1.0, 1.0),
            lambda: concentration_dependent_uptake_bound(1.0, 1.0, 0.0),
            lambda: concentration_dependent_uptake_bound(1.0, 1.0, 1.0, -0.1),
        )
        for call in invalid_calls:
            with self.subTest(call=call):
                with self.assertRaises((TypeError, ValueError)):
                    call()


class TestMetabolicOrchestrator(unittest.TestCase):
    def configured_arguments(self):
        model = load_metabolic_model()
        priors = load_metabolic_priors()
        return {
            'model': model,
            'priors': priors,
            'phenotype_config': build_neutral_phenotype_config(
                np.full((1, 3), 2.0, dtype=np.float64),
                np.ones((1, 3), dtype=np.float64),
                model,
                priors,
                8,
            ),
        }

    def make_batch(self):
        return TestPythonMetabolicSolver().make_batch(
            agent_count=1,
            environment_count=1,
            concentration=2.0e-4,
        )

    def configured_species_arguments(self):
        model = load_metabolic_model()
        priors = load_metabolic_priors()
        # Own-model multipliers represent one high-capacity and one non-utilizing MDCF
        # profile without asserting measured kinetics; see DOCUMENTATION.md section 10.4.
        profiles = (
            {'aff_calorie': '1', 'aff_glycan': '2', 'aff_mdcf': '4'},
            {'aff_calorie': '1', 'aff_glycan': '1', 'aff_mdcf': '0'},
        )
        return {
            'model': model,
            'priors': priors,
            'phenotype_config': build_species_phenotype_config(
                profiles,
                model,
                priors,
                8,
            ),
        }

    def make_species_batch(self):
        batch = TestPythonMetabolicSolver().make_batch(
            agent_count=2,
            environment_count=2,
            concentration=0.0,
        )
        batch['phenotype_id'][:] = (0, 1)
        # Equal MDCF pools isolate phenotype capacity as the only uptake difference;
        # see DOCUMENTATION.md section 10.4.
        batch['concentrations_m'][:, 2] = 2.0e-4
        return batch

    def test_species_affinities_scale_vmax_without_changing_km(self):
        config = self.configured_species_arguments()['phenotype_config']
        base_vmax = np.asarray(SPECIES_BASE_UPTAKE_VMAX_MOL_GDW_S, dtype=np.float64)
        expected_vmax = np.asarray(((1.0, 2.0, 4.0), (1.0, 1.0, 0.0))) * base_vmax
        expected_km = np.tile(np.asarray(SPECIES_UPTAKE_KM_MOLAR), (2, 1))

        np.testing.assert_array_equal(config['uptake_vmax'], expected_vmax)
        np.testing.assert_array_equal(config['uptake_km'], expected_km)
        self.assertEqual(config['uptake_vmax'][1, 2], 0.0)
        for name in ('uptake_vmax', 'uptake_km'):
            self.assertEqual(config[name].dtype, np.float64)
            self.assertTrue(config[name].flags.c_contiguous)

    def test_python_species_batch_preserves_exact_mdcf_non_utilization(self):
        orchestrator = MetabolicOrchestrator(
            mode='fba',
            backend='python',
            **self.configured_species_arguments(),
        )
        batch = self.make_species_batch()
        initial_concentrations = batch['concentrations_m'].copy()
        orchestrator.advance_batch(
            **batch,
            dt_seconds=0.1,
            depletion_safety=0.9,
            max_substeps=8,
            cache_enabled=False,
        )

        self.assertGreater(batch['substrate_consumed_mol'][0, 2], 0.0)
        self.assertEqual(batch['substrate_consumed_mol'][1, 2], 0.0)
        self.assertLess(batch['concentrations_m'][0, 2], initial_concentrations[0, 2])
        self.assertEqual(batch['concentrations_m'][1, 2], initial_concentrations[1, 2])

    def test_mode_and_backend_validation(self):
        with self.assertRaisesRegex(ValueError, 'legacy.*shadow.*fba'):
            MetabolicOrchestrator(mode='invalid')
        with self.assertRaisesRegex(ValueError, 'auto.*python.*cpp'):
            MetabolicOrchestrator(backend='invalid')

    def test_default_legacy_is_true_noop_and_never_imports_native(self):
        with patch('bacteria.metabolism.importlib.import_module') as import_module:
            orchestrator = MetabolicOrchestrator()
            self.assertIsNone(orchestrator.solve_flux())
            self.assertIsNone(orchestrator.advance_batch())
            self.assertIsNone(orchestrator.native_available)

        import_module.assert_not_called()
        self.assertEqual(orchestrator.mode, 'legacy')
        self.assertEqual(orchestrator.backend, 'auto')
        self.assertIsNone(orchestrator.backend_used)
        self.assertIsNone(orchestrator.model_revision)
        self.assertIsNone(orchestrator.prior_revision)

    def test_shadow_and_fba_require_explicit_configuration(self):
        for mode in ('shadow', 'fba'):
            with self.subTest(mode=mode):
                with self.assertRaisesRegex(ValueError, 'model.*priors.*phenotype_config'):
                    MetabolicOrchestrator(mode=mode)

    def test_python_backend_never_imports_and_solve_is_lazy(self):
        orchestrator = MetabolicOrchestrator(
            mode='shadow',
            backend='python',
            **self.configured_arguments(),
        )
        with patch('bacteria.metabolism.importlib.import_module') as importer:
            self.assertIsNone(orchestrator.backend_used)
            flux, cache_hit = orchestrator.solve_flux(
                concentrations_m=np.ones(3, dtype=np.float64)
            )
        importer.assert_not_called()
        self.assertFalse(cache_hit)
        self.assertEqual(flux[6], 3.0)
        self.assertEqual(orchestrator.backend_used, 'python')
        self.assertFalse(orchestrator.native_available)

    def test_shadow_restores_state_and_fba_commits_diagnostics(self):
        for mode in ('shadow', 'fba'):
            with self.subTest(mode=mode):
                orchestrator = MetabolicOrchestrator(
                    mode=mode,
                    backend='python',
                    **self.configured_arguments(),
                )
                batch = self.make_batch()
                initial_concentrations = batch['concentrations_m'].copy()
                initial_biomass = batch['biomass_gdw'].copy()
                summary = orchestrator.advance_batch(
                    **batch,
                    dt_seconds=0.1,
                    depletion_safety=0.9,
                    max_substeps=16,
                    cache_enabled=True,
                )
                self.assertGreaterEqual(summary[0], 1)
                self.assertTrue(np.all(batch['substrate_consumed_mol'] > 0.0))
                if mode == 'shadow':
                    np.testing.assert_array_equal(batch['concentrations_m'], initial_concentrations)
                    np.testing.assert_array_equal(batch['biomass_gdw'], initial_biomass)
                else:
                    self.assertTrue(np.all(batch['concentrations_m'] < initial_concentrations))
                    self.assertTrue(np.all(batch['biomass_gdw'] > initial_biomass))

    def test_auto_falls_back_only_for_missing_native_and_propagates_runtime(self):
        arguments = self.configured_arguments()
        fallback = MetabolicOrchestrator(mode='shadow', backend='auto', **arguments)
        missing_module = ModuleNotFoundError(
            "No module named 'bacteria_c_kernel'",
            name='bacteria_c_kernel',
        )
        with patch(
            'bacteria.metabolism.importlib.import_module',
            side_effect=missing_module,
        ):
            fallback.solve_flux(concentrations_m=np.ones(3, dtype=np.float64))
        self.assertEqual(fallback.backend_used, 'python')
        self.assertFalse(fallback.native_available)

        broken_dependency = ModuleNotFoundError(
            "No module named 'metabolic_dependency'",
            name='metabolic_dependency',
        )
        broken = MetabolicOrchestrator(mode='shadow', backend='auto', **arguments)
        with patch(
            'bacteria.metabolism.importlib.import_module',
            side_effect=broken_dependency,
        ):
            with self.assertRaisesRegex(ModuleNotFoundError, 'metabolic_dependency'):
                broken.solve_flux(concentrations_m=np.ones(3, dtype=np.float64))

        class FailingKernel:
            def __init__(self, **kwargs):
                pass

            def configure_phenotypes(self, *args):
                pass

            def solve_flux(self, *args):
                raise RuntimeError('native solve failed')

        native_module = type('NativeModule', (), {'MetabolicKernel': FailingKernel})()
        propagating = MetabolicOrchestrator(mode='shadow', backend='auto', **arguments)
        with patch('bacteria.metabolism.importlib.import_module', return_value=native_module):
            with self.assertRaisesRegex(RuntimeError, 'native solve failed'):
                propagating.solve_flux(concentrations_m=np.ones(3, dtype=np.float64))
        self.assertEqual(propagating.backend_used, 'cpp')
        self.assertTrue(propagating.native_available)

    def test_native_availability_property_respects_backend(self):
        arguments = self.configured_arguments()
        python_orchestrator = MetabolicOrchestrator(
            mode='shadow',
            backend='python',
            **arguments,
        )
        self.assertFalse(python_orchestrator.native_available)

        auto_orchestrator = MetabolicOrchestrator(mode='shadow', backend='auto', **arguments)
        self.assertIsNone(auto_orchestrator.native_available)


@unittest.skipUnless(NATIVE_METABOLIC_AVAILABLE, 'native MetabolicKernel is not built')
class TestNativeMetabolicKernel(unittest.TestCase):
    def constructor_arguments(self, cache_capacity=8):
        model = load_metabolic_model()
        return {
            'model_reaction_offsets': model.model_reaction_offsets,
            'model_metabolite_offsets': model.model_metabolite_offsets,
            'stoich_row_offsets': model.stoich_row_offsets,
            'stoich_column_indices': model.stoich_column_indices,
            'stoich_values': model.stoich_values,
            'lower_bounds': model.lower_bounds,
            'upper_bounds': model.upper_bounds,
            'objective_coefficients': model.objective_coefficients,
            'exchange_reaction_indices': model.exchange_reaction_indices,
            'biomass_reaction_indices': model.biomass_reaction_indices,
            'acetate_reaction_indices': model.acetate_reaction_indices,
            'cache_tolerance_molar': model.cache_tolerance_molar,
            'cache_capacity': cache_capacity,
        }

    def make_kernel(self, cache_capacity=8, phenotype_count=1):
        kernel = bacteria_c_kernel.MetabolicKernel(
            **self.constructor_arguments(cache_capacity)
        )
        self.configure_kernel(kernel, phenotype_count)
        return kernel

    def configure_kernel(self, kernel, phenotype_count=1):
        kernel.configure_phenotypes(
            np.zeros(phenotype_count, dtype=np.int32),
            np.full((phenotype_count, 3), 2.0, dtype=np.float64),
            np.ones((phenotype_count, 3), dtype=np.float64),
            np.ones((phenotype_count, 3), dtype=np.float64),
            np.ones((phenotype_count, 3), dtype=np.float64),
            np.zeros((phenotype_count, 3), dtype=np.float64),
            1,
            1,
        )

    def make_batch(self, agent_count=3, environment_count=2, concentration=8.0e-5):
        return TestPythonMetabolicSolver().make_batch(
            agent_count=agent_count,
            environment_count=environment_count,
            concentration=concentration,
        )

    def test_constructor_configuration_and_cache_management(self):
        kernel = bacteria_c_kernel.MetabolicKernel(**self.constructor_arguments())

        self.assertEqual(kernel.cache_stats(), (0, 0, 0, 0, 0))
        self.configure_kernel(kernel)
        kernel.clear_cache()
        self.assertEqual(kernel.cache_stats(), (0, 0, 0, 0, 0))

    def test_constructor_and_configuration_validate_arrays(self):
        arguments = self.constructor_arguments()
        arguments['stoich_values'] = arguments['stoich_values'].astype(np.float32)
        with self.assertRaisesRegex(TypeError, 'stoich_values.*dtype'):
            bacteria_c_kernel.MetabolicKernel(**arguments)

        kernel = self.make_kernel()
        with self.assertRaisesRegex(ValueError, 'shape'):
            kernel.configure_phenotypes(
                np.asarray((0,), dtype=np.int32),
                np.ones((1, 2), dtype=np.float64),
                np.ones((1, 3), dtype=np.float64),
                np.ones((1, 3), dtype=np.float64),
                np.ones((1, 3), dtype=np.float64),
                np.zeros((1, 3), dtype=np.float64),
                1,
                1,
            )

    def test_native_solve_matches_python_full_partial_zero_and_cache(self):
        kernel = self.make_kernel()
        _, python_kernel = TestPythonMetabolicSolver().make_kernel()
        cases = (
            np.ones(3, dtype=np.float64),
            np.asarray((1.0, 1.0 / 3.0, 0.0), dtype=np.float64),
            np.zeros(3, dtype=np.float64),
        )
        for concentrations in cases:
            with self.subTest(concentrations=concentrations):
                native_flux, native_hit = kernel.solve_flux(
                    0, concentrations, cache_enabled=False
                )
                python_flux, python_hit = python_kernel.solve_flux(
                    0, concentrations, cache_enabled=False
                )
                self.assertEqual(native_hit, python_hit)
                np.testing.assert_allclose(
                    native_flux, python_flux, rtol=0.0, atol=1.0e-10
                )

        kernel.clear_cache()
        concentration = np.ones(3, dtype=np.float64)
        _, first_hit = kernel.solve_flux(0, concentration)
        _, second_hit = kernel.solve_flux(0, concentration)
        self.assertFalse(first_hit)
        self.assertTrue(second_hit)
        self.assertEqual(kernel.cache_stats(), (1, 1, 1, 0, 1))

    def test_native_concatenated_models_and_infeasible_status(self):
        kernel = bacteria_c_kernel.MetabolicKernel(**duplicated_model_arguments())
        kernel.configure_phenotypes(
            np.asarray((0, 1), dtype=np.int32),
            np.full((2, 3), 2.0, dtype=np.float64),
            np.ones((2, 3), dtype=np.float64),
            np.ones((2, 3), dtype=np.float64),
            np.ones((2, 3), dtype=np.float64),
            np.zeros((2, 3), dtype=np.float64),
            1,
            1,
        )
        first, _ = kernel.solve_flux(0, np.ones(3, dtype=np.float64), False)
        second, _ = kernel.solve_flux(1, np.ones(3, dtype=np.float64), False)
        np.testing.assert_array_equal(first[8:], 0.0)
        np.testing.assert_array_equal(second[:8], 0.0)
        np.testing.assert_allclose(first[:8], second[8:], rtol=0.0, atol=1.0e-10)

        arguments = self.constructor_arguments()
        arguments['lower_bounds'] = np.array(arguments['lower_bounds'], copy=True)
        arguments['lower_bounds'][0] = 0.5
        infeasible = bacteria_c_kernel.MetabolicKernel(**arguments)
        self.configure_kernel(infeasible)
        with self.assertRaisesRegex(RuntimeError, 'infeasible'):
            infeasible.solve_flux(0, np.zeros(3, dtype=np.float64), False)

    def test_native_cache_adjacent_bins_lru_phenotype_and_disabled(self):
        tolerance = CACHE_TOLERANCE_MOLAR
        kernel = self.make_kernel(cache_capacity=2, phenotype_count=2)
        capped_a = np.full(3, 100000.75 * tolerance, dtype=np.float64)
        capped_b = np.full(3, 100001.25 * tolerance, dtype=np.float64)

        kernel.solve_flux(0, capped_a)
        _, adjacent_hit = kernel.solve_flux(0, capped_b)
        _, phenotype_hit = kernel.solve_flux(1, capped_b)
        self.assertTrue(adjacent_hit)
        self.assertFalse(phenotype_hit)
        self.assertEqual(kernel.cache_stats(), (2, 1, 2, 0, 2))

        kernel.clear_cache()
        unsaturated_a = np.full(3, 2.25 * tolerance, dtype=np.float64)
        unsaturated_b = np.full(3, 2.75 * tolerance, dtype=np.float64)
        first_flux, first_hit = kernel.solve_flux(0, unsaturated_a)
        cached_flux, within_tolerance_hit = kernel.solve_flux(0, unsaturated_b)
        self.assertFalse(first_hit)
        self.assertTrue(within_tolerance_hit)
        np.testing.assert_array_equal(cached_flux, first_flux)
        self.assertEqual(kernel.cache_stats(), (1, 1, 1, 0, 1))

        kernel.clear_cache()
        values = tuple(np.full(3, value, dtype=np.float64) for value in (0.1, 0.2, 0.3))
        kernel.solve_flux(0, values[0])
        kernel.solve_flux(0, values[1])
        kernel.solve_flux(0, values[0])
        kernel.solve_flux(0, values[2])
        self.assertEqual(kernel.cache_stats(), (2, 1, 3, 1, 3))
        _, evicted_hit = kernel.solve_flux(0, values[1])
        self.assertFalse(evicted_hit)
        self.assertEqual(kernel.cache_stats(), (2, 1, 4, 2, 4))

        kernel.clear_cache()
        kernel.solve_flux(0, values[0], cache_enabled=False)
        kernel.solve_flux(0, values[0], cache_enabled=False)
        self.assertEqual(kernel.cache_stats(), (0, 0, 0, 0, 2))

    def test_native_adaptive_batch_matches_python_and_conserves_mass(self):
        native = self.make_kernel()
        _, python_kernel = TestPythonMetabolicSolver().make_kernel()
        native_batch = self.make_batch()
        native_batch['environment_id'][:] = (0, 0, 1)
        python_batch = {name: value.copy() for name, value in native_batch.items()}
        initial_concentrations = native_batch['concentrations_m'].copy()

        native_summary = native.advance_batch(
            **native_batch,
            dt_seconds=0.75,
            depletion_safety=0.8,
            max_substeps=32,
            cache_enabled=True,
        )
        python_summary = python_kernel.advance_batch(
            **python_batch,
            dt_seconds=0.75,
            depletion_safety=0.8,
            max_substeps=32,
            cache_enabled=True,
        )

        self.assertEqual(native_summary[:5], python_summary[:5])
        self.assertGreater(native_summary[3], native_summary[1])
        self.assertAlmostEqual(native_summary[5], python_summary[5], places=12)
        for name in (
            'concentrations_m',
            'biomass_gdw',
            'growth_rate_s',
            'biomass_delta_gdw',
            'substrate_consumed_mol',
            'acetate_secreted_mol',
        ):
            np.testing.assert_allclose(
                native_batch[name], python_batch[name], rtol=0.0, atol=1.0e-12
            )
        for environment in range(2):
            members = native_batch['environment_id'] == environment
            depleted = (
                initial_concentrations[environment]
                - native_batch['concentrations_m'][environment]
            ) * native_batch['environment_volume_l'][environment]
            np.testing.assert_allclose(
                np.sum(native_batch['substrate_consumed_mol'][members], axis=0),
                depleted,
                rtol=0.0,
                atol=1.0e-18,
            )

    def test_native_batch_dt_zero_validation_alias_and_rollback(self):
        kernel = self.make_kernel()
        zero_batch = self.make_batch(agent_count=1, environment_count=1)
        summary = kernel.advance_batch(
            **zero_batch,
            dt_seconds=0.0,
            depletion_safety=0.8,
            max_substeps=2,
            cache_enabled=True,
        )
        self.assertEqual(summary, (0, 0, 0, 0, 0, 0.0))
        self.assertEqual(kernel.cache_stats(), (0, 0, 0, 0, 0))

        readonly_batch = self.make_batch(agent_count=1, environment_count=1)
        readonly_batch['concentrations_m'].setflags(write=False)
        with self.assertRaisesRegex(ValueError, 'concentrations_m.*writable'):
            kernel.advance_batch(
                **readonly_batch,
                dt_seconds=0.1,
                depletion_safety=0.8,
                max_substeps=2,
                cache_enabled=True,
            )

        alias_batch = self.make_batch(agent_count=1, environment_count=1)
        alias_batch['growth_rate_s'] = alias_batch['biomass_gdw']
        with self.assertRaisesRegex(ValueError, 'must not overlap'):
            kernel.advance_batch(
                **alias_batch,
                dt_seconds=0.1,
                depletion_safety=0.8,
                max_substeps=2,
                cache_enabled=True,
            )

        failing_batch = self.make_batch(agent_count=2, environment_count=1, concentration=5.0e-5)
        snapshots = {name: value.copy() for name, value in failing_batch.items()}
        with self.assertRaisesRegex(RuntimeError, 'max_substeps'):
            kernel.advance_batch(
                **failing_batch,
                dt_seconds=1.0,
                depletion_safety=0.5,
                max_substeps=1,
                cache_enabled=False,
            )
        for name, snapshot in snapshots.items():
            np.testing.assert_array_equal(failing_batch[name], snapshot)

    def test_native_solve_and_batch_reject_malformed_arrays_and_ranges(self):
        kernel = self.make_kernel()
        with self.assertRaisesRegex(TypeError, 'concentrations.*dtype'):
            kernel.solve_flux(0, np.ones(3, dtype=np.float32))
        with self.assertRaisesRegex(ValueError, 'C-contiguous'):
            kernel.solve_flux(0, np.ones(6, dtype=np.float64)[::2])
        with self.assertRaisesRegex(ValueError, 'nonnegative'):
            kernel.solve_flux(0, np.asarray((-1.0, 0.0, 0.0), dtype=np.float64))

        batch = self.make_batch(agent_count=1, environment_count=1)
        batch['environment_id'][0] = 1
        snapshots = {name: value.copy() for name, value in batch.items()}
        with self.assertRaisesRegex(ValueError, 'environment_id'):
            kernel.advance_batch(
                **batch,
                dt_seconds=0.1,
                depletion_safety=0.8,
                max_substeps=2,
                cache_enabled=True,
            )
        for name, snapshot in snapshots.items():
            np.testing.assert_array_equal(batch[name], snapshot)

    def test_cpp_orchestrator_shadow_restores_and_fba_commits(self):
        arguments = TestMetabolicOrchestrator().configured_arguments()
        for mode in ('shadow', 'fba'):
            with self.subTest(mode=mode):
                orchestrator = MetabolicOrchestrator(
                    mode=mode,
                    backend='cpp',
                    **arguments,
                )
                batch = self.make_batch(agent_count=1, environment_count=1, concentration=2.0e-4)
                initial_concentrations = batch['concentrations_m'].copy()
                initial_biomass = batch['biomass_gdw'].copy()
                orchestrator.advance_batch(
                    **batch,
                    dt_seconds=0.1,
                    depletion_safety=0.9,
                    max_substeps=8,
                    cache_enabled=True,
                )
                self.assertEqual(orchestrator.backend_used, 'cpp')
                self.assertTrue(orchestrator.native_available)
                self.assertTrue(np.all(batch['substrate_consumed_mol'] > 0.0))
                if mode == 'shadow':
                    np.testing.assert_array_equal(batch['concentrations_m'], initial_concentrations)
                    np.testing.assert_array_equal(batch['biomass_gdw'], initial_biomass)
                else:
                    self.assertTrue(np.all(batch['concentrations_m'] < initial_concentrations))
                    self.assertTrue(np.all(batch['biomass_gdw'] > initial_biomass))

    def test_cpp_species_batch_matches_python_capacity_gating(self):
        fixture = TestMetabolicOrchestrator()
        arguments = fixture.configured_species_arguments()
        python_batch = fixture.make_species_batch()
        native_batch = {name: value.copy() for name, value in python_batch.items()}
        python_orchestrator = MetabolicOrchestrator(
            mode='fba',
            backend='python',
            **arguments,
        )
        native_orchestrator = MetabolicOrchestrator(
            mode='fba',
            backend='cpp',
            **arguments,
        )

        for orchestrator, batch in (
            (python_orchestrator, python_batch),
            (native_orchestrator, native_batch),
        ):
            orchestrator.advance_batch(
                **batch,
                dt_seconds=0.1,
                depletion_safety=0.9,
                max_substeps=8,
                cache_enabled=False,
            )

        for name in (
            'concentrations_m',
            'biomass_gdw',
            'growth_rate_s',
            'biomass_delta_gdw',
            'substrate_consumed_mol',
            'acetate_secreted_mol',
        ):
            np.testing.assert_allclose(native_batch[name], python_batch[name], rtol=1.0e-12)
        self.assertEqual(native_batch['substrate_consumed_mol'][1, 2], 0.0)


if __name__ == '__main__':
    unittest.main()