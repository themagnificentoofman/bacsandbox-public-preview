import gc
import importlib
import itertools
import os
import statistics
import time
import tracemalloc
import unittest
import weakref

import numpy as np

from bacteria.metabolism import (
    CACHE_TOLERANCE_MOLAR,
    MAX_CACHE_CAPACITY,
    PythonMetabolicKernel,
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
RUN_METABOLIC_SCALE = os.environ.get('BACSANDBOX_RUN_METABOLIC_SCALE') == '1'
RUN_FBA_BENCHMARK = os.environ.get('BACSANDBOX_RUN_FBA_BENCHMARK') == '1'
NATIVE_METABOLIC_STRESS_DEFAULT_ITERATIONS = 100
MUTABLE_BATCH_ARRAYS = (
    'concentrations_m',
    'biomass_gdw',
    'growth_rate_s',
    'biomass_delta_gdw',
    'substrate_consumed_mol',
    'acetate_secreted_mol',
)


def available_kernel_types():
    kernel_types = [('python', PythonMetabolicKernel)]
    if NATIVE_METABOLIC_AVAILABLE:
        kernel_types.append(('native', bacteria_c_kernel.MetabolicKernel))
    return tuple(kernel_types)


def kernel_arguments(cache_capacity=8):
    model = load_metabolic_model()
    return model, {
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


def configure_kernel(
    kernel,
    model,
    phenotype_count=1,
    model_revision=None,
    prior_revision=1,
):
    # Synthetic capacity 2 M/gDW/s and half-saturation 1 M keep the reduced
    # fixture unsaturated at micromolar inputs; see DOCUMENTATION.md §10.3.
    kernel.configure_phenotypes(
        np.zeros(phenotype_count, dtype=np.int32),
        np.full((phenotype_count, 3), 2.0, dtype=np.float64),
        np.ones((phenotype_count, 3), dtype=np.float64),
        np.ones((phenotype_count, 3), dtype=np.float64),
        np.ones((phenotype_count, 3), dtype=np.float64),
        np.zeros((phenotype_count, 3), dtype=np.float64),
        model.model_revision if model_revision is None else model_revision,
        prior_revision,
    )


def make_kernel(kernel_type, cache_capacity=8, phenotype_count=1):
    model, arguments = kernel_arguments(cache_capacity)
    kernel = kernel_type(**arguments)
    configure_kernel(kernel, model, phenotype_count)
    return model, kernel


def make_batch(agent_count=2, environment_count=1, concentration_m=5.0e-5):
    # A 50 uM shared pool, 1 uL volume, and 0.5 ugDW agents force more than one
    # depletion-limited substep without negative mass; see DOCUMENTATION.md §10.3.
    return {
        'phenotype_id': np.zeros(agent_count, dtype=np.int32),
        'environment_id': np.arange(agent_count, dtype=np.int32) % environment_count,
        'concentrations_m': np.full(
            (environment_count, 3), concentration_m, dtype=np.float64
        ),
        'environment_volume_l': np.full(environment_count, 1.0e-6, dtype=np.float64),
        'biomass_gdw': np.full(agent_count, 5.0e-7, dtype=np.float64),
        'growth_rate_s': np.full(agent_count, -1.0, dtype=np.float64),
        'biomass_delta_gdw': np.full(agent_count, -2.0, dtype=np.float64),
        'substrate_consumed_mol': np.full((agent_count, 3), -3.0, dtype=np.float64),
        'acetate_secreted_mol': np.full(agent_count, -4.0, dtype=np.float64),
    }


def copy_batch(batch):
    return {name: np.array(value, copy=True) for name, value in batch.items()}


def dense_stoichiometry(model):
    matrix = np.zeros(
        (int(model.model_metabolite_offsets[-1]), model.lower_bounds.size),
        dtype=np.float64,
    )
    for row in range(matrix.shape[0]):
        start = int(model.stoich_row_offsets[row])
        stop = int(model.stoich_row_offsets[row + 1])
        matrix[row, model.stoich_column_indices[start:stop]] = model.stoich_values[start:stop]
    return matrix


def expected_reduced_flux(concentrations_m, inhibition_multiplier=1.0):
    concentrations = np.asarray(concentrations_m, dtype=np.float64)
    # Michaelis-Menten uptake 2C/(1+C), posterior inhibition, and the finite unit
    # exchange cap define each synthetic reference flux; see DOCUMENTATION.md §10.3.
    uptake = np.minimum(
        1.0,
        2.0 * concentrations / (1.0 + concentrations) * inhibition_multiplier,
    )
    # The reduced fixture routes every uptake unit to biomass and exports one-quarter
    # acetate bookkeeping unit, so Sx remains zero; see DOCUMENTATION.md §10.3.
    return np.concatenate((uptake, uptake, (np.sum(uptake), 0.25 * np.sum(uptake))))


def configure_neutral_and_inhibitory(kernel, model):
    # Phenotype 0 is neutral; phenotype 1 uses alpha=2, beta=3, and a 0.5 cap,
    # giving active fraction 0.8 without stochastic sampling; see DOCUMENTATION.md §10.3.
    kernel.configure_phenotypes(
        np.zeros(2, dtype=np.int32),
        np.full((2, 3), 2.0, dtype=np.float64),
        np.ones((2, 3), dtype=np.float64),
        np.asarray(((1.0, 1.0, 1.0), (2.0, 2.0, 2.0)), dtype=np.float64),
        np.asarray(((1.0, 1.0, 1.0), (3.0, 3.0, 3.0)), dtype=np.float64),
        np.asarray(((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)), dtype=np.float64),
        model.model_revision,
        7,
    )


def degenerate_model_arguments(cache_capacity=8):
    # Two equivalent generic-calorie routes create a bounded alternate optimum;
    # fixed pivot ordering must choose one reproducibly; see DOCUMENTATION.md §10.3.
    stoichiometry = np.asarray((
        (1.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, -1.0, 0.0),
        (0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25, 0.0, -1.0),
    ), dtype=np.float64)
    row_offsets = [0]
    column_indices = []
    values = []
    for row in stoichiometry:
        columns = np.flatnonzero(row)
        column_indices.extend(int(column) for column in columns)
        values.extend(float(row[column]) for column in columns)
        row_offsets.append(len(values))
    return {
        'model_reaction_offsets': np.asarray((0, 9), dtype=np.int64),
        'model_metabolite_offsets': np.asarray((0, 5), dtype=np.int64),
        'stoich_row_offsets': np.asarray(row_offsets, dtype=np.int64),
        'stoich_column_indices': np.asarray(column_indices, dtype=np.int32),
        'stoich_values': np.asarray(values, dtype=np.float64),
        'lower_bounds': np.zeros(9, dtype=np.float64),
        'upper_bounds': np.asarray((1.0,) * 7 + (3.0, 3.0), dtype=np.float64),
        'objective_coefficients': np.asarray((0.0,) * 7 + (1.0, 0.0), dtype=np.float64),
        'exchange_reaction_indices': np.asarray(((0, 1, 2),), dtype=np.int32),
        'biomass_reaction_indices': np.asarray((7,), dtype=np.int32),
        'acetate_reaction_indices': np.asarray((8,), dtype=np.int32),
        'cache_tolerance_molar': CACHE_TOLERANCE_MOLAR,
        'cache_capacity': cache_capacity,
    }


class TestMetabolicReferenceSolutions(unittest.TestCase):
    def test_reduced_model_ownership_and_neutral_production_priors(self):
        model = load_metabolic_model()
        priors = load_metabolic_priors()

        self.assertEqual(model.provenance, 'Own (BacSandbox model)')
        self.assertIn('not a genome-scale model', model.scope_statement.lower())
        self.assertEqual(model.model_count, 1)
        self.assertEqual(model.lower_bounds.size, 8)
        self.assertEqual(int(model.model_metabolite_offsets[-1]), 5)
        self.assertTrue(all(row.max_inhibition == 0.0 for row in priors.rows))

    def test_individual_and_mixed_neutral_and_inhibitory_flux_vectors(self):
        # One-third molar gives half of the unit exchange cap with the synthetic
        # 2C/(1+C) uptake law; see DOCUMENTATION.md §10.3.
        cases = {
            'generic': np.asarray((1.0 / 3.0, 0.0, 0.0), dtype=np.float64),
            'glycan': np.asarray((0.0, 1.0 / 3.0, 0.0), dtype=np.float64),
            'mdcf': np.asarray((0.0, 0.0, 1.0 / 3.0), dtype=np.float64),
            'mixed': np.asarray((1.0, 1.0 / 3.0, 0.1), dtype=np.float64),
        }
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                model, kernel = make_kernel(kernel_type)
                configure_neutral_and_inhibitory(kernel, model)
                stoichiometry = dense_stoichiometry(model)
                for case_name, concentrations in cases.items():
                    for phenotype, multiplier in ((0, 1.0), (1, 0.8)):
                        with self.subTest(case=case_name, phenotype=phenotype):
                            flux, cache_hit = kernel.solve_flux(
                                phenotype, concentrations, cache_enabled=False
                            )
                            expected = expected_reduced_flux(concentrations, multiplier)
                            self.assertFalse(cache_hit)
                            np.testing.assert_allclose(flux, expected, rtol=0.0, atol=1.0e-10)
                            np.testing.assert_allclose(flux[:3], expected[:3], rtol=0.0, atol=1.0e-10)
                            self.assertAlmostEqual(flux[6], expected[6], delta=1.0e-10)
                            self.assertAlmostEqual(flux[7], expected[7], delta=1.0e-10)
                            self.assertLessEqual(
                                float(np.max(np.abs(stoichiometry @ flux))),
                                1.0e-10,
                            )

    def test_seeded_120_case_python_native_reference_sweep(self):
        if not NATIVE_METABOLIC_AVAILABLE:
            self.skipTest('native MetabolicKernel is not built')
        python_model, python_kernel = make_kernel(PythonMetabolicKernel)
        native_model, native_kernel = make_kernel(bacteria_c_kernel.MetabolicKernel)
        configure_neutral_and_inhibitory(python_kernel, python_model)
        configure_neutral_and_inhibitory(native_kernel, native_model)
        tolerance = CACHE_TOLERANCE_MOLAR
        threshold = 1.0e-4
        # Fixed values straddle cache and adaptive boundaries, while the seeded log
        # sweep reaches zero, saturation, and the finite model cap; see DOCUMENTATION.md §10.3.
        edge_values = np.asarray((
            0.0,
            0.5 * tolerance,
            0.999 * tolerance,
            1.001 * tolerance,
            threshold - 0.1 * tolerance,
            threshold,
            threshold + 0.1 * tolerance,
            1.0 / 3.0,
            1.0,
            1.0e6,
        ), dtype=np.float64)
        fixed_cases = [
            np.asarray((edge_values[index % 10], edge_values[(index + 3) % 10], edge_values[(index + 7) % 10]))
            for index in range(20)
        ]
        random = np.random.default_rng(3103)
        random_cases = 10.0 ** random.uniform(-9.0, 6.0, size=(100, 3))
        random_cases[random.random((100, 3)) < 0.08] = 0.0
        cases = fixed_cases + [row for row in random_cases]
        self.assertEqual(len(cases), 120)

        for index, concentrations in enumerate(cases):
            phenotype = index % 2
            multiplier = 1.0 if phenotype == 0 else 0.8
            python_flux, python_hit = python_kernel.solve_flux(
                phenotype, concentrations, cache_enabled=False
            )
            native_flux, native_hit = native_kernel.solve_flux(
                phenotype, concentrations, cache_enabled=False
            )
            expected = expected_reduced_flux(concentrations, multiplier)
            self.assertFalse(python_hit)
            self.assertEqual(native_hit, python_hit)
            np.testing.assert_allclose(python_flux, expected, rtol=0.0, atol=1.0e-9)
            np.testing.assert_allclose(native_flux, expected, rtol=0.0, atol=1.0e-9)
            np.testing.assert_allclose(native_flux, python_flux, rtol=0.0, atol=1.0e-10)

        self.assertEqual(python_kernel.cache_stats(), (0, 0, 0, 0, 120))
        self.assertEqual(native_kernel.cache_stats(), python_kernel.cache_stats())

    def test_degenerate_pivot_is_repeatable_and_matches_native(self):
        results = {}
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                kernel = kernel_type(**degenerate_model_arguments())
                model = load_metabolic_model()
                configure_kernel(kernel, model)
                concentration = np.asarray((1.0, 0.0, 0.0), dtype=np.float64)
                repeated = [
                    kernel.solve_flux(0, concentration, cache_enabled=False)[0]
                    for _ in range(4)
                ]
                for flux in repeated[1:]:
                    np.testing.assert_array_equal(flux, repeated[0])
                # Bland-style index ordering deterministically selects route 0 among
                # equivalent bounded routes; see DOCUMENTATION.md §10.3.
                np.testing.assert_allclose(repeated[0][3:5], (1.0, 0.0), rtol=0.0, atol=1.0e-10)
                results[backend] = repeated[0]
        if 'native' in results:
            np.testing.assert_allclose(results['native'], results['python'], rtol=0.0, atol=1.0e-10)

    def test_explicit_infeasible_and_nonfinite_upper_bound_rejection(self):
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                model, arguments = kernel_arguments()
                arguments['lower_bounds'] = np.array(arguments['lower_bounds'], copy=True)
                arguments['lower_bounds'][0] = 0.5
                infeasible = kernel_type(**arguments)
                configure_kernel(infeasible, model)
                with self.assertRaisesRegex(RuntimeError, 'infeasible'):
                    infeasible.solve_flux(0, np.zeros(3, dtype=np.float64), False)

                # Finite reaction bounds make unbounded models unsupported by contract;
                # rejecting infinity is intentional; see DOCUMENTATION.md §10.3.
                nonfinite = kernel_arguments()[1]
                nonfinite['upper_bounds'] = np.array(nonfinite['upper_bounds'], copy=True)
                nonfinite['upper_bounds'][0] = np.inf
                with self.assertRaisesRegex(ValueError, 'finite'):
                    kernel_type(**nonfinite)


class TestMetabolicCacheCorrectness(unittest.TestCase):
    def test_unsaturated_same_adjacent_and_over_tolerance_bins(self):
        tolerance = CACHE_TOLERANCE_MOLAR
        concentrations = (
            np.full(3, 2.25 * tolerance, dtype=np.float64),
            np.full(3, 2.75 * tolerance, dtype=np.float64),
            np.full(3, 3.10 * tolerance, dtype=np.float64),
            np.full(3, 3.250001 * tolerance, dtype=np.float64),
        )
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                _, kernel = make_kernel(kernel_type)
                first_flux, first_hit = kernel.solve_flux(0, concentrations[0])
                same_flux, same_hit = kernel.solve_flux(0, concentrations[1])
                adjacent_flux, adjacent_hit = kernel.solve_flux(0, concentrations[2])
                over_flux, over_hit = kernel.solve_flux(0, concentrations[3])

                self.assertFalse(first_hit)
                self.assertTrue(same_hit)
                self.assertTrue(adjacent_hit)
                self.assertFalse(over_hit)
                np.testing.assert_array_equal(same_flux, first_flux)
                np.testing.assert_array_equal(adjacent_flux, first_flux)
                self.assertFalse(np.array_equal(over_flux, first_flux))
                self.assertEqual(kernel.cache_stats(), (2, 2, 2, 0, 2))

    def test_phenotype_and_revision_invalidation(self):
        concentration = np.full(3, 0.2, dtype=np.float64)
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                model, kernel = make_kernel(kernel_type, phenotype_count=2)
                kernel.solve_flux(0, concentration)
                _, cross_phenotype_hit = kernel.solve_flux(1, concentration)
                self.assertFalse(cross_phenotype_hit)
                self.assertEqual(kernel.cache_stats(), (2, 0, 2, 0, 2))

                configure_kernel(
                    kernel,
                    model,
                    phenotype_count=2,
                    model_revision=model.model_revision + 1,
                    prior_revision=1,
                )
                self.assertEqual(kernel.cache_stats(), (0, 0, 0, 0, 0))
                kernel.solve_flux(0, concentration)
                configure_kernel(
                    kernel,
                    model,
                    phenotype_count=2,
                    model_revision=model.model_revision + 1,
                    prior_revision=2,
                )
                self.assertEqual(kernel.cache_stats(), (0, 0, 0, 0, 0))

    def test_python_lru_tie_uses_oldest_insertion(self):
        _, kernel = make_kernel(PythonMetabolicKernel, cache_capacity=2)
        concentrations = [
            np.full(3, value, dtype=np.float64) for value in (0.1, 0.2, 0.3)
        ]
        kernel.solve_flux(0, concentrations[0])
        kernel.solve_flux(0, concentrations[1])
        entries = sorted(kernel._cache.items(), key=lambda item: item[1].insertion_sequence)
        oldest_key = entries[0][0]
        newest_key = entries[1][0]
        for entry in kernel._cache.values():
            entry.last_used_sequence = 100
        kernel.solve_flux(0, concentrations[2])

        self.assertNotIn(oldest_key, kernel._cache)
        self.assertIn(newest_key, kernel._cache)
        self.assertEqual(kernel.cache_stats(), (2, 0, 3, 1, 3))

    def test_capacity_one_high_churn_and_disabled_cache(self):
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                _, kernel = make_kernel(kernel_type, cache_capacity=1)
                # Four cache bins between requests force deterministic capacity-one
                # eviction without tolerance overlap; see DOCUMENTATION.md §10.3.
                states = [
                    np.full(3, 0.1 + 4.0 * CACHE_TOLERANCE_MOLAR * index, dtype=np.float64)
                    for index in range(10)
                ]
                for state in states:
                    _, cache_hit = kernel.solve_flux(0, state)
                    self.assertFalse(cache_hit)
                _, final_hit = kernel.solve_flux(0, states[-1])
                self.assertTrue(final_hit)
                self.assertEqual(kernel.cache_stats(), (1, 1, 10, 9, 10))

                kernel.clear_cache()
                for _ in range(3):
                    _, cache_hit = kernel.solve_flux(0, states[0], cache_enabled=False)
                    self.assertFalse(cache_hit)
                self.assertEqual(kernel.cache_stats(), (0, 0, 0, 0, 3))

    def test_batch_dedup_disabled_parity_and_warmed_hits(self):
        agent_count = 96
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                _, cached_kernel = make_kernel(kernel_type, cache_capacity=8)
                _, uncached_kernel = make_kernel(kernel_type, cache_capacity=8)
                # A large environment and tiny biomass keep this one-step fixture away
                # from depletion, isolating batch cache accounting; see DOCUMENTATION.md §10.3.
                initial = make_batch(agent_count=agent_count, concentration_m=0.25)
                initial['environment_volume_l'][:] = 1.0e-3
                initial['biomass_gdw'][:] = 1.0e-9
                cached_batch = copy_batch(initial)
                uncached_batch = copy_batch(initial)
                common = {
                    'dt_seconds': 0.01,
                    'low_substrate_threshold_m': 1.0e-4,
                    'depletion_safety': 0.9,
                    'max_substeps': 8,
                }
                cached_summary = cached_kernel.advance_batch(
                    **cached_batch, **common, cache_enabled=True
                )
                uncached_summary = uncached_kernel.advance_batch(
                    **uncached_batch, **common, cache_enabled=False
                )

                self.assertEqual(cached_summary[:4], (1, 1, 0, agent_count))
                self.assertEqual(uncached_summary[:4], (1, agent_count, 0, 0))
                for name in MUTABLE_BATCH_ARRAYS:
                    np.testing.assert_allclose(
                        cached_batch[name], uncached_batch[name], rtol=0.0, atol=1.0e-12
                    )

                warmed_batch = copy_batch(initial)
                warmed_summary = cached_kernel.advance_batch(
                    **warmed_batch, **common, cache_enabled=True
                )
                self.assertEqual(warmed_summary[:4], (1, 0, agent_count, 0))
                self.assertEqual(cached_kernel.cache_stats(), (1, agent_count, agent_count, 0, 1))


class TestMetabolicAdaptiveIntegration(unittest.TestCase):
    def test_glycan_mdcf_threshold_is_strictly_less_than(self):
        threshold = 1.0e-4
        cases = (
            ('above', 1.1e-4, 1),
            ('exact', threshold, 1),
            ('below', 0.9e-4, None),
        )
        backend_results = {}
        for backend, kernel_type in available_kernel_types():
            backend_results[backend] = {}
            for name, concentration, expected_substeps in cases:
                with self.subTest(backend=backend, boundary=name):
                    _, kernel = make_kernel(kernel_type)
                    batch = make_batch(agent_count=1, concentration_m=concentration)
                    # Generic calorie stays abundant so glycan/MDCF alone determine
                    # whether the 0.1 mM safety threshold applies; see DOCUMENTATION.md §10.3.
                    batch['concentrations_m'][0, 0] = 1.0
                    summary = kernel.advance_batch(
                        **batch,
                        dt_seconds=0.9,
                        low_substrate_threshold_m=threshold,
                        depletion_safety=0.8,
                        max_substeps=8,
                        cache_enabled=False,
                    )
                    if expected_substeps is None:
                        self.assertGreater(summary[0], 1)
                        self.assertLess(summary[5], 0.9)
                    else:
                        self.assertEqual(summary[0], expected_substeps)
                        self.assertEqual(summary[5], 0.9)
                    backend_results[backend][name] = (summary, batch)

        if 'native' in backend_results:
            for name, _, _ in cases:
                python_summary, python_batch = backend_results['python'][name]
                native_summary, native_batch = backend_results['native'][name]
                self.assertEqual(native_summary[:5], python_summary[:5])
                self.assertAlmostEqual(native_summary[5], python_summary[5], delta=1.0e-12)
                for field in MUTABLE_BATCH_ARRAYS:
                    np.testing.assert_allclose(
                        native_batch[field], python_batch[field], rtol=0.0, atol=1.0e-12
                    )

    def test_thousand_agent_sudden_nutritional_crash_conserves_mass(self):
        agent_count = 1000
        environment_count = 4
        initial = make_batch(
            agent_count=agent_count,
            environment_count=environment_count,
            concentration_m=2.0e-4,
        )
        # Dense 0.5 ugDW agents in 1 uL pools exhaust 0.2 mM substrate early in a
        # one-second tick, exercising the crash substep; see DOCUMENTATION.md §10.3.
        initial['biomass_gdw'][:] = 5.0e-7
        backend_results = {}
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                _, kernel = make_kernel(kernel_type, cache_capacity=16)
                batch = copy_batch(initial)
                starting_concentrations = batch['concentrations_m'].copy()
                summary = kernel.advance_batch(
                    **batch,
                    dt_seconds=1.0,
                    low_substrate_threshold_m=1.0e-4,
                    depletion_safety=0.9,
                    max_substeps=8,
                    cache_enabled=True,
                )

                self.assertEqual(summary[0], 2)
                self.assertGreater(summary[5], 0.0)
                self.assertLess(summary[5], 1.0)
                self.assertTrue(np.all(batch['concentrations_m'] >= 0.0))
                self.assertTrue(np.all(np.isfinite(batch['biomass_gdw'])))
                for environment in range(environment_count):
                    members = batch['environment_id'] == environment
                    depleted_mol = (
                        starting_concentrations[environment]
                        - batch['concentrations_m'][environment]
                    ) * batch['environment_volume_l'][environment]
                    np.testing.assert_allclose(
                        np.sum(batch['substrate_consumed_mol'][members], axis=0),
                        depleted_mol,
                        rtol=0.0,
                        atol=1.0e-15,
                    )
                backend_results[backend] = (summary, batch)

        if 'native' in backend_results:
            python_summary, python_batch = backend_results['python']
            native_summary, native_batch = backend_results['native']
            self.assertEqual(native_summary[:5], python_summary[:5])
            self.assertAlmostEqual(native_summary[5], python_summary[5], delta=1.0e-12)
            for field in MUTABLE_BATCH_ARRAYS:
                np.testing.assert_allclose(
                    native_batch[field], python_batch[field], rtol=0.0, atol=1.0e-5
                )

    def test_multiple_environments_zero_concentration_and_zero_biomass(self):
        initial = make_batch(agent_count=6, environment_count=3, concentration_m=2.0e-4)
        initial['concentrations_m'][0] = 0.0
        initial['concentrations_m'][1] = (1.0e-4, 0.0, 2.0e-4)
        initial['biomass_gdw'][::2] = 0.0
        backend_results = {}
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                _, kernel = make_kernel(kernel_type)
                batch = copy_batch(initial)
                starting_concentrations = batch['concentrations_m'].copy()
                starting_biomass = batch['biomass_gdw'].copy()
                summary = kernel.advance_batch(
                    **batch,
                    dt_seconds=0.2,
                    low_substrate_threshold_m=1.0e-4,
                    depletion_safety=0.8,
                    max_substeps=16,
                    cache_enabled=True,
                )
                self.assertGreaterEqual(summary[0], 1)
                np.testing.assert_array_equal(batch['concentrations_m'][0], 0.0)
                np.testing.assert_array_equal(batch['biomass_gdw'][::2], starting_biomass[::2])
                np.testing.assert_array_equal(batch['biomass_delta_gdw'][::2], 0.0)
                np.testing.assert_array_equal(batch['substrate_consumed_mol'][::2], 0.0)
                for environment in range(3):
                    members = batch['environment_id'] == environment
                    depleted_mol = (
                        starting_concentrations[environment]
                        - batch['concentrations_m'][environment]
                    ) * batch['environment_volume_l'][environment]
                    np.testing.assert_allclose(
                        np.sum(batch['substrate_consumed_mol'][members], axis=0),
                        depleted_mol,
                        rtol=0.0,
                        atol=1.0e-18,
                    )
                backend_results[backend] = (summary, batch)

        if 'native' in backend_results:
            self.assertEqual(backend_results['native'][0][:5], backend_results['python'][0][:5])
            for field in MUTABLE_BATCH_ARRAYS:
                np.testing.assert_allclose(
                    backend_results['native'][1][field],
                    backend_results['python'][1][field],
                    rtol=0.0,
                    atol=1.0e-12,
                )

    def test_outputs_reset_per_call_and_accumulate_across_substeps(self):
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                _, kernel = make_kernel(kernel_type)
                batch = make_batch(agent_count=3, concentration_m=8.0e-5)
                initial_concentrations = batch['concentrations_m'].copy()
                initial_biomass = batch['biomass_gdw'].copy()
                summary = kernel.advance_batch(
                    **batch,
                    dt_seconds=0.75,
                    low_substrate_threshold_m=1.0e-4,
                    depletion_safety=0.8,
                    max_substeps=32,
                    cache_enabled=True,
                )
                self.assertGreater(summary[0], 1)
                np.testing.assert_allclose(
                    batch['biomass_delta_gdw'],
                    batch['biomass_gdw'] - initial_biomass,
                    rtol=0.0,
                    atol=1.0e-18,
                )
                depleted_mol = (
                    initial_concentrations[0] - batch['concentrations_m'][0]
                ) * batch['environment_volume_l'][0]
                np.testing.assert_allclose(
                    np.sum(batch['substrate_consumed_mol'], axis=0),
                    depleted_mol,
                    rtol=0.0,
                    atol=1.0e-18,
                )

                second_initial_biomass = batch['biomass_gdw'].copy()
                batch['growth_rate_s'].fill(91.0)
                batch['biomass_delta_gdw'].fill(92.0)
                batch['substrate_consumed_mol'].fill(93.0)
                batch['acetate_secreted_mol'].fill(94.0)
                kernel.advance_batch(
                    **batch,
                    dt_seconds=0.01,
                    low_substrate_threshold_m=1.0e-4,
                    depletion_safety=0.8,
                    max_substeps=8,
                    cache_enabled=True,
                )
                np.testing.assert_allclose(
                    batch['biomass_delta_gdw'],
                    batch['biomass_gdw'] - second_initial_biomass,
                    rtol=0.0,
                    atol=1.0e-18,
                )
                self.assertTrue(np.all(batch['growth_rate_s'] < 91.0))
                self.assertTrue(np.all(batch['substrate_consumed_mol'] < 93.0))
                self.assertTrue(np.all(batch['acetate_secreted_mol'] < 94.0))

    def test_dt_zero_and_internal_substeps_do_not_mutate_input_scalar(self):
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                _, kernel = make_kernel(kernel_type)
                zero_batch = make_batch(agent_count=2)
                state_snapshots = {
                    name: zero_batch[name].copy()
                    for name in ('concentrations_m', 'biomass_gdw')
                }
                summary = kernel.advance_batch(
                    **zero_batch,
                    dt_seconds=0.0,
                    low_substrate_threshold_m=1.0e-4,
                    depletion_safety=0.8,
                    max_substeps=8,
                    cache_enabled=True,
                )
                self.assertEqual(summary, (0, 0, 0, 0, 0, 0.0))
                for name, expected in state_snapshots.items():
                    np.testing.assert_array_equal(zero_batch[name], expected)
                for name in (
                    'growth_rate_s',
                    'biomass_delta_gdw',
                    'substrate_consumed_mol',
                    'acetate_secreted_mol',
                ):
                    np.testing.assert_array_equal(zero_batch[name], 0.0)

                adaptive_batch = make_batch(agent_count=1, concentration_m=0.9e-4)
                adaptive_batch['concentrations_m'][0, 0] = 1.0
                global_dt_seconds = 0.9
                adaptive_summary = kernel.advance_batch(
                    **adaptive_batch,
                    dt_seconds=global_dt_seconds,
                    low_substrate_threshold_m=1.0e-4,
                    depletion_safety=0.8,
                    max_substeps=8,
                    cache_enabled=False,
                )
                self.assertGreater(adaptive_summary[0], 1)
                self.assertEqual(global_dt_seconds, 0.9)


class MetabolicBatchValidationMixin:
    def assert_rejected_without_writes(
        self,
        kernel,
        batch,
        overrides=None,
        exception_types=(TypeError, ValueError, OverflowError),
    ):
        snapshots = {name: batch[name].tobytes() for name in MUTABLE_BATCH_ARRAYS}
        arguments = {
            'dt_seconds': 0.1,
            'low_substrate_threshold_m': 1.0e-4,
            'depletion_safety': 0.8,
            'max_substeps': 8,
            'cache_enabled': True,
        }
        if overrides:
            arguments.update(overrides)
        with self.assertRaises(exception_types):
            kernel.advance_batch(**batch, **arguments)
        for name, expected in snapshots.items():
            self.assertEqual(batch[name].tobytes(), expected, name)


class TestMetabolicBatchABI(MetabolicBatchValidationMixin, unittest.TestCase):
    def test_exact_dtypes_and_dimensions_are_table_driven(self):
        for backend, kernel_type in available_kernel_types():
            _, kernel = make_kernel(kernel_type)
            template = make_batch()
            for name, original in template.items():
                with self.subTest(backend=backend, field=name, contract='dtype'):
                    batch = copy_batch(template)
                    wrong_dtype = np.float32 if original.dtype == np.float64 else np.int64
                    batch[name] = original.astype(wrong_dtype)
                    self.assert_rejected_without_writes(kernel, batch)

                with self.subTest(backend=backend, field=name, contract='dimension'):
                    batch = copy_batch(template)
                    if original.ndim == 1:
                        batch[name] = original.reshape((1, original.size))
                    else:
                        batch[name] = original.reshape(-1)
                    self.assert_rejected_without_writes(kernel, batch)

    def test_noncontiguous_arrays_and_n_e_p_mismatches(self):
        for backend, kernel_type in available_kernel_types():
            _, kernel = make_kernel(kernel_type)
            template = make_batch(environment_count=2)
            for name, original in template.items():
                with self.subTest(backend=backend, field=name, contract='C-contiguous'):
                    expanded_shape = list(original.shape)
                    expanded_shape[-1] *= 2
                    expanded = np.zeros(expanded_shape, dtype=original.dtype)
                    slices = (slice(None),) * (original.ndim - 1) + (slice(None, None, 2),)
                    batch = copy_batch(template)
                    batch[name] = expanded[slices]
                    self.assertFalse(batch[name].flags.c_contiguous)
                    self.assert_rejected_without_writes(kernel, batch)

            mismatch_batches = []
            batch = copy_batch(template)
            batch['environment_id'] = np.zeros(3, dtype=np.int32)
            mismatch_batches.append(('N environment_id', batch))
            batch = copy_batch(template)
            batch['environment_volume_l'] = np.ones(3, dtype=np.float64)
            mismatch_batches.append(('E volume', batch))
            batch = copy_batch(template)
            batch['substrate_consumed_mol'] = np.zeros((3, 3), dtype=np.float64)
            mismatch_batches.append(('N consumed', batch))
            for name, batch in mismatch_batches:
                with self.subTest(backend=backend, mismatch=name):
                    self.assert_rejected_without_writes(kernel, batch)

            model, fresh = make_kernel(kernel_type)
            with self.subTest(backend=backend, mismatch='P configuration'):
                with self.assertRaises((TypeError, ValueError)):
                    fresh.configure_phenotypes(
                        np.zeros(2, dtype=np.int32),
                        np.ones((1, 3), dtype=np.float64),
                        np.ones((2, 3), dtype=np.float64),
                        np.ones((2, 3), dtype=np.float64),
                        np.ones((2, 3), dtype=np.float64),
                        np.zeros((2, 3), dtype=np.float64),
                        model.model_revision,
                        1,
                    )

    def test_all_six_mutable_arrays_require_writable_storage(self):
        for backend, kernel_type in available_kernel_types():
            _, kernel = make_kernel(kernel_type)
            for name in MUTABLE_BATCH_ARRAYS:
                with self.subTest(backend=backend, field=name):
                    batch = make_batch()
                    batch[name].setflags(write=False)
                    self.assert_rejected_without_writes(kernel, batch)

    def test_every_unsafe_array_pair_and_partial_overlap_is_rejected(self):
        template = make_batch()
        names = tuple(template)
        pairs = (
            pair
            for pair in itertools.combinations(names, 2)
            if pair[0] in MUTABLE_BATCH_ARRAYS or pair[1] in MUTABLE_BATCH_ARRAYS
        )
        pairs = tuple(pairs)
        self.assertEqual(len(pairs), 33)
        for backend, kernel_type in available_kernel_types():
            _, kernel = make_kernel(kernel_type)
            for left_name, right_name in pairs:
                with self.subTest(backend=backend, left=left_name, right=right_name):
                    batch = copy_batch(template)
                    left_source = batch[left_name]
                    right_source = batch[right_name]
                    shared = np.zeros(max(left_source.nbytes, right_source.nbytes), dtype=np.uint8)
                    batch[left_name] = np.ndarray(
                        left_source.shape, dtype=left_source.dtype, buffer=shared
                    )
                    batch[right_name] = np.ndarray(
                        right_source.shape, dtype=right_source.dtype, buffer=shared
                    )
                    self.assertTrue(np.shares_memory(batch[left_name], batch[right_name]))
                    self.assert_rejected_without_writes(kernel, batch)

            with self.subTest(backend=backend, overlap='partial'):
                batch = make_batch()
                shared = np.ones(3, dtype=np.float64)
                batch['biomass_gdw'] = shared[:2]
                batch['growth_rate_s'] = shared[1:]
                self.assertTrue(np.shares_memory(batch['biomass_gdw'], batch['growth_rate_s']))
                self.assert_rejected_without_writes(kernel, batch)

    def test_invalid_ids_nonfinite_and_negative_values_validate_before_write(self):
        cases = (
            ('phenotype negative', 'phenotype_id', 0, -1),
            ('phenotype high', 'phenotype_id', 0, 1),
            ('environment negative', 'environment_id', 0, -1),
            ('environment high', 'environment_id', 0, 1),
            ('concentration negative', 'concentrations_m', (0, 0), -1.0),
            ('concentration nan', 'concentrations_m', (0, 0), np.nan),
            ('concentration inf', 'concentrations_m', (0, 0), np.inf),
            ('volume zero', 'environment_volume_l', 0, 0.0),
            ('volume negative', 'environment_volume_l', 0, -1.0),
            ('volume nan', 'environment_volume_l', 0, np.nan),
            ('volume inf', 'environment_volume_l', 0, np.inf),
            ('biomass negative', 'biomass_gdw', 0, -1.0),
            ('biomass nan', 'biomass_gdw', 0, np.nan),
            ('biomass inf', 'biomass_gdw', 0, np.inf),
            ('growth nan', 'growth_rate_s', 0, np.nan),
            ('growth inf', 'growth_rate_s', 0, np.inf),
            ('delta nan', 'biomass_delta_gdw', 0, np.nan),
            ('consumed inf', 'substrate_consumed_mol', (0, 0), np.inf),
            ('acetate nan', 'acetate_secreted_mol', 0, np.nan),
        )
        for backend, kernel_type in available_kernel_types():
            _, kernel = make_kernel(kernel_type)
            for case_name, field, index, value in cases:
                with self.subTest(backend=backend, case=case_name):
                    batch = make_batch()
                    batch[field][index] = value
                    self.assert_rejected_without_writes(kernel, batch)

    def test_invalid_scalars_validate_before_write(self):
        cases = (
            ('dt negative', {'dt_seconds': -1.0}),
            ('dt nan', {'dt_seconds': np.nan}),
            ('dt inf', {'dt_seconds': np.inf}),
            ('threshold negative', {'low_substrate_threshold_m': -1.0}),
            ('threshold nan', {'low_substrate_threshold_m': np.nan}),
            ('threshold inf', {'low_substrate_threshold_m': np.inf}),
            ('safety zero', {'depletion_safety': 0.0}),
            ('safety negative', {'depletion_safety': -0.1}),
            ('safety high', {'depletion_safety': 1.1}),
            ('safety nan', {'depletion_safety': np.nan}),
            ('safety inf', {'depletion_safety': np.inf}),
            ('substeps zero', {'max_substeps': 0}),
            ('substeps negative', {'max_substeps': -1}),
            ('cache flag integer', {'cache_enabled': 1}),
        )
        for backend, kernel_type in available_kernel_types():
            _, kernel = make_kernel(kernel_type)
            for case_name, overrides in cases:
                with self.subTest(backend=backend, case=case_name):
                    self.assert_rejected_without_writes(
                        kernel, make_batch(), overrides=overrides
                    )

    def test_constructor_capacity_and_revision_limits(self):
        for backend, kernel_type in available_kernel_types():
            for capacity in (0, MAX_CACHE_CAPACITY + 1, True):
                with self.subTest(backend=backend, capacity=capacity):
                    arguments = kernel_arguments()[1]
                    arguments['cache_capacity'] = capacity
                    with self.assertRaises((TypeError, ValueError, OverflowError)):
                        kernel_type(**arguments)

            model, kernel = make_kernel(kernel_type)
            for revision in (-1, int(np.iinfo(np.uint64).max) + 1):
                with self.subTest(backend=backend, revision=revision):
                    with self.assertRaises((TypeError, ValueError, OverflowError)):
                        configure_kernel(
                            kernel,
                            model,
                            model_revision=revision,
                            prior_revision=1,
                        )


class TestMetabolicTransactions(unittest.TestCase):
    def test_max_substeps_rolls_back_arrays_and_cache(self):
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                _, kernel = make_kernel(kernel_type)
                # A pre-existing entry proves rollback restores cache ownership rather
                # than merely clearing failed work; see DOCUMENTATION.md §10.3.
                kernel.solve_flux(0, np.full(3, 0.2, dtype=np.float64))
                batch = make_batch()
                snapshots = {
                    name: batch[name].tobytes() for name in MUTABLE_BATCH_ARRAYS
                }
                starting_stats = kernel.cache_stats()
                with self.assertRaisesRegex(RuntimeError, 'max_substeps'):
                    kernel.advance_batch(
                        **batch,
                        dt_seconds=1.0,
                        low_substrate_threshold_m=1.0e-4,
                        depletion_safety=0.5,
                        max_substeps=1,
                        cache_enabled=True,
                    )

                for name, expected in snapshots.items():
                    self.assertEqual(batch[name].tobytes(), expected, name)
                self.assertEqual(kernel.cache_stats(), starting_stats)

    def test_infeasible_after_a_successful_substep_rolls_back_everything(self):
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                model, arguments = kernel_arguments()
                arguments['lower_bounds'] = np.array(arguments['lower_bounds'], copy=True)
                # A 0.5 generic-uptake floor is feasible at 1 M but infeasible after
                # depletion to zero, forcing a post-mutation failure; see DOCUMENTATION.md §10.3.
                arguments['lower_bounds'][0] = 0.5

                probe = kernel_type(**arguments)
                configure_kernel(probe, model)
                initial = make_batch(agent_count=1, concentration_m=0.0)
                initial['concentrations_m'][0] = (1.0, 0.0, 0.0)
                initial['biomass_gdw'][0] = 1.0e-5
                successful = copy_batch(initial)
                probe.advance_batch(
                    **successful,
                    dt_seconds=0.05,
                    low_substrate_threshold_m=1.0e-4,
                    depletion_safety=0.9,
                    max_substeps=8,
                    cache_enabled=True,
                )
                self.assertLess(successful['concentrations_m'][0, 0], 1.0)
                self.assertGreater(successful['biomass_gdw'][0], initial['biomass_gdw'][0])

                failing = kernel_type(**arguments)
                configure_kernel(failing, model)
                batch = copy_batch(initial)
                snapshots = {name: batch[name].tobytes() for name in MUTABLE_BATCH_ARRAYS}
                starting_stats = failing.cache_stats()
                with self.assertRaisesRegex(RuntimeError, 'infeasible'):
                    failing.advance_batch(
                        **batch,
                        dt_seconds=0.2,
                        low_substrate_threshold_m=1.0e-4,
                        depletion_safety=0.9,
                        max_substeps=8,
                        cache_enabled=True,
                    )
                for name, expected in snapshots.items():
                    self.assertEqual(batch[name].tobytes(), expected, name)
                self.assertEqual(failing.cache_stats(), starting_stats)


class TestMetabolicBufferOwnership(unittest.TestCase):
    def test_grow_shrink_regrow_batches_and_temporary_results_are_collectible(self):
        population_pattern = (2, 128, 5, 256, 3, 192, 7)
        for backend, kernel_type in available_kernel_types():
            with self.subTest(backend=backend):
                _, kernel = make_kernel(kernel_type, cache_capacity=16)
                references = []
                for agent_count in population_pattern:
                    batch = make_batch(agent_count=agent_count, concentration_m=0.2)
                    batch['environment_volume_l'][:] = 1.0e-2
                    batch['biomass_gdw'][:] = 1.0e-9
                    flux, _ = kernel.solve_flux(
                        0, np.full(3, 0.2, dtype=np.float64), cache_enabled=True
                    )
                    references.append(weakref.ref(flux))
                    references.extend(weakref.ref(batch[name]) for name in MUTABLE_BATCH_ARRAYS)
                    kernel.advance_batch(
                        **batch,
                        dt_seconds=0.001,
                        low_substrate_threshold_m=1.0e-4,
                        depletion_safety=0.9,
                        max_substeps=4,
                        cache_enabled=True,
                    )
                    del flux
                    del batch

                gc.collect()
                retained = sum(reference() is not None for reference in references)
                self.assertEqual(retained, 0, msg=f'{backend} retained {retained} temporary arrays')

    def test_python_retained_allocation_is_bounded_after_warmup(self):
        _, kernel = make_kernel(PythonMetabolicKernel, cache_capacity=16)
        population_pattern = (4, 64, 12, 128, 7)

        def run_iteration(iteration):
            agent_count = population_pattern[iteration % len(population_pattern)]
            # A large pool and short tick isolate allocation behavior from adaptive
            # depletion depth; see DOCUMENTATION.md §10.3.
            batch = make_batch(agent_count=agent_count, concentration_m=0.15)
            batch['environment_volume_l'][:] = 1.0e-2
            batch['biomass_gdw'][:] = 1.0e-9
            kernel.advance_batch(
                **batch,
                dt_seconds=0.001,
                low_substrate_threshold_m=1.0e-4,
                depletion_safety=0.9,
                max_substeps=4,
                cache_enabled=True,
            )

        for iteration in range(20):
            run_iteration(iteration)
        gc.collect()
        tracemalloc.start()
        try:
            for iteration in range(20):
                run_iteration(iteration)
            gc.collect()
            baseline_bytes = tracemalloc.get_traced_memory()[0]
            samples = []
            for iteration in range(150):
                run_iteration(iteration)
                if (iteration + 1) % 25 == 0:
                    gc.collect()
                    samples.append(tracemalloc.get_traced_memory()[0])
            retained_growth = samples[-1] - baseline_bytes
            self.assertLessEqual(
                retained_growth,
                1_048_576,
                msg=f'Python metabolic retention grew by {retained_growth} bytes',
            )
        finally:
            tracemalloc.stop()


@unittest.skipUnless(NATIVE_METABOLIC_AVAILABLE, 'native MetabolicKernel is not built')
class TestNativeMetabolicStress(unittest.TestCase):
    """Valgrind target; Python reachability checks do not prove native leak freedom."""

    def test_repeated_cache_churn_and_eviction(self):
        iteration_count = int(os.environ.get(
            'BACSANDBOX_NATIVE_METABOLIC_STRESS_ITERATIONS',
            str(NATIVE_METABOLIC_STRESS_DEFAULT_ITERATIONS),
        ))
        self.assertGreater(iteration_count, 0)
        self.assertLessEqual(iteration_count, 100_000)
        cache_capacity = 8
        _, kernel = make_kernel(
            bacteria_c_kernel.MetabolicKernel,
            cache_capacity=cache_capacity,
            phenotype_count=2,
        )
        solve_concentrations = np.empty(3, dtype=np.float64)
        batch = make_batch(agent_count=12, environment_count=3, concentration_m=0.1)
        batch['phenotype_id'][:] = np.arange(12, dtype=np.int32) % 2
        batch['environment_volume_l'][:] = 1.0e-3
        addresses = {name: value.ctypes.data for name, value in batch.items()}
        start = time.perf_counter()

        for iteration in range(iteration_count):
            # Four-bin strides force churn beyond capacity while remaining unsaturated;
            # see DOCUMENTATION.md §10.3.
            base = 0.05 + 4.0 * CACHE_TOLERANCE_MOLAR * (iteration % 32)
            solve_concentrations[:] = (
                base,
                base + 2.0 * CACHE_TOLERANCE_MOLAR,
                base + 4.0 * CACHE_TOLERANCE_MOLAR,
            )
            flux, _ = kernel.solve_flux(
                iteration % 2, solve_concentrations, cache_enabled=True
            )
            self.assertTrue(np.all(np.isfinite(flux)))
            self.assertTrue(np.all(flux >= 0.0))

            for environment in range(3):
                batch['concentrations_m'][environment] = (
                    solve_concentrations
                    + environment * 8.0 * CACHE_TOLERANCE_MOLAR
                )
            batch['biomass_gdw'].fill(1.0e-9)
            batch['growth_rate_s'].fill(-1.0)
            batch['biomass_delta_gdw'].fill(-1.0)
            batch['substrate_consumed_mol'].fill(-1.0)
            batch['acetate_secreted_mol'].fill(-1.0)
            summary = kernel.advance_batch(
                **batch,
                dt_seconds=0.001,
                low_substrate_threshold_m=1.0e-4,
                depletion_safety=0.9,
                max_substeps=4,
                cache_enabled=True,
            )
            self.assertGreaterEqual(summary[0], 1)
            self.assertTrue(np.all(np.isfinite(batch['concentrations_m'])))
            self.assertTrue(np.all(batch['concentrations_m'] >= 0.0))
            self.assertTrue(np.all(np.isfinite(batch['biomass_gdw'])))
            self.assertLessEqual(kernel.cache_stats()[0], cache_capacity)

        elapsed = time.perf_counter() - start
        for name, address in addresses.items():
            self.assertEqual(batch[name].ctypes.data, address, msg=name)
        self.assertGreater(kernel.cache_stats()[3], 0)
        print(
            f'Native Metabolic Stress: {iteration_count} iterations in {elapsed:.6f}s; '
            f'cache={kernel.cache_stats()}'
        )


def configure_scale_phenotypes(kernel, model):
    # Three synthetic parameter rows span neutral and inhibited uptake without
    # representing measured phenotypes; see DOCUMENTATION.md §10.3.
    kernel.configure_phenotypes(
        np.zeros(3, dtype=np.int32),
        np.asarray(((2.0, 1.8, 1.6), (1.5, 2.2, 1.9), (2.4, 1.7, 2.1)), dtype=np.float64),
        np.asarray(((1.0, 0.8, 1.2), (0.7, 1.1, 0.9), (1.3, 0.6, 1.0)), dtype=np.float64),
        np.asarray(((1.0, 1.0, 1.0), (2.0, 3.0, 4.0), (3.0, 2.0, 1.0)), dtype=np.float64),
        np.asarray(((1.0, 1.0, 1.0), (3.0, 2.0, 1.0), (2.0, 3.0, 4.0)), dtype=np.float64),
        np.asarray(((0.0, 0.0, 0.0), (0.4, 0.3, 0.2), (0.2, 0.4, 0.3)), dtype=np.float64),
        model.model_revision,
        19,
    )


def scale_batch():
    agent_count = 1000
    environment_count = 8
    agent_indices = np.arange(agent_count, dtype=np.int32)
    environment_indices = np.arange(environment_count, dtype=np.float64)[:, None]
    # Environment spacing exceeds cache tolerance, while 1 mL pools and nanogram
    # biomass keep each 0.05 s tick bounded; see DOCUMENTATION.md §10.3.
    concentrations = np.asarray((0.0020, 0.0025, 0.0030), dtype=np.float64) + (
        environment_indices * np.asarray((2.0e-4, 2.5e-4, 3.0e-4), dtype=np.float64)
    )
    return {
        'phenotype_id': agent_indices % 3,
        'environment_id': (agent_indices // 3) % environment_count,
        'concentrations_m': np.ascontiguousarray(concentrations),
        'environment_volume_l': np.full(environment_count, 1.0e-3, dtype=np.float64),
        'biomass_gdw': (1.0e-9 * (1.0 + 0.1 * (agent_indices % 5))).astype(np.float64),
        'growth_rate_s': np.zeros(agent_count, dtype=np.float64),
        'biomass_delta_gdw': np.zeros(agent_count, dtype=np.float64),
        'substrate_consumed_mol': np.zeros((agent_count, 3), dtype=np.float64),
        'acetate_secreted_mol': np.zeros(agent_count, dtype=np.float64),
    }


@unittest.skipUnless(RUN_METABOLIC_SCALE, 'set BACSANDBOX_RUN_METABOLIC_SCALE=1 to run')
@unittest.skipUnless(NATIVE_METABOLIC_AVAILABLE, 'native MetabolicKernel is not built')
class TestMetabolicScale(unittest.TestCase):
    def test_thousand_agents_repeated_adaptive_ticks(self):
        tick_count = 12
        model, python_kernel = make_kernel(PythonMetabolicKernel, cache_capacity=128)
        _, native_kernel = make_kernel(
            bacteria_c_kernel.MetabolicKernel, cache_capacity=128
        )
        configure_scale_phenotypes(python_kernel, model)
        configure_scale_phenotypes(native_kernel, model)
        python_batch = scale_batch()
        native_batch = copy_batch(python_batch)
        python_duration = 0.0
        native_duration = 0.0
        python_hits = 0
        native_hits = 0
        python_solves = 0
        native_solves = 0

        for tick in range(tick_count):
            start = time.perf_counter()
            python_summary = python_kernel.advance_batch(
                **python_batch,
                dt_seconds=0.05,
                low_substrate_threshold_m=1.0e-4,
                depletion_safety=0.9,
                max_substeps=8,
                cache_enabled=True,
            )
            python_duration += time.perf_counter() - start
            start = time.perf_counter()
            native_summary = native_kernel.advance_batch(
                **native_batch,
                dt_seconds=0.05,
                low_substrate_threshold_m=1.0e-4,
                depletion_safety=0.9,
                max_substeps=8,
                cache_enabled=True,
            )
            native_duration += time.perf_counter() - start

            self.assertEqual(native_summary[:5], python_summary[:5], msg=f'tick={tick}')
            self.assertAlmostEqual(native_summary[5], python_summary[5], delta=1.0e-12)
            for name in MUTABLE_BATCH_ARRAYS:
                np.testing.assert_allclose(
                    native_batch[name],
                    python_batch[name],
                    rtol=0.0,
                    atol=1.0e-5,
                    err_msg=f'{name} tick={tick}',
                )
            python_solves += python_summary[1]
            native_solves += native_summary[1]
            python_hits += python_summary[2]
            native_hits += native_summary[2]

        self.assertGreater(python_solves, 0)
        self.assertGreater(python_hits, 0)
        self.assertEqual(native_solves, python_solves)
        self.assertEqual(native_hits, python_hits)
        print(
            'Metabolic Scale Gate: '
            f'agents=1000; ticks={tick_count}; '
            f'Python={python_duration:.6f}s; Native={native_duration:.6f}s; '
            f'Python cache hits={python_hits}, solves={python_solves}; '
            f'Native cache hits={native_hits}, solves={native_solves}; '
            'fixture=synthetic reduced model, not a biological prediction/full Simulation'
        )


@unittest.skipUnless(RUN_FBA_BENCHMARK, 'set BACSANDBOX_RUN_FBA_BENCHMARK=1 to run')
@unittest.skipUnless(NATIVE_METABOLIC_AVAILABLE, 'native MetabolicKernel is not built')
class TestMetabolicCacheBenchmark(unittest.TestCase):
    def test_native_cache_on_off_median(self):
        state_count = 64
        repeat_count = 32
        sample_count = 5
        _, kernel = make_kernel(
            bacteria_c_kernel.MetabolicKernel, cache_capacity=128
        )
        # Four-bin spacing prevents tolerance aliasing; repeating the precomputed states
        # isolates native cache lookup from LP work; see DOCUMENTATION.md §10.3.
        states = tuple(
            np.asarray((
                0.02 + 4.0 * CACHE_TOLERANCE_MOLAR * index,
                0.03 + 8.0 * CACHE_TOLERANCE_MOLAR * (index % 16),
                0.04 + 12.0 * CACHE_TOLERANCE_MOLAR * (index % 8),
            ), dtype=np.float64)
            for index in range(state_count)
        )
        workload = states * repeat_count

        def run_sample(cache_enabled):
            kernel.clear_cache()
            output = np.empty((len(workload), 8), dtype=np.float64)
            start = time.perf_counter()
            for index, concentrations in enumerate(workload):
                output[index] = kernel.solve_flux(
                    0, concentrations, cache_enabled=cache_enabled
                )[0]
            return time.perf_counter() - start, output

        run_sample(False)
        run_sample(True)
        uncached_durations = []
        cached_durations = []
        uncached_output = None
        cached_output = None
        for _ in range(sample_count):
            duration, uncached_output = run_sample(False)
            uncached_durations.append(duration)
            duration, cached_output = run_sample(True)
            cached_durations.append(duration)

        np.testing.assert_allclose(cached_output, uncached_output, rtol=0.0, atol=1.0e-10)
        self.assertEqual(len(uncached_durations), sample_count)
        self.assertEqual(len(cached_durations), sample_count)
        self.assertEqual(
            kernel.cache_stats(),
            (state_count, len(workload) - state_count, state_count, 0, state_count),
        )
        uncached_median = statistics.median(uncached_durations)
        cached_median = statistics.median(cached_durations)
        self.assertGreater(uncached_median, 0.0)
        self.assertGreater(cached_median, 0.0)
        speedup = uncached_median / cached_median
        print(
            f'FBA Caching Speedup: {speedup:.3f}x '
            f'(Uncached median: {uncached_median:.6f}s; Cached median: {cached_median:.6f}s)'
        )


if __name__ == '__main__':
    unittest.main()