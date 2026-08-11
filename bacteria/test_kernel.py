import gc
import os
import random
import statistics
import time
import unittest
import weakref
from unittest.mock import Mock, patch

import numpy as np

from bacteria.constants import (
    ECM_ACETATE_DECAY_S,
    ECM_ACETATE_DIFFUSION_UM2_S,
    ECM_ALPHA_MAX,
    ECM_BAC_DECAY_S,
    ECM_BAC_DIFFUSION_UM2_S,
    ECM_CELL_SIZE_UM,
    ECM_MIN_CONCENTRATION,
    ECM_PHAGE_DECAY_S,
    ECM_PHAGE_DIFFUSION_UM2_S,
    ECM_QS_DECAY_S,
    ECM_QS_DIFFUSION_UM2_S,
    FPS,
)
from bacteria.ecm import ECM


try:
    import bacteria_c_kernel
except ImportError:
    bacteria_c_kernel = None


NATIVE_AVAILABLE = bacteria_c_kernel is not None
# A test step is one normal real-time frame expressed in seconds; see DOCUMENTATION.md §4.9.
TEST_DT_SECONDS = 1.0 / FPS
LAYER_NAMES = ('qs', 'bac', 'acetate', 'phage')
LAYER_DIFFUSION = (
    ECM_QS_DIFFUSION_UM2_S,
    ECM_BAC_DIFFUSION_UM2_S,
    ECM_ACETATE_DIFFUSION_UM2_S,
    ECM_PHAGE_DIFFUSION_UM2_S,
)
LAYER_DECAY = (
    ECM_QS_DECAY_S,
    ECM_BAC_DECAY_S,
    ECM_ACETATE_DECAY_S,
    ECM_PHAGE_DECAY_S,
)
# This concentration remains well above the sparse-grid floor during long parity runs; see DOCUMENTATION.md §4.9.
ACTIVE_CONCENTRATION = 25.0
# These concentrations straddle the prune floor without equality ambiguity; see DOCUMENTATION.md §4.9.
BELOW_FLOOR_CONCENTRATION = ECM_MIN_CONCENTRATION * 0.75
ABOVE_FLOOR_CONCENTRATION = ECM_MIN_CONCENTRATION * 1.25
# The fixed fixture seed and 100-step horizon make long-run sparse diffusion reproducible; see DOCUMENTATION.md §4.9.
PARITY_SEED = 20260808
PARITY_TICK_COUNT = 100
PARITY_CELLS_PER_LAYER = 48
# Native and Python concentrations may differ by floating-point rounding but not model-scale error; see DOCUMENTATION.md §4.9.
PARITY_ABS_TOLERANCE = 1.0e-5
# Repeated short batches exercise ownership without making the ordinary suite memory-heavy; see DOCUMENTATION.md §4.9.
NATIVE_STRESS_DEFAULT_ITERATIONS = 250
# The opt-in workload is large enough to time sparse-grid processing while remaining bounded; see DOCUMENTATION.md §4.9.
BENCHMARK_SEED = 90210
BENCHMARK_GRID_SIDE = 40
BENCHMARK_TICK_COUNT = 8
BENCHMARK_REPEAT_COUNT = 3
RUN_BENCHMARK = os.environ.get('BACSANDBOX_RUN_BENCHMARK') == '1'


def missing_native_module():
    return ModuleNotFoundError(
        "No module named 'bacteria_c_kernel'",
        name='bacteria_c_kernel',
    )


def parity_fixture():
    rng = random.Random(PARITY_SEED)
    island_centers = ((-48, -37), (0, 0), (53, -29), (-22, 61), (71, 58))
    fixture = {}
    for layer_index, layer_name in enumerate(LAYER_NAMES):
        rows = [
            ((-240 - layer_index, 180 + layer_index), BELOW_FLOOR_CONCENTRATION),
            ((240 + layer_index, -180 - layer_index), ABOVE_FLOOR_CONCENTRATION),
        ]
        occupied = {coordinate for coordinate, _ in rows}

        # Adjacent crosses force overlapping neighbor contributions inside separated signed-coordinate islands; see DOCUMENTATION.md §4.9.
        for center_x, center_y in island_centers:
            for offset_x, offset_y in ((0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)):
                coordinate = (center_x + offset_x, center_y + offset_y)
                occupied.add(coordinate)
                concentration = ACTIVE_CONCENTRATION * (1.0 + 0.2 * layer_index) + rng.uniform(0.0, 10.0)
                rows.append((coordinate, concentration))

        # Seeded local scatter supplies many active cells while preserving finite island bounds; see DOCUMENTATION.md §4.9.
        while len(rows) < PARITY_CELLS_PER_LAYER:
            center_x, center_y = island_centers[rng.randrange(len(island_centers))]
            coordinate = (center_x + rng.randint(-4, 4), center_y + rng.randint(-4, 4))
            if coordinate in occupied:
                continue
            occupied.add(coordinate)
            concentration = ACTIVE_CONCENTRATION + rng.uniform(0.0, 20.0)
            rows.append((coordinate, concentration))

        fixture[layer_name] = rows
    return fixture


def benchmark_fixture():
    rng = random.Random(BENCHMARK_SEED)
    half_side = BENCHMARK_GRID_SIDE // 2
    fixture = {}
    for layer_index, layer_name in enumerate(LAYER_NAMES):
        rows = []
        for row in range(BENCHMARK_GRID_SIDE):
            for column in range(BENCHMARK_GRID_SIDE):
                coordinate = (column - half_side, row - half_side)
                concentration = ACTIVE_CONCENTRATION * (1.0 + 0.1 * layer_index) + rng.random()
                rows.append((coordinate, concentration))
        fixture[layer_name] = rows
    return fixture


def make_ecm(backend, fixture):
    ecm = ECM(cell_size=ECM_CELL_SIZE_UM, backend=backend)
    for layer_name in LAYER_NAMES:
        for coordinate, concentration in fixture[layer_name]:
            ecm._layers[layer_name][coordinate] = concentration
    return ecm


def snapshot_layers(ecm):
    return {layer_name: dict(ecm._layers[layer_name]) for layer_name in LAYER_NAMES}


def assert_layer_maps_close(test_case, expected, actual, context):
    for layer_name in LAYER_NAMES:
        expected_layer = expected[layer_name]
        actual_layer = actual[layer_name]
        test_case.assertEqual(
            set(expected_layer),
            set(actual_layer),
            msg=f'{context}, layer={layer_name}: key sets differ',
        )
        for coordinate, expected_value in expected_layer.items():
            difference = abs(expected_value - actual_layer[coordinate])
            test_case.assertLessEqual(
                difference,
                PARITY_ABS_TOLERANCE,
                msg=(
                    f'{context}, layer={layer_name}, cell={coordinate}: '
                    f'{expected_value!r} != {actual_layer[coordinate]!r}'
                ),
            )


def valid_native_arguments(row_count=1):
    coordinates = np.column_stack((
        np.arange(row_count, dtype=np.int64),
        -np.arange(row_count, dtype=np.int64),
    ))
    return {
        'coordinates': np.ascontiguousarray(coordinates, dtype=np.int64).reshape((-1, 2)),
        'values': np.full(row_count, ACTIVE_CONCENTRATION, dtype=np.float64),
        'offsets': np.asarray((0, row_count, row_count, row_count, row_count), dtype=np.int64),
        'diffusion': np.asarray(LAYER_DIFFUSION, dtype=np.float64),
        'decay': np.asarray(LAYER_DECAY, dtype=np.float64),
        'dt': TEST_DT_SECONDS,
        'cell_size': float(ECM_CELL_SIZE_UM),
        'alpha_max': ECM_ALPHA_MAX,
        'min_concentration': ECM_MIN_CONCENTRATION,
    }


def native_result_for_layers(
    layers,
    diffusion=LAYER_DIFFUSION,
    decay=LAYER_DECAY,
    dt=TEST_DT_SECONDS,
    alpha_max=ECM_ALPHA_MAX,
):
    coordinates = []
    values = []
    offsets = [0]
    for layer_name in LAYER_NAMES:
        for coordinate, concentration in layers[layer_name].items():
            coordinates.append(coordinate)
            values.append(concentration)
        offsets.append(len(values))

    return bacteria_c_kernel.process_layers(
        np.asarray(coordinates, dtype=np.int64).reshape((-1, 2)),
        np.asarray(values, dtype=np.float64),
        np.asarray(offsets, dtype=np.int64),
        np.asarray(diffusion, dtype=np.float64),
        np.asarray(decay, dtype=np.float64),
        dt,
        float(ECM_CELL_SIZE_UM),
        alpha_max,
        ECM_MIN_CONCENTRATION,
    )


def validated_native_layers(result):
    rebuilt = ECM._validate_native_result(result, np)
    return dict(zip(LAYER_NAMES, rebuilt))


def python_result_for_layers(layers, dt=TEST_DT_SECONDS):
    ecm = ECM(cell_size=ECM_CELL_SIZE_UM, backend='python')
    for layer_name in LAYER_NAMES:
        ecm._layers[layer_name].update(layers[layer_name])
    ecm.tick(dt)
    return snapshot_layers(ecm)


def native_single_layer(grid, diffusion, decay, dt, alpha_max=ECM_ALPHA_MAX):
    layers = {layer_name: {} for layer_name in LAYER_NAMES}
    layers['qs'].update(grid)
    result = native_result_for_layers(
        layers,
        diffusion=(diffusion, 0.0, 0.0, 0.0),
        decay=(decay, 0.0, 0.0, 0.0),
        dt=dt,
        alpha_max=alpha_max,
    )
    return validated_native_layers(result)['qs']


def python_single_layer(grid, diffusion, decay, dt):
    reference = dict(grid)
    ecm = ECM(cell_size=ECM_CELL_SIZE_UM, backend='python')
    ecm._process_layer_python(reference, diffusion, decay, dt)
    return reference


class TestECMBackendDispatch(unittest.TestCase):
    def test_invalid_backend_rejected(self):
        with self.assertRaisesRegex(ValueError, "auto.*python.*cpp"):
            ECM(backend='invalid')

    def test_python_backend_never_imports_native_module(self):
        ecm = ECM(backend='python')
        ecm.add(0.0, 0.0, 10.0, layer='qs')

        with patch('bacteria.ecm.importlib.import_module') as import_module:
            ecm.tick(TEST_DT_SECONDS)

        import_module.assert_not_called()
        self.assertEqual(ecm.backend_used, 'python')
        self.assertIsNone(ecm.native_available)

    def test_auto_backend_skips_empty_ticks_then_caches_fallback(self):
        ecm = ECM(backend='auto')
        with patch(
            'bacteria.ecm.importlib.import_module',
            side_effect=missing_native_module(),
        ) as import_module:
            ecm.tick(TEST_DT_SECONDS)
            import_module.assert_not_called()

            ecm.add(0.0, 0.0, 10.0, layer='qs')
            ecm.tick(TEST_DT_SECONDS)
            ecm.tick(TEST_DT_SECONDS)

        import_module.assert_called_once_with('bacteria_c_kernel')
        self.assertEqual(ecm.backend_used, 'python')
        self.assertFalse(ecm.native_available)

    def test_cpp_backend_reports_unavailable_extension(self):
        ecm = ECM(backend='cpp')
        ecm.add(0.0, 0.0, 10.0, layer='qs')

        with patch(
            'bacteria.ecm.importlib.import_module',
            side_effect=missing_native_module(),
        ):
            with self.assertRaisesRegex(
                ImportError,
                "bacteria_c_kernel.*unavailable.*backend='python'",
            ):
                ecm.tick(TEST_DT_SECONDS)

    def test_auto_backend_propagates_native_processing_failure(self):
        native_module = Mock()
        native_module.process_layers.side_effect = RuntimeError('native processing failed')
        ecm = ECM(backend='auto')
        ecm.add(0.0, 0.0, 10.0, layer='qs')

        with patch('bacteria.ecm.importlib.import_module', return_value=native_module):
            with self.assertRaisesRegex(RuntimeError, 'native processing failed'):
                ecm.tick(TEST_DT_SECONDS)

        self.assertIsNone(ecm.backend_used)
        self.assertTrue(ecm.native_available)

    def test_missing_native_symbol_falls_back_only_in_auto_mode(self):
        incomplete_module = object()
        auto = ECM(backend='auto')
        auto.add(0.0, 0.0, ACTIVE_CONCENTRATION, layer='qs')
        with patch('bacteria.ecm.importlib.import_module', return_value=incomplete_module):
            auto.tick(TEST_DT_SECONDS)
        self.assertEqual(auto.backend_used, 'python')
        self.assertFalse(auto.native_available)

        explicit = ECM(backend='cpp')
        explicit.add(0.0, 0.0, ACTIVE_CONCENTRATION, layer='qs')
        with patch('bacteria.ecm.importlib.import_module', return_value=incomplete_module):
            with self.assertRaisesRegex(ImportError, 'lacks process_layers'):
                explicit.tick(TEST_DT_SECONDS)

    @unittest.skipUnless(NATIVE_AVAILABLE, 'bacteria_c_kernel is not built')
    def test_native_tick_preserves_layer_dictionary_identities(self):
        ecm = ECM(backend='cpp')
        original_layers = {layer_name: ecm._layers[layer_name] for layer_name in LAYER_NAMES}
        for layer_index, layer_name in enumerate(LAYER_NAMES):
            ecm._layers[layer_name][(layer_index, -layer_index)] = ACTIVE_CONCENTRATION

        ecm.tick(TEST_DT_SECONDS)

        self.assertEqual(ecm.backend_used, 'cpp')
        for layer_name in LAYER_NAMES:
            self.assertIs(ecm._layers[layer_name], original_layers[layer_name])


@unittest.skipUnless(NATIVE_AVAILABLE, 'bacteria_c_kernel is not built')
class TestKernelParity(unittest.TestCase):
    def test_seeded_four_layer_parity_after_every_tick(self):
        fixture = parity_fixture()
        python_ecm = make_ecm('python', fixture)
        native_ecm = make_ecm('cpp', fixture)

        for tick_index in range(1, PARITY_TICK_COUNT + 1):
            python_ecm.tick(TEST_DT_SECONDS)
            native_ecm.tick(TEST_DT_SECONDS)
            assert_layer_maps_close(
                self,
                snapshot_layers(python_ecm),
                snapshot_layers(native_ecm),
                context=f'tick={tick_index}',
            )

        for layer_name in LAYER_NAMES:
            self.assertTrue(python_ecm._layers[layer_name], msg=f'layer={layer_name} became inactive')


@unittest.skipUnless(NATIVE_AVAILABLE, 'bacteria_c_kernel is not built')
class TestNativeKernelEdges(unittest.TestCase):
    def assert_grid_close(self, expected, actual, context):
        expected_layers = {layer_name: {} for layer_name in LAYER_NAMES}
        actual_layers = {layer_name: {} for layer_name in LAYER_NAMES}
        expected_layers['qs'] = expected
        actual_layers['qs'] = actual
        assert_layer_maps_close(self, expected_layers, actual_layers, context)

    def assert_native_rejected(self, overrides, exception_type, message_pattern):
        arguments = valid_native_arguments()
        arguments.update(overrides)
        with self.assertRaisesRegex(exception_type, message_pattern):
            bacteria_c_kernel.process_layers(**arguments)

    def test_empty_layer_segments_match_python_reference(self):
        layers = {layer_name: {} for layer_name in LAYER_NAMES}
        layers['bac'] = {(-3, 4): ACTIVE_CONCENTRATION, (7, -2): ACTIVE_CONCENTRATION * 1.5}

        result = native_result_for_layers(layers)
        native_layers = validated_native_layers(result)
        python_layers = python_result_for_layers(layers)

        self.assertEqual(int(result[2][0]), int(result[2][1]))
        self.assertEqual(int(result[2][2]), int(result[2][3]))
        self.assertEqual(int(result[2][3]), int(result[2][4]))
        assert_layer_maps_close(self, python_layers, native_layers, context='empty segments')

        empty_result = native_result_for_layers({layer_name: {} for layer_name in LAYER_NAMES})
        self.assertEqual(empty_result[0].shape, (0, 2))
        self.assertEqual(empty_result[1].shape, (0,))
        np.testing.assert_array_equal(empty_result[2], np.zeros(5, dtype=np.int64))

    def test_alpha_max_clamp_matches_python_reference(self):
        grid = {(0, 0): ACTIVE_CONCENTRATION}
        # This coefficient makes raw alpha equal 2, so the documented stability ceiling must control the stencil; see DOCUMENTATION.md §4.9.
        diffusion = 2.0 * (ECM_CELL_SIZE_UM ** 2) / TEST_DT_SECONDS

        expected = python_single_layer(grid, diffusion, 0.0, TEST_DT_SECONDS)
        actual = native_single_layer(grid, diffusion, 0.0, TEST_DT_SECONDS)

        self.assert_grid_close(expected, actual, context='alpha max clamp')

    def test_decay_clamp_at_large_dt_matches_python_reference(self):
        grid = {(0, 0): ACTIVE_CONCENTRATION}
        # Two decay lifetimes make 1 - decay*dt negative, so the nonnegative decay clamp removes the source; see DOCUMENTATION.md §4.9.
        large_dt = 2.0 / ECM_QS_DECAY_S

        expected = python_single_layer(grid, ECM_QS_DIFFUSION_UM2_S, ECM_QS_DECAY_S, large_dt)
        actual = native_single_layer(grid, ECM_QS_DIFFUSION_UM2_S, ECM_QS_DECAY_S, large_dt)

        self.assert_grid_close(expected, actual, context='decay clamp')
        self.assertEqual(actual, {})

    def test_floor_skip_overlap_and_final_prune_match_python_reference(self):
        # A short unit-decay step prunes the near-floor source while adjacent active sources overlap at (1, 0); see DOCUMENTATION.md §4.9.
        dt = 0.02
        decay = 1.0
        grid = {
            (-30, -30): BELOW_FLOOR_CONCENTRATION,
            (30, 30): ECM_MIN_CONCENTRATION * 1.01,
            (0, 0): ACTIVE_CONCENTRATION,
            (2, 0): ACTIVE_CONCENTRATION,
        }

        expected = python_single_layer(grid, ECM_QS_DIFFUSION_UM2_S, decay, dt)
        actual = native_single_layer(grid, ECM_QS_DIFFUSION_UM2_S, decay, dt)

        self.assert_grid_close(expected, actual, context='floor and overlap')
        self.assertNotIn((-30, -30), actual)
        self.assertNotIn((30, 30), actual)
        self.assertGreater(actual[(1, 0)], ECM_MIN_CONCENTRATION)

    def test_signed_coordinates_are_unbounded_and_match_reference(self):
        grid = {
            (-5, -7): ACTIVE_CONCENTRATION,
            (6, -8): ACTIVE_CONCENTRATION * 1.5,
            (-11, 9): ACTIVE_CONCENTRATION * 2.0,
        }

        expected = python_single_layer(grid, ECM_QS_DIFFUSION_UM2_S, 0.0, TEST_DT_SECONDS)
        actual = native_single_layer(grid, ECM_QS_DIFFUSION_UM2_S, 0.0, TEST_DT_SECONDS)

        self.assert_grid_close(expected, actual, context='signed coordinates')
        self.assertIn((-6, -7), actual)
        self.assertIn((-11, 10), actual)

    def test_int64_neighbor_overflow_raises(self):
        int64_bounds = (np.iinfo(np.int64).min, np.iinfo(np.int64).max)
        for coordinate_x in int64_bounds:
            with self.subTest(coordinate_x=int(coordinate_x)):
                with self.assertRaisesRegex(OverflowError, 'int64 range'):
                    native_single_layer(
                        {(int(coordinate_x), 0): ACTIVE_CONCENTRATION},
                        ECM_QS_DIFFUSION_UM2_S,
                        0.0,
                        TEST_DT_SECONDS,
                    )

    def test_wrong_dtypes_rejected(self):
        cases = (
            ('coordinates', np.asarray(((0, 0),), dtype=np.int32)),
            ('values', np.asarray((ACTIVE_CONCENTRATION,), dtype=np.float32)),
            ('offsets', np.asarray((0, 1, 1, 1, 1), dtype=np.int32)),
            ('diffusion', np.asarray(LAYER_DIFFUSION, dtype=np.float32)),
            ('decay', np.asarray(LAYER_DECAY, dtype=np.float32)),
        )
        for argument_name, invalid_array in cases:
            with self.subTest(argument_name=argument_name):
                self.assert_native_rejected(
                    {argument_name: invalid_array},
                    TypeError,
                    f'{argument_name}.*invalid dtype',
                )

    def test_wrong_shapes_rejected(self):
        cases = (
            ('coordinates', np.asarray((0, 0), dtype=np.int64), r'coordinates.*\(N, 2\)'),
            ('coordinates', np.zeros((1, 3), dtype=np.int64), r'coordinates.*\(N, 2\)'),
            ('values', np.asarray(((ACTIVE_CONCENTRATION,),), dtype=np.float64), r'values.*\(N,\)'),
            ('values', np.asarray((ACTIVE_CONCENTRATION, ACTIVE_CONCENTRATION), dtype=np.float64), r'values.*\(N,\)'),
            ('offsets', np.asarray((0, 1, 1, 1), dtype=np.int64), r'offsets.*\(5,\)'),
            ('diffusion', np.ones(3, dtype=np.float64), r'diffusion.*\(4,\)'),
            ('decay', np.ones((2, 2), dtype=np.float64), r'decay.*\(4,\)'),
        )
        for argument_name, invalid_array, message_pattern in cases:
            with self.subTest(argument_name=argument_name, shape=invalid_array.shape):
                self.assert_native_rejected({argument_name: invalid_array}, ValueError, message_pattern)

    def test_noncontiguous_arrays_rejected(self):
        cases = (
            ('coordinates', np.arange(8, dtype=np.int64).reshape(2, 4)[:, ::2]),
            ('values', np.arange(4, dtype=np.float64)[::2]),
            ('offsets', np.asarray((0, 9, 2, 9, 2, 9, 2, 9, 2, 9), dtype=np.int64)[::2]),
            ('diffusion', np.arange(8, dtype=np.float64)[::2]),
            ('decay', np.arange(8, dtype=np.float64)[::2]),
        )
        for argument_name, invalid_array in cases:
            with self.subTest(argument_name=argument_name):
                self.assertFalse(invalid_array.flags.c_contiguous)
                arguments = valid_native_arguments(row_count=2)
                arguments[argument_name] = invalid_array
                with self.assertRaisesRegex(ValueError, f'{argument_name}.*C-contiguous'):
                    bacteria_c_kernel.process_layers(**arguments)

    def test_invalid_offsets_rejected(self):
        cases = (
            (np.asarray((1, 2, 2, 2, 2), dtype=np.int64), 'start at zero'),
            (np.asarray((0, 1, 1, 1, 1), dtype=np.int64), 'end at N'),
            (np.asarray((0, 2, 1, 2, 2), dtype=np.int64), 'monotonic'),
            (np.asarray((0, -1, 1, 2, 2), dtype=np.int64), 'nonnegative'),
        )
        for offsets, message_pattern in cases:
            with self.subTest(offsets=offsets.tolist()):
                arguments = valid_native_arguments(row_count=2)
                arguments['offsets'] = offsets
                with self.assertRaisesRegex(ValueError, message_pattern):
                    bacteria_c_kernel.process_layers(**arguments)

    def test_nonfinite_concentrations_rejected(self):
        for invalid_value in (float('nan'), float('inf'), float('-inf')):
            with self.subTest(value=invalid_value):
                self.assert_native_rejected(
                    {'values': np.asarray((invalid_value,), dtype=np.float64)},
                    ValueError,
                    'finite concentrations',
                )

    def test_nonfinite_and_negative_parameters_rejected(self):
        for argument_name in ('diffusion', 'decay'):
            for invalid_value in (-1.0, float('nan'), float('inf')):
                with self.subTest(argument_name=argument_name, value=invalid_value):
                    values = np.asarray(LAYER_DIFFUSION if argument_name == 'diffusion' else LAYER_DECAY)
                    values = values.astype(np.float64)
                    values[0] = invalid_value
                    self.assert_native_rejected(
                        {argument_name: values},
                        ValueError,
                        'finite and nonnegative',
                    )

        for argument_name in ('dt', 'alpha_max', 'min_concentration'):
            for invalid_value in (-1.0, float('nan'), float('inf')):
                with self.subTest(argument_name=argument_name, value=invalid_value):
                    self.assert_native_rejected(
                        {argument_name: invalid_value},
                        ValueError,
                        'finite and nonnegative',
                    )

        for invalid_value in (0.0, -1.0, float('nan'), float('inf')):
            with self.subTest(cell_size=invalid_value):
                self.assert_native_rejected(
                    {'cell_size': invalid_value},
                    ValueError,
                    'finite and positive',
                )


@unittest.skipUnless(NATIVE_AVAILABLE, 'bacteria_c_kernel is not built')
class TestNativeAllocationStress(unittest.TestCase):
    def test_repeated_results_are_collectible(self):
        iteration_count = int(os.environ.get(
            'BACSANDBOX_NATIVE_STRESS_ITERATIONS',
            str(NATIVE_STRESS_DEFAULT_ITERATIONS),
        ))
        references = []

        for iteration in range(iteration_count):
            # Changing row counts and signed coordinates exercises owned result buffers across sparse batches; see DOCUMENTATION.md §4.9.
            row_count = 1 + iteration % 7
            arguments = valid_native_arguments(row_count=row_count)
            arguments['coordinates'][:, 0] += iteration
            arguments['coordinates'][:, 1] -= iteration
            result = bacteria_c_kernel.process_layers(**arguments)
            references.extend(weakref.ref(array) for array in result)
            del result

        gc.collect()
        retained = sum(reference() is not None for reference in references)
        self.assertEqual(retained, 0, msg=f'{retained} native result arrays remained reachable')


@unittest.skipUnless(RUN_BENCHMARK, 'set BACSANDBOX_RUN_BENCHMARK=1 to run')
@unittest.skipUnless(NATIVE_AVAILABLE, 'bacteria_c_kernel is not built')
class TestKernelBenchmark(unittest.TestCase):
    def run_once(self, backend, fixture, tick_count):
        ecm = make_ecm(backend, fixture)
        start = time.perf_counter()
        for _ in range(tick_count):
            ecm.tick(TEST_DT_SECONDS)
        elapsed = time.perf_counter() - start
        return elapsed, snapshot_layers(ecm)

    def timed_repeats(self, backend, fixture):
        durations = []
        final_layers = None
        for _ in range(BENCHMARK_REPEAT_COUNT):
            duration, final_layers = self.run_once(backend, fixture, BENCHMARK_TICK_COUNT)
            durations.append(duration)
        return statistics.median(durations), final_layers

    def test_cpp_speedup(self):
        fixture = benchmark_fixture()
        self.run_once('python', fixture, 1)
        self.run_once('cpp', fixture, 1)

        python_seconds, python_layers = self.timed_repeats('python', fixture)
        cpp_seconds, cpp_layers = self.timed_repeats('cpp', fixture)

        assert_layer_maps_close(self, python_layers, cpp_layers, context='benchmark output')
        speedup = python_seconds / cpp_seconds
        print(
            f'C++ Speedup Factor: {speedup:.2f}x '
            f'(Python median: {python_seconds:.6f}s; C++ median: {cpp_seconds:.6f}s)'
        )


if __name__ == '__main__':
    unittest.main()