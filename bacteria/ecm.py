import importlib
import os
import sys

if __package__ is None and not hasattr(sys, 'frozen'):
    path = os.path.realpath(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(path)))
    __package__ = 'bacteria'
    
from .constants import (
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
)


_BACKENDS = frozenset(('auto', 'python', 'cpp'))
_NATIVE_MODULE_NAME = 'bacteria_c_kernel'
# Fixed ordering keeps each sparse layer aligned with its diffusion and decay parameters; see DOCUMENTATION.md §4.9.
_LAYER_ORDER = ('qs', 'bac', 'acetate', 'phage')
_LAYER_DIFFUSION = (
    ECM_QS_DIFFUSION_UM2_S,
    ECM_BAC_DIFFUSION_UM2_S,
    ECM_ACETATE_DIFFUSION_UM2_S,
    ECM_PHAGE_DIFFUSION_UM2_S,
)
_LAYER_DECAY = (
    ECM_QS_DECAY_S,
    ECM_BAC_DECAY_S,
    ECM_ACETATE_DECAY_S,
    ECM_PHAGE_DECAY_S,
)


class _NativeModuleUnavailable(ImportError):
    pass


class ECM:
    """
    Extracellular Matrix managing diffusion and degradation of environmental molecules.
    Uses a coarse spatial hash (dictionary grid) for fast lookups and low overhead.
    """
    def __init__(self, cell_size=ECM_CELL_SIZE_UM, backend='auto'):
        if backend not in _BACKENDS:
            raise ValueError("backend must be one of 'auto', 'python', or 'cpp'")

        self.cell_size = cell_size
        self._backend = backend
        self._backend_used = None
        self._native_import_attempted = False
        self._native_module = None
        self._native_import_error = None
        self.qs_inducers = {}
        self.bacteriocins = {}
        # Acetate cross-feeding metabolite (M2): overflow carbon source diffusing producers -> consumers.
        self.acetate = {}
        # Free bacteriophage virion field (M5): diffusible lytic-phage concentration driving SIR epidemics.
        self.phage = {}
        # Single source of truth mapping layer name -> backing grid, so add/sample never drift apart.
        self._layers = {
            'qs': self.qs_inducers,
            'bac': self.bacteriocins,
            'acetate': self.acetate,
            'phage': self.phage,
        }

    @property
    def backend_used(self):
        """Backend used by the most recent non-empty successful tick, if any."""
        return self._backend_used

    @property
    def native_available(self):
        """Whether the native module loaded, or None before any native load attempt."""
        if not self._native_import_attempted:
            return None
        return self._native_module is not None

    def add(self, x, y, amount, layer='qs'):
        """
        Add molecules to a specific layer at (x, y).
        """
        cx = int(x // self.cell_size)
        cy = int(y // self.cell_size)
        target = self._layers.get(layer, self.qs_inducers)
        target[(cx, cy)] = target.get((cx, cy), 0.0) + amount

    def sample(self, x, y, layer='qs'):
        """
        Sample the local concentration of a layer at (x, y).
        """
        cx = int(x // self.cell_size)
        cy = int(y // self.cell_size)
        target = self._layers.get(layer, self.qs_inducers)
        return target.get((cx, cy), 0.0)

    def tick(self, dt):
        if not any(self._layers[layer_name] for layer_name in _LAYER_ORDER):
            return

        if self._backend == 'python':
            self._tick_python(dt)
            self._backend_used = 'python'
            return

        try:
            native_module = self._load_native_backend()
        except _NativeModuleUnavailable:
            if self._backend == 'cpp':
                raise
            self._tick_python(dt)
            self._backend_used = 'python'
            return

        self._tick_cpp(native_module, dt)
        self._backend_used = 'cpp'

    def _tick_python(self, dt):
        self._process_layer(self.qs_inducers, ECM_QS_DIFFUSION_UM2_S, ECM_QS_DECAY_S, dt)
        self._process_layer(self.bacteriocins, ECM_BAC_DIFFUSION_UM2_S, ECM_BAC_DECAY_S, dt)
        self._process_layer(self.acetate, ECM_ACETATE_DIFFUSION_UM2_S, ECM_ACETATE_DECAY_S, dt)
        self._process_layer(self.phage, ECM_PHAGE_DIFFUSION_UM2_S, ECM_PHAGE_DECAY_S, dt)

    def _load_native_backend(self):
        if self._native_import_attempted:
            if self._native_import_error is not None:
                raise self._native_import_error
            return self._native_module

        self._native_import_attempted = True
        try:
            self._native_module = importlib.import_module(_NATIVE_MODULE_NAME)
        except ModuleNotFoundError as exc:
            if exc.name != _NATIVE_MODULE_NAME:
                self._native_import_error = exc
                raise
            error = _NativeModuleUnavailable(
                "The 'bacteria_c_kernel' extension is unavailable; build it or use backend='python'."
            )
            self._native_import_error = error
            raise error from exc
        except Exception as exc:
            self._native_import_error = exc
            raise

        return self._native_module

    def _tick_cpp(self, native_module, dt):
        import numpy as np

        total_size = sum(len(self._layers[layer_name]) for layer_name in _LAYER_ORDER)
        coordinates = np.empty((total_size, 2), dtype=np.int64)
        values = np.empty(total_size, dtype=np.float64)
        layer_offsets = np.empty(len(_LAYER_ORDER) + 1, dtype=np.int64)
        layer_offsets[0] = 0
        start = 0
        for layer_index, layer_name in enumerate(_LAYER_ORDER):
            grid = self._layers[layer_name]
            stop = start + len(grid)
            if stop > start:
                coordinates[start:stop] = tuple(grid)
                values[start:stop] = tuple(grid.values())
            layer_offsets[layer_index + 1] = stop
            start = stop

        diffusion = np.ascontiguousarray(_LAYER_DIFFUSION, dtype=np.float64)
        decay = np.ascontiguousarray(_LAYER_DECAY, dtype=np.float64)

        result = native_module.process_layers(
            coordinates,
            values,
            layer_offsets,
            diffusion,
            decay,
            dt,
            self.cell_size,
            ECM_ALPHA_MAX,
            ECM_MIN_CONCENTRATION,
        )
        rebuilt_layers = self._validate_native_result(result, np)

        for layer_name, rebuilt_layer in zip(_LAYER_ORDER, rebuilt_layers):
            grid = self._layers[layer_name]
            grid.clear()
            grid.update(rebuilt_layer)

    @staticmethod
    def _validate_native_result(result, np):
        if not isinstance(result, tuple) or len(result) != 3:
            raise TypeError('bacteria_c_kernel.process_layers must return a three-array tuple')

        coordinates, values, offsets = result
        arrays = (
            ('coordinates', coordinates, np.dtype(np.int64)),
            ('values', values, np.dtype(np.float64)),
            ('offsets', offsets, np.dtype(np.int64)),
        )
        for name, array, expected_dtype in arrays:
            if not isinstance(array, np.ndarray):
                raise TypeError(f'native {name} result must be a NumPy array')
            if array.dtype != expected_dtype:
                raise TypeError(f'native {name} result has an invalid dtype')
            if not array.flags.c_contiguous:
                raise ValueError(f'native {name} result must be C-contiguous')

        if coordinates.ndim != 2 or coordinates.shape[1] != 2:
            raise ValueError('native coordinates result must have shape (N, 2)')
        if values.ndim != 1 or values.shape[0] != coordinates.shape[0]:
            raise ValueError('native values result must have shape (N,)')
        if offsets.ndim != 1 or offsets.shape[0] != len(_LAYER_ORDER) + 1:
            raise ValueError('native offsets result must have shape (5,)')
        if offsets[0] != 0 or offsets[-1] != values.shape[0] or np.any(offsets[1:] < offsets[:-1]):
            raise ValueError('native offsets result must be monotonic from zero through N')
        if not np.all(np.isfinite(values)) or np.any(values < ECM_MIN_CONCENTRATION):
            raise ValueError('native values result must contain finite, retained concentrations')

        rebuilt_layers = []
        for layer_index in range(len(_LAYER_ORDER)):
            start = int(offsets[layer_index])
            stop = int(offsets[layer_index + 1])
            coordinate_rows = map(tuple, coordinates[start:stop].tolist())
            rebuilt_layer = dict(zip(coordinate_rows, values[start:stop].tolist()))
            if len(rebuilt_layer) != stop - start:
                raise ValueError('native coordinates result contains a duplicate within a layer')
            rebuilt_layers.append(rebuilt_layer)

        return rebuilt_layers

    def _process_layer_python(self, grid, diffusion_um2_s, decay_s, dt):
        if not grid:
            return

        # This explicit Fickian stencil caps alpha = D * dt / dx^2 below 0.25 to preserve nonnegative mass; see DOCUMENTATION.md §4.9.
        alpha = min(ECM_ALPHA_MAX, max(0.0, diffusion_um2_s * dt / (self.cell_size ** 2)))
        decay_multiplier = max(0.0, min(1.0, 1.0 - decay_s * dt))
        next_grid = {}

        for k, val in list(grid.items()):
            if val < ECM_MIN_CONCENTRATION:
                continue

            decayed = val * decay_multiplier
            retained = decayed * (1.0 - 4.0 * alpha)
            neighbor_share = decayed * alpha

            if retained >= ECM_MIN_CONCENTRATION:
                next_grid[k] = next_grid.get(k, 0.0) + retained

            cx, cy = k
            if neighbor_share >= ECM_MIN_CONCENTRATION:
                for nx, ny in ((cx - 1, cy), (cx + 1, cy), (cx, cy - 1), (cx, cy + 1)):
                    next_grid[(nx, ny)] = next_grid.get((nx, ny), 0.0) + neighbor_share

        grid.clear()
        grid.update({k: v for k, v in next_grid.items() if v >= ECM_MIN_CONCENTRATION})

    def _process_layer(self, grid, diffusion_um2_s, decay_s, dt):
        return self._process_layer_python(grid, diffusion_um2_s, decay_s, dt)
