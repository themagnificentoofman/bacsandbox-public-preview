from dataclasses import dataclass
import importlib
import math
from numbers import Integral, Real
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
    SUBSTRATE_NAMES,
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


@dataclass(frozen=True)
class SubstrateLayerSpec:
    """Configuration for one continuous molar substrate field."""

    name: str
    diffusion_um2_s: float
    decay_rate_s_inv: float
    min_concentration_m: float
    provenance: str


@dataclass(frozen=True)
class SubstrateLayerSnapshot:
    """Immutable copy of substrate registration and all field values."""

    layer_specs: tuple
    layer_order: tuple
    layer_items: tuple


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
        # FBA substrate fields are deliberately separate from the fixed four-layer
        # signal ABI consumed by bacteria_c_kernel.process_layers.
        self._substrate_layer_specs = ()
        self._substrate_layers = {}

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

    @property
    def substrate_layer_names(self):
        """Return the registered molar substrate fields in canonical ABI order."""
        return tuple(spec.name for spec in self._substrate_layer_specs)

    def register_substrate_layers(self, specs):
        """Register exactly the three canonical continuous substrate fields."""
        specs = self._validate_substrate_specs(specs)

        if self._substrate_layer_specs:
            if specs != self._substrate_layer_specs:
                raise RuntimeError('substrate layers are already registered with different specs')
            return

        self._substrate_layer_specs = specs
        self._substrate_layers = {spec.name: {} for spec in specs}

    def world_to_cell(self, x_um, y_um):
        """Map finite world coordinates in micrometers to one sparse grid key."""
        for name, value in (('x_um', x_um), ('y_um', y_um)):
            if isinstance(value, bool) or not isinstance(value, Real):
                raise TypeError(f'{name} must be a real number')
            if not math.isfinite(float(value)):
                raise ValueError(f'{name} must be finite')
        return (
            int(math.floor(float(x_um) / self.cell_size)),
            int(math.floor(float(y_um) / self.cell_size)),
        )

    def deposit_substrates_m(self, coordinate, increment_m):
        """Add one canonical molar concentration vector at a grid coordinate."""
        import numpy as np

        coordinate = self._validate_substrate_coordinate(coordinate)
        values = self._validate_substrate_matrix(
            increment_m,
            expected_shape=(len(SUBSTRATE_NAMES),),
            name='increment_m',
        )
        updated = []
        for layer_index, spec in enumerate(self._substrate_layer_specs):
            current = self._substrate_layers[spec.name].get(coordinate, 0.0)
            value = current + float(values[layer_index])
            if not np.isfinite(value):
                raise OverflowError('deposited substrate concentration is not finite')
            if 0.0 < value < spec.min_concentration_m:
                raise ValueError(
                    f'deposited {spec.name} concentration is below its retention floor'
                )
            updated.append(value)
        for spec, value in zip(self._substrate_layer_specs, updated):
            if value >= spec.min_concentration_m:
                self._substrate_layers[spec.name][coordinate] = value

    def sample_substrate_batch_m(self, coordinates):
        """Return fresh C-contiguous ``(E, 3)`` concentrations for unique cells."""
        import numpy as np

        coordinates = self._validate_substrate_coordinates(coordinates)
        result = np.zeros((len(coordinates), len(SUBSTRATE_NAMES)), dtype=np.float64)
        for environment, coordinate in enumerate(coordinates):
            for substrate, spec in enumerate(self._substrate_layer_specs):
                result[environment, substrate] = self._substrate_layers[spec.name].get(
                    coordinate,
                    0.0,
                )
        return result

    def replace_substrate_batch_m(self, coordinates, concentrations_m):
        """Atomically replace sampled cells from a validated ``(E, 3)`` matrix."""
        coordinates = self._validate_substrate_coordinates(coordinates)
        values = self._validate_substrate_matrix(
            concentrations_m,
            expected_shape=(len(coordinates), len(SUBSTRATE_NAMES)),
            name='concentrations_m',
        )
        replacements = []
        for substrate, spec in enumerate(self._substrate_layer_specs):
            replacement = tuple(float(value) for value in values[:, substrate])
            if any(0.0 < value < spec.min_concentration_m for value in replacement):
                raise ValueError(
                    f'replacement {spec.name} concentration is below its retention floor'
                )
            replacements.append(replacement)
        for substrate, spec in enumerate(self._substrate_layer_specs):
            layer = self._substrate_layers[spec.name]
            for coordinate, value in zip(coordinates, replacements[substrate]):
                if value == 0.0:
                    layer.pop(coordinate, None)
                else:
                    layer[coordinate] = value

    def snapshot_substrate_layers(self):
        """Capture registration and every cell because diffusion may touch neighbors."""
        return SubstrateLayerSnapshot(
            layer_specs=self._substrate_layer_specs,
            layer_order=self.substrate_layer_names,
            layer_items=tuple(
                tuple(sorted(self._substrate_layers[name].items()))
                for name in self.substrate_layer_names
            ),
        )

    def restore_substrate_layers(self, snapshot):
        """Restore registration and values, preserving dict identities when possible."""
        if not isinstance(snapshot, SubstrateLayerSnapshot):
            raise TypeError('snapshot must be a SubstrateLayerSnapshot')
        if not snapshot.layer_specs:
            if snapshot.layer_order or snapshot.layer_items:
                raise ValueError('unregistered substrate snapshot must not contain layers')
            for layer in self._substrate_layers.values():
                layer.clear()
            self._substrate_layer_specs = ()
            self._substrate_layers = {}
            return

        specs = self._validate_substrate_specs(snapshot.layer_specs)
        if snapshot.layer_order != tuple(spec.name for spec in specs):
            raise ValueError('snapshot substrate order does not match its layer specs')
        if len(snapshot.layer_items) != len(specs):
            raise ValueError('snapshot substrate values do not match its layer specs')
        rebuilt = []
        for spec, items in zip(specs, snapshot.layer_items):
            layer = {}
            for coordinate, value in items:
                coordinate = self._coerce_substrate_coordinate(coordinate)
                if not math.isfinite(value) or value < spec.min_concentration_m:
                    raise ValueError('snapshot contains an invalid substrate concentration')
                if coordinate in layer:
                    raise ValueError('snapshot contains a duplicate substrate coordinate')
                layer[coordinate] = float(value)
            rebuilt.append(layer)
        if self._substrate_layer_specs == specs:
            for spec, layer in zip(specs, rebuilt):
                target = self._substrate_layers[spec.name]
                target.clear()
                target.update(layer)
        else:
            self._substrate_layer_specs = specs
            self._substrate_layers = {
                spec.name: layer for spec, layer in zip(specs, rebuilt)
            }

    def clear_substrate_layers(self):
        """Remove all continuous substrate values while retaining registration."""
        for layer in self._substrate_layers.values():
            layer.clear()

    def tick_substrates(self, dt_seconds):
        """Diffuse and decay the three FBA substrate fields on the Python path."""
        self._require_substrate_layers()
        if isinstance(dt_seconds, bool) or not isinstance(dt_seconds, Real):
            raise TypeError('dt_seconds must be a real number')
        dt_seconds = float(dt_seconds)
        if not math.isfinite(dt_seconds) or dt_seconds < 0.0:
            raise ValueError('dt_seconds must be finite and nonnegative')
        for spec in self._substrate_layer_specs:
            self._process_layer_python(
                self._substrate_layers[spec.name],
                float(spec.diffusion_um2_s),
                float(spec.decay_rate_s_inv),
                dt_seconds,
                min_concentration=spec.min_concentration_m,
            )

    def _require_substrate_layers(self):
        if not self._substrate_layer_specs:
            raise RuntimeError('continuous substrate layers are not registered')

    def _validate_substrate_coordinate(self, coordinate):
        self._require_substrate_layers()
        return self._coerce_substrate_coordinate(coordinate)

    @staticmethod
    def _coerce_substrate_coordinate(coordinate):
        if not isinstance(coordinate, tuple) or len(coordinate) != 2:
            raise TypeError('substrate coordinate must be a two-integer tuple')
        if any(isinstance(value, bool) or not isinstance(value, Integral) for value in coordinate):
            raise TypeError('substrate coordinate must be a two-integer tuple')
        return int(coordinate[0]), int(coordinate[1])

    @staticmethod
    def _validate_substrate_specs(specs):
        specs = tuple(specs)
        if len(specs) != len(SUBSTRATE_NAMES) or not all(
            isinstance(spec, SubstrateLayerSpec) for spec in specs
        ):
            raise TypeError('substrate layer specs must be SubstrateLayerSpec records')
        if tuple(spec.name for spec in specs) != SUBSTRATE_NAMES:
            raise ValueError(f'substrate layer names must exactly match {SUBSTRATE_NAMES}')
        for spec in specs:
            for name, value in (
                ('diffusion_um2_s', spec.diffusion_um2_s),
                ('decay_rate_s_inv', spec.decay_rate_s_inv),
                ('min_concentration_m', spec.min_concentration_m),
            ):
                if isinstance(value, bool) or not isinstance(value, Real):
                    raise TypeError(f'{spec.name}.{name} must be a real number')
                if not math.isfinite(float(value)) or float(value) < 0.0:
                    raise ValueError(f'{spec.name}.{name} must be finite and nonnegative')
            if not isinstance(spec.provenance, str) or not spec.provenance.strip():
                raise ValueError(f'{spec.name}.provenance must be a non-empty string')
        return specs

    def _validate_substrate_coordinates(self, coordinates):
        coordinates = tuple(self._validate_substrate_coordinate(value) for value in coordinates)
        if len(set(coordinates)) != len(coordinates):
            raise ValueError('substrate coordinates must be unique')
        return coordinates

    def _validate_substrate_matrix(self, values, expected_shape, name):
        import numpy as np

        self._require_substrate_layers()
        if not isinstance(values, np.ndarray) or values.dtype != np.dtype(np.float64):
            raise TypeError(f'{name} must be a float64 NumPy array')
        if values.shape != expected_shape:
            raise ValueError(f'{name} must have shape {expected_shape}')
        if not values.flags.c_contiguous:
            raise ValueError(f'{name} must be C-contiguous')
        if not np.all(np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError(f'{name} must contain finite nonnegative concentrations')
        return values

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

        if not callable(getattr(self._native_module, 'process_layers', None)):
            error = _NativeModuleUnavailable(
                "The 'bacteria_c_kernel' extension lacks process_layers; "
                "rebuild it or use backend='python'."
            )
            self._native_module = None
            self._native_import_error = error
            raise error

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

    def _process_layer_python(
        self,
        grid,
        diffusion_um2_s,
        decay_s,
        dt,
        min_concentration=ECM_MIN_CONCENTRATION,
    ):
        if not grid:
            return

        # This explicit Fickian stencil caps alpha = D * dt / dx^2 below 0.25 to preserve nonnegative mass; see DOCUMENTATION.md §4.9.
        alpha = min(ECM_ALPHA_MAX, max(0.0, diffusion_um2_s * dt / (self.cell_size ** 2)))
        decay_multiplier = max(0.0, min(1.0, 1.0 - decay_s * dt))
        next_grid = {}

        for k, val in list(grid.items()):
            if val < min_concentration:
                continue

            decayed = val * decay_multiplier
            retained = decayed * (1.0 - 4.0 * alpha)
            neighbor_share = decayed * alpha

            if retained >= min_concentration:
                next_grid[k] = next_grid.get(k, 0.0) + retained

            cx, cy = k
            if neighbor_share >= min_concentration:
                for nx, ny in ((cx - 1, cy), (cx + 1, cy), (cx, cy - 1), (cx, cy + 1)):
                    next_grid[(nx, ny)] = next_grid.get((nx, ny), 0.0) + neighbor_share

        grid.clear()
        grid.update({k: v for k, v in next_grid.items() if v >= min_concentration})

    def _process_layer(self, grid, diffusion_um2_s, decay_s, dt):
        return self._process_layer_python(grid, diffusion_um2_s, decay_s, dt)
