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

class ECM:
    """
    Extracellular Matrix managing diffusion and degradation of environmental molecules.
    Uses a coarse spatial hash (dictionary grid) for fast lookups and low overhead.
    """
    def __init__(self, cell_size=ECM_CELL_SIZE_UM):
        self.cell_size = cell_size
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
        self._process_layer(self.qs_inducers, ECM_QS_DIFFUSION_UM2_S, ECM_QS_DECAY_S, dt)
        self._process_layer(self.bacteriocins, ECM_BAC_DIFFUSION_UM2_S, ECM_BAC_DECAY_S, dt)
        self._process_layer(self.acetate, ECM_ACETATE_DIFFUSION_UM2_S, ECM_ACETATE_DECAY_S, dt)
        self._process_layer(self.phage, ECM_PHAGE_DIFFUSION_UM2_S, ECM_PHAGE_DECAY_S, dt)

    def _process_layer(self, grid, diffusion_um2_s, decay_s, dt):
        if not grid:
            return

        # Explicit Fickian 2D stencil: alpha = D * dt / dx^2; alpha <= 0.25 preserves nonnegative mass.
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
