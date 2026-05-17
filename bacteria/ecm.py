import math
import os
import sys

if __package__ is None and not hasattr(sys, 'frozen'):
    path = os.path.realpath(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(path)))
    __package__ = 'bacteria'
    
from .constants import FPS

class ECM:
    """
    Extracellular Matrix managing diffusion and degradation of environmental molecules.
    Uses a coarse spatial hash (dictionary grid) for fast lookups and low overhead.
    """
    def __init__(self, cell_size=50):
        self.cell_size = cell_size
        self.qs_inducers = {}
        self.bacteriocins = {}
        
        self.chunk_size = 200
        self._qs_keys = []
        self._qs_idx = 0
        self._bac_keys = []
        self._bac_idx = 0

    def add(self, x, y, amount, layer='qs'):
        """
        Add molecules to a specific layer at (x, y).
        """
        cx = int(x // self.cell_size)
        cy = int(y // self.cell_size)
        target = self.qs_inducers if layer == 'qs' else self.bacteriocins
        target[(cx, cy)] = target.get((cx, cy), 0.0) + amount

    def sample(self, x, y, layer='qs'):
        """
        Sample the local concentration of a layer at (x, y).
        """
        cx = int(x // self.cell_size)
        cy = int(y // self.cell_size)
        target = self.qs_inducers if layer == 'qs' else self.bacteriocins
        return target.get((cx, cy), 0.0)

    def tick(self, dt):
        self._process_layer('qs', self.qs_inducers, 0.3, 0.05, dt)
        self._process_layer('bac', self.bacteriocins, 0.1, 0.1, dt)

    def _process_layer(self, layer_name, grid, base_diffuse, base_decay, dt):
        is_qs = layer_name == 'qs'
        keys = self._qs_keys if is_qs else self._bac_keys
        idx = self._qs_idx if is_qs else self._bac_idx
        
        if idx >= len(keys) or len(keys) == 0:
            keys = list(grid.keys())
            idx = 0
            if is_qs:
                self._qs_keys = keys
            else:
                self._bac_keys = keys
        
        if not keys:
            return
            
        sweep_frames = max(1.0, len(keys) / self.chunk_size)
        eff_dt = dt * sweep_frames
        
        diffuse_rate = base_diffuse * eff_dt
        decay_rate = base_decay * eff_dt
        
        chunk = keys[idx:idx + self.chunk_size]
        new_idx = idx + self.chunk_size
        
        for k in chunk:
            val = grid.get(k, 0.0)
            if val < 0.001:
                if val <= 0 and k in grid:
                    del grid[k]
                continue
                
            decayed = val * (1.0 - decay_rate)
            diffusion_amt = decayed * diffuse_rate
            retained = decayed - (diffusion_amt * 4.0)
            
            if retained < 0:
                retained = 0
            
            grid[k] = retained
            
            cx, cy = k
            for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
                grid[(nx, ny)] = grid.get((nx, ny), 0.0) + diffusion_amt
                
        if is_qs:
            self._qs_idx = new_idx
        else:
            self._bac_idx = new_idx
