#!/usr/bin/env python3
"""Verify that BacSandbox's integrated optimization features are present.

This file used to rewrite source files with one-off string patches. Those patches
are now part of the maintained codebase, so running this script should never edit
project files or reintroduce stale mechanics.
"""

from bacteria.core_sim import Bacterium, FoodNode, FoodSpatialHash, Simulation
from bacteria.ecm import ECM


def verify_integrated_optimizations():
    """Return a checklist proving the optimization hooks are available."""
    return {
        'food_node_pool': hasattr(FoodNode, 'acquire') and hasattr(FoodNode, 'release'),
        'bacterium_pool': hasattr(Bacterium, 'acquire') and hasattr(Bacterium, 'release'),
        'food_spatial_hash': FoodSpatialHash is not None,
        'simulation_viewport_update': 'viewport_bounds' in Simulation.update.__code__.co_varnames,
        'ecm_full_tick_update': not hasattr(ECM(), 'chunk_size'),
    }


if __name__ == '__main__':
    checks = verify_integrated_optimizations()
    for name, ok in checks.items():
        print(f'{name}: {"ok" if ok else "missing"}')
    if not all(checks.values()):
        raise SystemExit(1)