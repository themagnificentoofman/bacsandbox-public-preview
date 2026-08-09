# Bacteria Sandbox

Bacteria Sandbox is an interactive 2D bacterial simulation built with Python and Pygame.

For the full technical model, equations, constants, CSV schemas, and implementation details, see [DOCUMENTATION.md](DOCUMENTATION.md).

## Setup

1. Install Python 3.
2. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

3. Run the app:

   ```bash
   python -m bacteria
   ```

### Optional Compiled Kernel

The normal requirements-only setup above needs no compiler and retains Python fallbacks for both ECM and contact physics. The optional unified extension accelerates ECM diffusion/decay and contact-pair processing when it is installed. To build it in editable mode with automatic OpenMP detection, run the command for your shell from the repository root.

PowerShell:

```powershell
$env:BACSANDBOX_OPENMP = "auto"
python -m pip install -e .
```

POSIX shells:

```bash
BACSANDBOX_OPENMP=auto python -m pip install -e .
```

Use `BACSANDBOX_OPENMP=on` to require OpenMP or `BACSANDBOX_OPENMP=off` to build deliberately without it. OpenMP builds on macOS require `libomp`. Backend and build details are documented in [DOCUMENTATION.md §10.1](DOCUMENTATION.md#101-hybrid-c-computational-kernel) and [§10.2](DOCUMENTATION.md#102-vectorized-physics-and-spatial-hash-kernel).

## Using The App

The Pygame interface includes:

* Grid: pan/zoom simulation view with minimap, scale indicator, food nodes, bacteria, and hover inspection.
* Graph: population history charts with import/export controls.
* Stats: live population, elapsed time, tickrate, births, deaths, and per-species summaries.
* Species: editable species parameter table with duplicate/delete and undo/redo support.
* Globals: editable environment settings for oxygen, carbon dioxide, pH, carrying capacity, and food drops.

## Project Files

* `species_data.csv`: species configuration loaded at startup and edited by the Species tab.
* `env_config.json`: global environment defaults loaded at startup and edited by the Globals tab.
* `simulation_stats.csv`: exported/imported simulation timeline, including elapsed seconds, tickrate, and cumulative global births/deaths.
* `DOCUMENTATION.md`: single source of truth for formulas, biological rationale, constants, UI behavior, architecture, and CSV details.
