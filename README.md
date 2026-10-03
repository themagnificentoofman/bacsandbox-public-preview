# Bacteria Sandbox

Bacteria Sandbox is an interactive 2D bacterial simulation built with Python and Pygame.

For the full technical model, equations, constants, CSV schemas, and implementation details, see [DOCUMENTATION.md](DOCUMENTATION.md).

## Setup

Use standard CPython 3.14. Development and CI use **3.14.7**, pinned in
[.python-version](.python-version); the package range is `>=3.14,<3.15`.
Free-threaded Python (`3.14t`) is outside this validation target.

1. Install Python 3.14.7 and create one project environment named `.venv`:

   macOS/Linux:

   ```bash
   python3.14 -m venv .venv
   source .venv/bin/activate
   ```

   Windows PowerShell:

   ```powershell
   py -3.14 -m venv .venv
   .\.venv\Scripts\Activate.ps1
   ```

2. Verify `python --version` reports `3.14.7`, then install dependencies:

   ```bash
   python -m pip install --upgrade pip
   python -m pip install -r requirements.txt
   python -m pip check
   ```

3. Run the app:

   ```bash
   python -m bacteria
   ```

Both install routes use `pygame-ce` (still imported as `pygame`) and the same
runtime dependencies from `requirements.txt`. Do not install the separate `pygame`
distribution alongside it. Native build dependencies are isolated by pip from
`pyproject.toml`. For an existing Python 3.12 installation, follow the
[migration and verification notes](DOCUMENTATION.md#21-python-314-migration-and-verification)
before rebuilding `.venv`. In VS Code, select this environment with **Python:
Select Interpreter**; an existing explicit selection overrides the workspace default.

### Optional Compiled Kernel

The normal requirements-only setup above needs no compiler and retains Python fallbacks for both ECM and contact physics. The optional unified extension accelerates ECM diffusion/decay and contact-pair processing and also exposes the reduced metabolic kernel through an opt-in direct API. The app's existing pg-based metabolism remains the default; the reduced API is not a UI mode or automatic `Simulation` integration and is not a genome-scale model. To build the extension in editable mode with automatic OpenMP detection, run the command for your shell from the repository root.

PowerShell:

```powershell
$env:BACSANDBOX_OPENMP = "auto"
python -m pip install -e .
```

POSIX shells:

```bash
BACSANDBOX_OPENMP=auto python -m pip install -e .
```

Use `BACSANDBOX_OPENMP=on` to require OpenMP or `BACSANDBOX_OPENMP=off` to build deliberately without it. OpenMP builds on macOS require `libomp`. Backend and build details are documented in [DOCUMENTATION.md §10.1](DOCUMENTATION.md#101-hybrid-c-computational-kernel), [§10.2](DOCUMENTATION.md#102-vectorized-physics-and-spatial-hash-kernel), and [§10.3](DOCUMENTATION.md#103-constraint-based-metabolic-caching-kernel).

For a macOS OpenMP build:

```bash
brew install libomp
CPPFLAGS="-I$(brew --prefix libomp)/include" \
LDFLAGS="-L$(brew --prefix libomp)/lib" \
ARCHFLAGS="-arch $(uname -m)" \
BACSANDBOX_OPENMP=on python -m pip install -e .
```

After building, verify the native extension and run the full suite:

```bash
python -m pip check
python -c "import bacteria_c_kernel as k; print('OpenMP:', k.openmp_enabled()); assert hasattr(k, 'MetabolicKernel')"
python -m unittest discover
```

## Using The App

The Pygame interface includes:

* Grid: pan/zoom simulation view with minimap, scale indicator, food nodes, bacteria, and hover inspection.
* Graph: population history charts with import/export controls.
* Stats: live population, elapsed time, tickrate, births, deaths, and per-species summaries.
* Species: editable species parameter table with duplicate/delete and undo/redo support.
* Globals: editable environment settings for oxygen, carbon dioxide, pH, carrying capacity, and food drops.

## Research Audits

Inspect species inputs before interpreting a run:

```bash
python audit_species_parameters.py --output parameter_audit.json
python -m unittest bacteria.test_parameter_audit bacteria.test_comparative_sweeps
```

The audit reports explicit/defaulted values and implementation diagnostics; it does
not infer measured traits. Comparative sweeps export four stable CSVs plus
`run_manifest.json`, recording realized founders, environment and analysis limits.
See [parameter auditing](DOCUMENTATION.md#51-parameter-inventory-and-missing-values)
and [comparative sweeps](DOCUMENTATION.md#105-comparative-spatial-control-sweep-harness).

For source-preserving Zhou benchmark preparation, see the
[data-readiness tools](research/data_readiness/README.md) and
[documented evidence limits](DOCUMENTATION.md#106-empirical-table-readiness).
The CSV audit and optional supplement extractors do not fit or validate the model.

## Project Files

* `species_data.csv`: species configuration loaded at startup and edited by the Species tab.
* `env_config.json`: global environment defaults loaded at startup and edited by the Globals tab.
* `simulation_stats.csv`: exported/imported simulation timeline, including elapsed seconds, tickrate, and cumulative global births/deaths.
* `DOCUMENTATION.md`: single source of truth for formulas, biological rationale, constants, UI behavior, architecture, and CSV details.
