# Bacteria Sandbox — public preview

Bacteria Sandbox (BacSandbox) is a Python/Pygame sandbox for exploring how individual simulated bacteria move, find nutrients, interact, and change population abundance over time. Optional C++ kernels accelerate selected calculations. This repository lets a reviewer run the software, inspect its architecture, and follow its development history.

**Start here:** [reviewer guide and architecture](PORTFOLIO.md) · [technical documentation](DOCUMENTATION.md) · [history and disclosure notes](HISTORY.md) · [verification](VERIFICATION.md)

This is a privacy-reviewed copy of the development repository. Historical commits retain their chronological development and ancestry, with the documented sanitization described in `HISTORY.md`. The runnable demo uses synthetic inputs; biological predictive accuracy has not been established.

## Run the demonstration

Use standard **CPython 3.14.7**, the source repository's pinned runtime. The package supports `>=3.14,<3.15`; free-threaded Python is outside the stated target. Run the following from the repository root:

```sh
python3.14 -m venv .venv
source .venv/bin/activate
python --version
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip check
python examples/portfolio_demo.py
```

On Windows PowerShell, create and activate the environment with:

```powershell
py -3.14 -m venv .venv
.\.venv\Scripts\Activate.ps1
```

Then use the same `python -m pip ...` and demo commands above. The requirements-only installation needs no C++ compiler. It installs `pygame-ce`, which is imported as `pygame`; do not also install the separate `pygame` distribution.

The demo runs the existing simulation engine for 20 simulated seconds, with two invented populations and illustrative inputs. It shows movement, nutrient patches, population counts, births, and deaths. Space pauses/resumes it; Escape exits. Cells are enlarged for visibility. It uses a fixed seed, fixed time steps, and explicit Python backends. Its populations are not measured biological species, and this presentation wrapper is separate from the registered research experiments.

For a bounded run without a display, optionally saving a new population timeline:

```sh
python examples/portfolio_demo.py --headless --seconds 5 --output demo-populations.csv
```

The script refuses to overwrite an existing output file. A matching seed reproduces the example within the same code/runtime; it is not a promise of bitwise parity across platforms, dependency versions, or backends.

## Explore the full application

```sh
python -m bacteria
```

Click **Play** to start. Use **Grid** to view agents and nutrients, **Graph** to inspect population trajectories, **Stats** for births/deaths and elapsed time, **Species** to inspect/edit model parameters, and **Globals** to change environment settings. The application reads and can save `species_data.csv` and `env_config.json` in the working directory. Its export action produces a simulation timeline; use a disposable checkout if you want to keep the shipped settings untouched.

The desktop app requires a graphical display. This repository is an executable source demonstration, not a hosted browser app. A short recording of it can accompany the Maker Portfolio; the reviewer should not have to infer motion from a still screenshot.

## Inspect the implementation

| Component | Entry point |
| --- | --- |
| Individual agents, lifecycle, spatial simulation | [`bacteria/core_sim.py`](bacteria/core_sim.py) |
| Desktop interface and editing controls | [`bacteria/app.py`](bacteria/app.py), [`bacteria/ui_components.py`](bacteria/ui_components.py) |
| Extracellular fields and contact physics | [`bacteria/ecm.py`](bacteria/ecm.py), [`bacteria/physics.py`](bacteria/physics.py) |
| Optional reduced metabolic model and live adapter | [`bacteria/metabolism.py`](bacteria/metabolism.py), [`bacteria/metabolic_integration.py`](bacteria/metabolic_integration.py) |
| Dietary input translation and model assets | [`bacteria/diet_translator.py`](bacteria/diet_translator.py), [`bacteria/models/`](bacteria/models/) |
| Optional C++17 kernels | [`bacteria/src/`](bacteria/src/) |
| Parameter audit and comparative runs | [`audit_species_parameters.py`](audit_species_parameters.py), [`run_comparative_sweeps.py`](run_comparative_sweeps.py) |
| Research preparation and separate workflows | [`research/`](research/) |

The default desktop model uses discrete nutrient nodes and glucose-equivalent reserve bookkeeping. The opt-in reduced metabolic API is a separate path; it does **not** become active merely because the C++ extension is installed, and it is **not** a genome-scale reconstruction.

## Optional native build and checks

A C++17 toolchain is required for the editable native build. From the repository root on macOS/Linux:

```sh
BACSANDBOX_OPENMP=auto python -m pip install -e .
python -c "import bacteria_c_kernel as k; print('OpenMP:', k.openmp_enabled()); assert hasattr(k, 'MetabolicKernel')"
```

On PowerShell:

```powershell
$env:BACSANDBOX_OPENMP = "auto"
python -m pip install -e .
```

`auto` probes OpenMP and falls back to a serial native build if unavailable. Use `on` to require it or `off` to disable it. See [the technical documentation](DOCUMENTATION.md#101-hybrid-c-computational-kernel) for toolchain details.

Run the simulator's test modules with:

```sh
python -m unittest discover -s bacteria -p 'test_*.py'
```

Native-only checks may skip when the extension is absent. Separate source-dependent research tests require optional analysis dependencies and the exact publisher supplement workbooks; those originals are not distributed here. See [the data-readiness instructions](research/data_readiness/README.md). Do not interpret skipped tests or successful software checks as completed biological validation.

## Scope, history, and credit

Zachary Zheng's project development, including AI-assisted work, is represented by the retained history. [The reviewer guide](PORTFOLIO.md#contributions-and-assistance) explains the contribution boundary; [the history notes](HISTORY.md) explain the public-copy changes. Publication of this source does not assign a new open-source license to the project or change third-party licenses.

Private mentor correspondence, internal storage links, operational checkpoints, and retained execution evidence are excluded as documented. Research code and historical technical documentation remain useful for inspection, but references to excluded records do not make those records part of this preview or establish the current state of an experiment.
