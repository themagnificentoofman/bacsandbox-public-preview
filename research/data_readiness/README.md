# Empirical table readiness

These optional tools prepare source-traceable observations for a proposed retrospective Zhou benchmark. They do not fit or evaluate the simulator. Scientific findings, bibliographic attribution and unresolved observation choices are in [DOCUMENTATION.md §10.6](../../DOCUMENTATION.md#106-empirical-table-readiness).

## Reproduce

Use the project's CPython 3.14.7 environment and install the optional analysis dependencies:

```sh
python -m pip install -r research/data_readiness/requirements-readiness.txt
python -m research.data_readiness.fetch_sources
python -m research.data_readiness.prepare_observations
python -m research.data_readiness.extract_s2_s13
python research/data_readiness/reproduce.py
python research/data_readiness/run_final_checks.py
```

The fetch command downloads seven public supplement files from the exact URLs in the two `provenance/*.json` manifests. It verifies size, SHA256 and MD5 before installing each file and refuses to overwrite differing existing bytes. `--verify-only` performs no network requests. If downloads are unavailable, place the unchanged originals under `raw/supplements/` (SI PDF and S1/S7/S8/S9) and `raw/supplements_s2_s13/` (S2/S13), then run:

```sh
python -m research.data_readiness.fetch_sources --verify-only
```

Originals and generated outputs are ignored by Git. Publisher papers, original XLSX files, copied author scripts, private correspondence and generated data exports are not included in this repository change. Source hashes and URLs are public provenance; obtain originals from their publisher/archive and retain their notices.

## Outputs and contracts

- `extract_zhou.py` writes S1C/S7/S8/S9 CSVs and audits under `results/extracted/`. Metadata inheritance is section-bounded and retains anchor cells. Measurement blanks, Excel errors, zeros and negative values remain explicit.
- `prepare_observations.py` rebuilds those inputs, then writes typed observations and descriptive diagnostics under `results/prepared/`. Observation keys preserve experimental context. S9A uses its provided log cells and retains the alternative linear values. It does not select an endpoint or success rule.
- `extract_s2_s13.py` writes genome annotations, target-identity records, identifier-collision diagnostics and diet chemistry under `results/s2_s13/`. It does not merge ambiguous genome entries, infer rates, convert food concentrations to dose, or join replicate labels across assays.
- `reproduce.py` inventories shipped configuration values in `results/parameter_ledger.csv`, retaining unknown empirical provenance.
- `run_final_checks.py` runs the generic parser, parameter/sweep and source-readiness suites, including `test_benchmark_mapping` and `test_feasibility`, together, followed by ordinary repository discovery, dependency checks and whitespace checks. It records the selected module list, counts and every skip reason under `results/`. Supply the pinned originals before running it. Standard project discovery covers the generic parser; source-specific tests are explicitly selected by this runner.

Extraction and preparation rebuild derived files. Original source bytes are read-only. The source-specific suites require openpyxl and verified originals; they are not silently replaced with synthetic measurements when originals are absent. Synthetic cases in parser/transfer tests exercise error handling only.

Optional published-data figures:

```sh
python -m research.data_readiness.plot_diagnostics
```

Use Matplotlib 3.10.8 to reproduce the reviewed figure. Plotted points retain every finite S1C/S9A observation; bars show available-value means, not confidence intervals or predictions. Mouse phases remain categorical source labels. Read the generated counts and the scientific limitations before interpreting a figure.

For an independently exported CSV, the standard-library audit accepts an explicit JSON contract:

```sh
python audit_benchmark_table.py table.csv --contract contract.json --metadata samples.csv --output new-audit.json
```

Required contract fields are `required_columns`, `key_columns`, `value_columns`, `missing_tokens`, `units`, `time_mapping` and `source_locator`. Include `join_columns` when supplying metadata. Existing audit outputs are refused. Exit 0 means structural checks passed; exit 2 reports structural/contract failure. Neither status establishes scientific readiness. S9B's stored Excel errors intentionally produce a failed structural audit while the source-preservation tests pass.

## Selected benchmark mapping

To generate the source-preserving ZHOU_BM001 mapping in a new directory:

```sh
python -m research.data_readiness.map_benchmark --output research/data_readiness/results/zhou_bm001_mapping
python -m unittest research.data_readiness.test_benchmark_mapping research.data_readiness.test_zhou_tables research.data_readiness.test_prepared_observations -v
```

This optional step requires the pinned S7/S8/S9 originals and `species_data.csv`.
The combined test command also requires the SI/S1 originals. Run from the repository root.
An existing output directory is refused. CSV and JSON outputs retain source cells,
metadata anchors, hashes, units, and distinct measurement states. `checksums.json`
covers the generated mapping files. A failed filesystem write may leave an
incomplete directory; verify its checksums and rerun to a new destination.

The exports contain the strain crosswalk, per-animal baseline observations,
selected target observations, endpoint eligibility, animal readiness, and a
categorical time crosswalk. A model label candidate is not a calibrated strain
profile. Inverse-log values are arithmetic quantities, not agent counts. Unresolved
logged zeros remain unconverted; the missing fecal endpoint remains unavailable.
Every animal retains `simulation_ready=false`. No simulator or parameter fitting
is invoked. Scientific scope and inclusion rules are documented in §10.6.

`mapping_audit.json` schema 2 scopes `simulation_ready=false` to the source-mapping
artifact: it does not assess a separately implemented executable or current project
completion. The former `remaining_requirements` field is replaced by
`external_execution_requirements`, which names contracts to verify outside this
mapping, plus `inference_limitations` and `abm_followup_only`. Agent rounding and
ecological seed sweeps are not continuous-reference prerequisites. Unknown grouping
limits inference rather than descriptive reporting, which has no biological
pass/fail threshold. Consumers of schema 1 must migrate these status fields; source
observations, units, identities and endpoint masks are unchanged.

## Minimal model feasibility

```sh
python -m research.data_readiness.audit_feasibility --output research/data_readiness/results/zhou_bm001_feasibility
python -m unittest research.data_readiness.test_feasibility -v
```

This optional audit uses the project's runtime dependencies plus the pinned
S7/S8/S9 originals. It invokes actual configuration interfaces and a lifecycle
policy in isolation; it never advances the ecological simulator or fits parameters.
It writes diet-interface probes, configured-profile provenance, complete-baseline
agent-resolution diagnostics, literal-day tick arithmetic, and input hashes to a
new directory using the mapping writer's non-overwrite/checksum contract.

Ambiguous baseline vectors are excluded only from the resolution diagnostic;
they remain in the benchmark cohort. Expected agent quotas are not founder counts.
The example budgets are diagnostic values, not selected simulation settings.
Time arithmetic is neither a selected biological clock nor a wall-time estimate.
The emitted `benchmark_execution_ready` flag remains false. Scientific interpretation,
minimal-model proposals and implementation boundaries are in §10.6.

The audit resolves distribution metadata from the imported pygame implementation:
`pygame-ce` for the supported environment, or `pygame` for a legacy installation.
It records module and distribution versions and discloses co-installation or
version disagreement instead of attributing a run to the other provider. Missing
metadata for the actual imported provider remains an error. The combined optional
source/mapping/feasibility suite (40 tests) and the read-only audit command were
verified with CPython 3.14.7, pygame-ce 2.5.8, NumPy 2.5.3 and openpyxl 3.1.5.
This is source and engineering verification; ordinary hosted CI does not execute
the optional workbook suites.

## Decisions before evaluation

ZHOU_BM001 preparation now selects a deterministic, continuous 18-label,
well-mixed reference and error-reporting-only evaluation. This reference is not
implemented or assessed by these source-preparation tools. Only iteration-1
unsupplemented control follow-ups may fit baseline target behavior; glucomannan
parameters require independent evidence
or predeclared uncertainty sets. Treated iteration-1 and all iteration-2 follow-up
values are evaluation-only, including for indirect choices of exposure, loss,
yield, zero handling or scenario selection. The iteration-1 between-diet contrast
has a calibrated control component and is not wholly independent validation.

The primary continuous reference has no agent founder rounding or ecological seed
sweep. Keep both source-derived timing scenarios and the exact 39/38 endpoint
masks. Report separate endpoint/iteration errors and comparator differences without
a biological pass/fail threshold. The executable protocol must separately specify
resource/biomass scales, independent response constraints or predeclared uncertainty
sets, zero handling, finite parameter bounds, and executable numerical checks.
These source-preparation outputs do not certify that freeze. Unknown cage
grouping limits inference but need not block a descriptive comparison. No fitting,
model prediction, mentor approval or empirical validation is implied by these tools.
