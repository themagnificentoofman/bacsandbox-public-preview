# Preview verification — October 26, 2026

This is engineering/release verification, not biological validation.

| Check | Result |
| --- | --- |
| Source identity | Development master `ce64ecd0fe8e206b5a4758ff61139a282056e22b` |
| Source-object integrity | All 417 distinct source blobs match their Git SHA |
| Filtered history | 83 source commits represented; all 15 merges and ordered parents retained |
| Historical attribution | Author/committer names and timestamp instants retained; contact addresses sanitized; coauthor/bot credit retained |
| Retained historical content | Exact bytes/modes, except two private execution paths in `operations/workflow/REPORTING.md` |
| Exclusions | Checked in every historical tree; excluded-only objects and original unsanitized commits absent from the public object database |
| Git object integrity | `git fsck --full --strict` passed |
| Core Python fallback tests | 337 discovered; 288 passed, 49 skipped; no failures |
| Core native OpenMP tests | 337 discovered; 332 passed, 5 skipped; no failures |
| Synthetic demo | Repeated seed-2026 runs produced identical CSV bytes: 5 simulated seconds, 300 steps, 47 final agents, 7 births, 0 deaths |
| GUI smoke checks | Synthetic viewer rendered with SDL dummy display; original App initialized and handled quit |
| Dependency compatibility | Installed requirements passed the dependency check |

Test host: Linux x86_64, standard CPython 3.14.7. Runtime versions included
pygame-ce 2.5.8, NumPy 2.5.4, pandas 3.1.0, and Matplotlib 3.11.2.
The native build used pybind11 3.1.0 and OpenMP.

The core suite command was `python -m unittest discover`. Separate research
suites requiring publisher supplement workbooks were not rerun for this
portfolio preparation. Optional scale/benchmark/native/environment tests can
skip; a skip is not a passing check. No registered experiment was resumed,
qualified, or scientifically frozen by this work. No fresh macOS or Windows
execution or interactive human GUI review is claimed here.

The public release adds presentation documentation and a synthetic wrapper.
Tests exercise the existing implementation; they do not authenticate every
historical scientific claim or establish biological predictive accuracy.

Independent filtered-history verification passed 19,053 checks. The private integrated operational suite passed 80 tests (one environment-dependent skip); the public operational suite did not produce a complete terminal report and is not claimed as passing. All final PR-head hosted checks passed before merging. The public presentation retains research infrastructure and excludes standalone execution reports.
