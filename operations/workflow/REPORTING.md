# Compact arm reporting

User instruction, October 23, 2026: implement token-efficient execution. The
admitted native coordinator owns the complete queue and all safety/verification
checks. The model handles exceptions, meaningful boundaries, interpretation, and
requested status. Routine successful arms require zero model decisions.

`report.py` uses only the standard library and is independent of `controller.py`.
Run it from a separate checkout; do not replace files bound into an active plan.
It never launches workers, acquires or removes locks, signals processes, invokes
verifiers, retries actions, edits plans/state, or opens compressed event streams.
It writes only its bounded stdout. No persistent monitor is installed by these
commands; the watcher is a foreground process with a finite timeout.

## Commands

Use the state and cohort paths from the latest admitted checkpoint:

```bash
python operations/workflow/report.py status --state /absolute/STATE.json --cohort /absolute/grid/C08
python operations/workflow/report.py watch --state /absolute/STATE.json --cohort /absolute/grid/C08 --interval 30 --timeout 900
python operations/workflow/report.py diagnose --state /absolute/STATE.json --cohort /absolute/grid/C08
```

Status is one compact JSON object, normally about 200–350 tokens. It identifies
the candidate, phase/operation, arm counts, fresh process observation, newest
progress-file age, free space versus the current plan bound, issues and next gate.
Counts are derived from registered pairs, not a hardcoded seed list or arbitrary
directory count. Small JSON reads are capped at 8 MiB. No entire state, manifests,
process arguments, per-arm inventories or source hashes are printed.

Three completion counts have deliberately different meanings:

- `markers`: completion files observed. Existence alone proves little.
- `metadata_checked`: result hashes, status and registered identity match each
  marker. Output streams and duplicate attempts are not revalidated here.
- `native_verified_cached`: registered C08 records in a hash-bound execution
  receipt already accepted by the controller's native verifier. This is cached
  certification, not a fresh recheck of scientific outputs. It remains zero
  while the cohort verifier is pending, even if individual arms returned.

Fresh shared `/proc` argument-token matches identify the admitted script/config,
seed and arm; namespace-local PIDs in `active.lock` are never trusted as liveness.
Claims are reported separately. Missing or incomplete process access means unknown
liveness. A saved in-flight state with no observed owner requires reconciliation,
not automatic retry. Zombie processes are excluded. A slow/stale progress warning
is an observation to investigate, not permission to kill or restart a worker.

`watch` prints an initial summary, then keeps polling inside Python without
printing routine arm completions, time changes or normal capacity fluctuations.
It returns at a confirmed issue, phase/operation change, all metadata-complete
arms, native verification change, terminal phase, or timeout. Two observations
confirm events; the initial observation counts. Default reporting interval is
30 seconds and stale-progress threshold 600 seconds, adjustable independently of
native monitoring. These settings never change native capacity, ownership or
scientific checks. An event returns a bounded diagnostic packet. Exit 2 means
attention is needed; exit 0 can mean timeout and never authorizes a scientific
gate. Inspect the event and phase, not just the exit code.

`diagnose` includes the operation, first affected arm/attempt when available,
retained evidence paths, side-effect uncertainty, and at most 1,200 characters
from the last 12 log lines (at most 8 KiB read). It prefers the affected arm log,
then the verifier log, then the action log. Read additional evidence only when
this packet points to it. Never treat an incomplete read as a completed action.

## Continuation and handoff

1. Read the latest checkpoint and pending operation; reuse completed evidence.
2. Request one compact status or let one foreground watcher wait for a boundary.
   Do not ask the model to poll once per arm or run redundant watchers.
3. For an issue, retain the packet and reconcile through the admitted workflow.
   The helper does not take ownership or change recovery/admission requirements.
4. Update the PR and Drive at milestones, blockers and pause boundaries.

For the October 23 successor, the recorded paths are:
`/path/to/private-workspace` and
`/path/to/private-cohort`. Resolve these against the latest
handoff on continuation. The scientific source remains `f8f196fccd879ce2d87daf1078a374ed20a86d64`.
This reporting change has no dependency on merging PR #12, nor does it complete
C08, the full eligibility decision, retention, qualification, freeze or main run.

The reporter and coordinator require a surviving execution environment. Silence
from a watcher is not evidence that workers remain alive. A stopped chat/session
does not become durable compute through reduced polling.

## Verification

```bash
python -m unittest discover -s operations/workflow -p 'test_*.py' -v
```

Reporting tests cover false completion, changed inputs/receipts, stale locks,
process identity and incomplete visibility, capacity, bounded diagnostics,
read-only operation, quiet routine updates, transient reads, failure and boundary
events. The dedicated Linux CI workflow runs these tests with the controller and
handoff suite. Native-file integration tests retain their existing environment
requirements; the reporter tests need no simulator or scientific data.
