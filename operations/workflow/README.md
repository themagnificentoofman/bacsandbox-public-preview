# Checkpoint-first BacSandbox workflows

Standing user instruction adopted October 22, 2026. Apply this policy in good
faith to implementation, execution, retention, verification and handoffs within
BacSandbox. It is the default way of working, not an amendment to an experiment.

## Continue from state

1. Read the latest durable checkpoint, its pending transaction and the operating
   procedure. Resolve observations against native receipts, not chat recollection.
2. Check relevant source/runtime/configuration identities, ownership, capacity and
   changes since that checkpoint. Reuse completed simulations, diagnostic traces,
   verified remote objects and reviewed tools. A new conversation is not itself
   grounds for another broad audit.
3. Advance authorized routine actions in bounded batches with one owner. Update
   the checkpoint at each verified operation, and stop at the recorded boundary
   or a real blocker. Record the reason if wider verification becomes necessary.
4. At milestones, blockers and pauses, publish a brief handoff and PR status
   pointing to the checkpoint. Preserve the exact tools as well as the state.
   A chat summary or Git SHA alone may omit temporary operational dependencies.

Implementation, execution, diagnostic review, remote retention, local cleanup,
scientific qualification and freeze are distinct statuses. Do not report one as
proof of another. Successful upload is not verified retention; retention does not
authorize removing every local input. Restore and verify canonical paths before
tools that require those paths run. An old release receipt does not authorize
deleting a file that later reappears there.

Use a measured, bounded archive queue: package one set, upload, read back, verify,
record receipts, release only eligible exact copies, then advance. Preserve the
existing reserve and admission formulas, including auxiliary and metadata space.
Do not add concurrency when it consumes the headroom required for verification.
Routine confirmations are unnecessary within existing authority; a changed scope,
scientific gate, or uncertain side effect must be resolved on its own facts.

## Controller contract

`controller.py` is a standard-library-only POSIX coordinator, separate from the
pinned scientific source. It runs explicit argument arrays for existing native
adapters. The adapters continue to own full scientific checks, exact remote byte
and identity checks, cleanup eligibility and current-byte/inode checks. The
controller neither approves a new cleanup plan nor substitutes receipt existence
for native validation.

It provides:

- an atomic, fsynced JSON state with operation history and a pinned plan;
- a nonblocking writer lock inherited by child commands;
- a durable `in_flight` record before each side effect;
- bounded sequential actions, capacity admission and a fixed stop boundary;
- success only after the native verification command and receipt contract pass;
- read-only cached status, with no process signals or remote writes;
- reconciliation of interrupted or externally completed work without re-execution.

Use a single shared state path for each workflow. The lock is local to that path
and host; it does not coordinate a legacy controller or another machine. Legacy
ownership must be explicitly handed off before activation. Commands must remain
foreground and retain the inherited descriptor; daemonizing adapters are outside
this contract. The controller does not supply a durable compute host or survive a
host loss. Store its checkpoint and receipts durably at pause/milestone boundaries.

### Plan format

The JSON plan has `schema_version: 1`, a stable `workflow_id`,
`authorization_basis`, `allowed_scopes`, `bindings`, `stop_after`, and ordered
`steps`. Bindings use absolute `path` and `sha256`. Bind the scientific identity
records, all adapter/helper sources, runtime identity evidence, and configuration
inputs. Do not include credentials or signed URLs. A plan is immutable after
initialization; activation or a changed plan requires an explicitly documented
migration, a new state, and preservation of the old state and receipts.

Each command step has the following fields:

| Field | Contract |
| --- | --- |
| `id`, `scope`, `kind` | Unique safe ID, authorized scope, `command` |
| `argv` | Absolute executable and explicit arguments; no implicit shell |
| `verify_argv` | Read-only native verifier that checks the actual side effects |
| `cwd` | Absolute working directory |
| `receipt_path` | Absolute native receipt path, separate from state/plan/logs |
| `receipt_require` | Required top-level receipt values, including native success status |
| `capacity_path`, `minimum_free_bytes` | Measured admission bound before the command; native per-object checks still apply |

A `gate` step instead has `id`, `scope`, `kind: gate`, and `reason`. It blocks;
there is no force/skip command. An optional plan-level `activation_blocker`
prevents all commands from executing. The checked-in transition plan uses this.

```bash
python operations/workflow/controller.py init --plan PLAN.json --state STATE.json
python operations/workflow/controller.py advance --state STATE.json --max-actions 16
python operations/workflow/controller.py status --state STATE.json
python operations/workflow/controller.py reconcile --state STATE.json
```

`advance` reuses completed receipts. A failure or uncertain completion stops the
batch. `reconcile` invokes only the native verifier; it never repeats the action.
If the native operation did not finish, preserve the original attempt and use the
existing recovery procedure to create a new attempt and amended plan. Missing or
ambiguous upload responses require remote identity reconciliation before any new
upload. Reusing a known upload must continue readback, not rebuild its archive.

Status is explicitly a cached observation. Before actual mutations, the controller
checks bound files and prior small receipts; the native adapters recheck relevant
current evidence. Plan authors must avoid expensive redundant full scans in each
adapter. This optimization never removes a required final gate.

### Tests and adoption

```bash
python -m unittest discover -s operations/workflow -p 'test_*.py' -v
```

Tests use disposable fixtures. They exercise interruption recovery, no duplicate
effects, concurrent ownership, changed bindings, verification failure, capacity
holds, bounded batches and stop boundaries. They do not run a scientific arm or
release production evidence.

The generic coordinator is implemented. Live C01–C08 command adapters are **not
activated**: step 1 must finish with the current procedure, then its native
commands/verifiers and original receipts must be bound into an admitted plan.
No existing upload/readback/cleanup format has been replaced mid-transaction.
The transition checkpoint and operational-source bundle make this pending work
explicit. See [RESUME.md](RESUME.md).
