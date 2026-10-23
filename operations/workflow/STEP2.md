# Step 2: operational coordinator handoff

The Step 1 checkpoint supersedes the original transition observation. The
successor imports C05, C06 and C07 retention and consolidation separately, plus
the verified C08 filesystem snapshot. Its next step is a fixed C08 recovery
authorization gate. Importing completed work has no simulation, upload, deletion,
resume-marker, lock-recovery, selection or freeze action.

`handoff.py` binds the original completion receipts, native helper sources,
registered scientific identity and exact retained-source union. It uses the
original stable-file verifier once per unchanged inode, preserving every path.
It reuses prior native semantic validation after a fresh source-hash closure;
it does not claim new remote readback or repeat semantic validation. Reconciliation
checks the bound proof and current file identity (device, inode, mode, link count,
size, mtime and ctime), cohort membership and replay pairs. Any change requires
explicit investigation and a new admission. Allocation block counts are not
content identities.

Ownership observation uses the original lifecycle reader. It requires continuity
of boot ID, the outer proc filesystem and its persistent PID-1 lifetime, an
error-free process scan, absence of the recorded worker group and known owners,
and absence of a resume marker. A short process list in an isolated namespace
does not establish closure. Absence does not identify the exit cause or time.
C08's existing `active.lock` files remain preserved and are not transferred or
removed by this operational handoff.

`activate.py` claims one append-only owner record in the previous owner's
workspace before initializing the coordinator. It binds the reviewed repository
commit, test evidence, admission and shared owner record into the immutable plan.
An identical interrupted claim can continue; a different state/implementation
cannot silently take ownership. The coordinator's POSIX lock remains the
single-writer mechanism for its state and foreground verifier children.

## Procedure

1. Pin the current Step 1 checkpoint, completion receipt and reviewed native
   helpers in an admission specification. Keep the registered virtual-environment
   executable path while binding its resolved interpreter bytes.
2. Run disposable tests with the native verifier available. Record their actual
   results and source hashes. Run the read-only admission into a new attempt
   directory. Preserve a failed attempt; do not overwrite its intent.
3. Commit the adapters on a separate branch and verify that remote commit.
   Record its source hashes in the implementation receipt.
4. Run `activate.py --admission … --implementation … --tests …`. It records the
   shared owner, creates a new immutable plan, initializes the coordinator and
   imports the seven existing receipts through their bound verifiers. Existing
   work can never fall through to re-execution: its action entrypoint refuses.
5. Save the state, admission, source/version inventory and short resume note;
   synchronize the Drive handoff and PR at this milestone.

Status is read-only:

```bash
python operations/workflow/controller.py status --state /absolute/path/STATE.json
```

For an interrupted import, use `controller.py reconcile` after inspecting the
saved in-flight record. `handoff.py advance-imports` handles only newly discovered
existing receipts automatically. It never automatically retries a failed action
or failed verifier. Batches are bounded at 1–16 imports.

## Next scientific boundary

A later C08 request needs an explicitly admitted successor plan. Recheck the
registered scientific source/runtime/configuration, at least 8 GiB free, current
process ownership and native lock-recovery conditions. Preserve every partial
attempt; use the same saved inputs and seeds for a required restart. The
filesystem snapshot and `progress.json` are not mid-arm execution checkpoints.
The import-only plan does not qualify future C08 execution/retention adapters.
All scientific selection, numerical qualification and freeze gates still apply.

The recorded operational checkpoint is host-specific. After a host or filesystem
identity change, restore retained evidence and explicitly re-admit; do not edit
fingerprints or remove the gate to make an old state run.
