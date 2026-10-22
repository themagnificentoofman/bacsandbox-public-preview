"""Serial, receipt-based workflow controller; does not import the simulator.

Native commands own scientific validation, remote readback and cleanup safety.
This POSIX controller adds ownership, crash recovery and bounded scheduling.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


class WorkflowError(RuntimeError):
    """An operation requires reconciliation or corrected prerequisites."""


def now():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic(path, value):
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(value, f, indent=2, sort_keys=True)
            f.write('\n')
            f.flush()
            os.fsync(f.fileno())
        os.replace(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


@contextmanager
def owner(state_path):
    """Children inherit this lock, so a dead parent cannot permit overlap."""
    path = Path(str(state_path) + '.lock')
    with path.open('a+') as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise WorkflowError('Another controller or owned child holds the lock') from exc
        f.seek(0)
        f.truncate()
        json.dump({'pid': os.getpid(), 'acquired_utc': now()}, f)
        f.flush()
        try:
            yield f.fileno()
        finally:
            # Closing our descriptor, rather than LOCK_UN, preserves the child lock.
            pass


def check_binding(binding):
    path = Path(binding['path'])
    if not path.is_absolute() or path.is_symlink() or not path.is_file():
        raise WorkflowError('Missing/nonregular binding: ' + str(path))
    if digest(path) != binding['sha256']:
        raise WorkflowError('Binding changed: ' + str(path))


def validate(plan):
    if plan.get('schema_version') != 1 or not plan.get('workflow_id'):
        raise WorkflowError('Unsupported plan')
    if not plan.get('authorization_basis') or not plan.get('allowed_scopes'):
        raise WorkflowError('Explicit scope and existing authorization basis required')
    ids = [s['id'] for s in plan['steps']]
    if len(ids) != len(set(ids)) or any(not re.fullmatch(r'[A-Za-z0-9_-]+', x) for x in ids):
        raise WorkflowError('Unique safe step IDs required')
    if plan['stop_after'] not in ids:
        raise WorkflowError('Explicit stop boundary required')
    if not plan.get('bindings'):
        raise WorkflowError('Bound operational/scientific identities required')
    for s in plan['steps']:
        if s['scope'] not in plan['allowed_scopes']:
            raise WorkflowError('Step exceeds authorized scope: ' + s['id'])
        if s.get('kind') == 'gate':
            if not s.get('reason'):
                raise WorkflowError('Gate reason required')
            continue
        if s.get('kind') != 'command':
            raise WorkflowError('Unsupported step kind')
        for key in ('argv', 'verify_argv'):
            a = s.get(key)
            if not isinstance(a, list) or not a or not all(isinstance(x, str) and x for x in a):
                raise WorkflowError('Explicit argv arrays required')
            if not Path(a[0]).is_absolute():
                raise WorkflowError('Absolute executable path required')
        for key in ('cwd', 'receipt_path', 'capacity_path'):
            if not Path(s[key]).is_absolute():
                raise WorkflowError('Absolute paths required')
        if not s.get('receipt_require') or type(s['minimum_free_bytes']) is not int or s['minimum_free_bytes'] < 0:
            raise WorkflowError('Native receipt contract and capacity bound required')


def initialize(plan_path, state_path):
    plan_path, state_path = Path(plan_path).resolve(), Path(state_path).resolve()
    with owner(state_path):
        if state_path.exists():
            raise WorkflowError('State exists; resume it')
        plan = read(plan_path)
        validate(plan)
        for binding in plan['bindings']:
            check_binding(binding)
        atomic(state_path, {
            'schema_version': 1, 'workflow_id': plan['workflow_id'],
            'plan': {'path': str(plan_path), 'sha256': digest(plan_path)},
            'controller_sha256': digest(__file__), 'next_index': 0,
            'phase': 'ready', 'completed': {}, 'events': [], 'updated_utc': now(),
        })
    return status(state_path)


def loaded(state_path):
    state = read(state_path)
    check_binding(state['plan'])
    if state['controller_sha256'] != digest(__file__):
        raise WorkflowError('Controller changed; create an explicitly migrated checkpoint')
    plan = read(state['plan']['path'])
    validate(plan)
    if plan['workflow_id'] != state['workflow_id']:
        raise WorkflowError('Workflow identity differs')
    return plan, state


def save(state_path, state, event, **detail):
    state['updated_utc'] = now()
    state['events'].append({'utc': state['updated_utc'], 'event': event, **detail})
    atomic(state_path, state)


def prerequisites(plan, state):
    for binding in plan['bindings']:
        check_binding(binding)
    # Small operational receipts, never repeated full scans of scientific outputs.
    for binding in state['completed'].values():
        check_binding(binding)


def invoke(argv, step, log, lock_fd):
    with Path(log).open('xb') as f:
        result = subprocess.run(argv, cwd=step['cwd'], stdout=f, stderr=subprocess.STDOUT,
                                pass_fds=(lock_fd,), check=False)
        f.flush()
        os.fsync(f.fileno())
    if result.returncode:
        raise WorkflowError('Native command failed; inspect ' + str(log))


def verify(step, log, lock_fd):
    # The native verifier must independently check side effects and identity, not
    # merely accept the existence of a success-shaped JSON file.
    invoke(step['verify_argv'], step, log, lock_fd)
    receipt = read(step['receipt_path'])
    if any(receipt.get(k) != v for k, v in step['receipt_require'].items()):
        raise WorkflowError('Native receipt does not match contract')
    return {'path': step['receipt_path'], 'sha256': digest(step['receipt_path'])}


def finish(state_path, plan, state, step, receipt, event):
    state['completed'][step['id']] = receipt
    state['next_index'] += 1
    state.pop('in_flight', None)
    state['phase'] = 'stopped' if step['id'] == plan['stop_after'] else 'ready'
    save(state_path, state, event, step=step['id'], receipt=receipt)


def advance(state_path, max_actions=16):
    if type(max_actions) is not int or not 1 <= max_actions <= 100:
        raise WorkflowError('Batch bound must be 1..100')
    state_path = Path(state_path).resolve()
    with owner(state_path) as lock_fd:
        plan, state = loaded(state_path)
        if state['phase'] in ('in_flight', 'failed'):
            raise WorkflowError('Ambiguous/failed operation: reconcile; do not retry')
        for _ in range(max_actions):
            if state['phase'] in ('stopped', 'blocked'):
                break
            if plan.get('activation_blocker'):
                state['phase'] = 'blocked'
                save(state_path, state, 'activation_blocked', reason=plan['activation_blocker'])
                break
            prerequisites(plan, state)
            step = plan['steps'][state['next_index']]
            if step['kind'] == 'gate':
                state['phase'] = 'blocked'
                save(state_path, state, 'gate', step=step['id'], reason=step['reason'])
                break
            free = shutil.disk_usage(step['capacity_path']).free
            if free < step['minimum_free_bytes']:
                save(state_path, state, 'capacity_hold', free_bytes=free,
                     required_bytes=step['minimum_free_bytes'])
                break
            # Never rerun an externally completed/partially attempted command.
            if Path(step['receipt_path']).exists():
                state['phase'] = 'in_flight'
                state['in_flight'] = step['id']
                save(state_path, state, 'existing_receipt_requires_reconciliation', step=step['id'])
                break
            logs = state_path.parent / (state_path.name + '.logs')
            logs.mkdir(exist_ok=True)
            state['phase'] = 'in_flight'
            state['in_flight'] = step['id']
            save(state_path, state, 'started', step=step['id'], free_bytes=free)
            try:
                invoke(step['argv'], step, logs / (step['id'] + '.run.log'), lock_fd)
                receipt = verify(step, logs / (step['id'] + '.verify.log'), lock_fd)
                prerequisites(plan, state)
                finish(state_path, plan, state, step, receipt, 'verified')
            except Exception as exc:
                state['phase'] = 'failed'
                save(state_path, state, 'failed', step=step['id'], error=str(exc))
                raise
    return status(state_path)


def reconcile(state_path):
    """Accept already-completed work only through its original native verifier."""
    state_path = Path(state_path).resolve()
    with owner(state_path) as lock_fd:
        plan, state = loaded(state_path)
        if state['phase'] not in ('in_flight', 'failed'):
            raise WorkflowError('No interrupted operation to reconcile')
        prerequisites(plan, state)
        step = plan['steps'][state['next_index']]
        if step['kind'] != 'command' or state.get('in_flight') != step['id']:
            raise WorkflowError('In-flight identity mismatch')
        logs = state_path.parent / (state_path.name + '.logs')
        logs.mkdir(exist_ok=True)
        log = logs / (step['id'] + '.reconcile-' + str(len(state['events'])) + '.log')
        try:
            receipt = verify(step, log, lock_fd)
            prerequisites(plan, state)
        except Exception as exc:
            save(state_path, state, 'reconciliation_failed', error=str(exc))
            raise
        finish(state_path, plan, state, step, receipt, 'reconciled_without_retry')
    return status(state_path)


def status(state_path):
    # No lock, mutation, process signal, expensive hash pass or remote request.
    state = read(state_path)
    return {k: state[k] for k in ('workflow_id', 'phase', 'next_index', 'completed', 'updated_utc')} | {
        'in_flight': state.get('in_flight'), 'last_event': state['events'][-1] if state['events'] else None,
        'status_is_cached_observation': True,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['init', 'status', 'advance', 'reconcile'])
    p.add_argument('--state', type=Path, required=True)
    p.add_argument('--plan', type=Path)
    p.add_argument('--max-actions', type=int, default=16)
    a = p.parse_args()
    try:
        if a.action == 'init':
            if a.plan is None:
                p.error('init requires --plan')
            result = initialize(a.plan, a.state)
        elif a.action == 'advance':
            result = advance(a.state, a.max_actions)
        else:
            result = globals()[a.action](a.state)
        print(json.dumps(result, indent=2))
    except (WorkflowError, OSError, ValueError, KeyError) as exc:
        p.exit(2, str(exc) + '\n')


if __name__ == '__main__':
    main()
