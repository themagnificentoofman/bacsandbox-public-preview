"""Read-only, bounded reporting for the admitted spatial workflow.

Never imports the simulator, takes ownership, runs verifiers, retries, or writes
to the workflow. Native verification and safety monitoring retain their gates.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time


MAX_JSON = 8 * 1024 * 1024
READ_ERRORS = (OSError, ValueError, KeyError, TypeError, IndexError, AttributeError, OverflowError)


def read(path):
    with Path(path).open('rb') as stream:
        raw = stream.read(MAX_JSON + 1)
    if len(raw) > MAX_JSON:
        raise ValueError('metadata exceeds reporting bound')
    return json.loads(raw), raw


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def clipped(value, limit=180):
    return str(value).replace('\x1b', '')[:limit]


def option(args, flag):
    try:
        return args[args.index(flag) + 1]
    except (ValueError, IndexError):
        return None


def processes(state_path, step, proc=Path('/proc')):
    """Fresh argument-token matches in the shared process view; never use lock PIDs.

    Namespace-local lock PIDs can refer to different processes. Restrict matches
    to the exact admitted script/config or controller/state arguments instead.
    """
    workers, owners, unreadable = set(), 0, False
    argv = step.get('argv', [])
    config = option(argv, '--config') or option(argv, '--spec')
    script = next((a for a in argv if a.endswith('.py')), None)
    try:
        entries = list(proc.iterdir())
    except OSError:
        return workers, owners, False
    for entry in entries:
        if not entry.name.isdigit():
            continue
        try:
            with (entry / 'cmdline').open('rb') as stream:
                raw = stream.read(16385)
            if len(raw) > 16384:
                unreadable = True
                continue
            args = [a.decode(errors='replace') for a in raw.split(b'\0') if a]
            if not args:
                continue
            if not Path(args[0]).name.startswith('python'):
                continue
            stat = (entry / 'stat').read_text().rsplit(')', 1)[1].split()
            if stat[0] in ('Z', 'X'):
                continue
            controller_script = next((a for a in args if Path(a).name == 'controller.py'), None)
            controller = (controller_script and option(args, controller_script) in ('advance', 'reconcile')
                          and option(args, '--state') == str(state_path))
            adapter = (script and script in args and config
                       and config in (option(args, '--config'), option(args, '--spec'))
                       and option(args, script) in (option(argv, script), 'arm'))
            if adapter and option(args, '--seed') and option(args, '--arm'):
                workers.add((int(option(args, '--seed')), option(args, '--arm')))
            elif controller or adapter:
                owners += 1
        except FileNotFoundError:
            continue  # Process exited during observation.
        except (OSError, ValueError, IndexError):
            unreadable = True
    return workers, owners, not unreadable


def snapshot(state_path, cohort, *, stale_seconds=600, proc=Path('/proc'), clock=time.time):
    state_path, cohort = Path(state_path).resolve(), Path(cohort).resolve()
    now = clock()
    state, _ = read(state_path)
    plan, raw = read(state['plan']['path'])
    if (hashlib.sha256(raw).hexdigest() != state['plan']['sha256']
            or plan['workflow_id'] != state['workflow_id']):
        raise ValueError('plan identity differs from saved state')
    manifest, _ = read(cohort / 'manifest.json')
    unsigned = dict(manifest)
    unsigned.pop('manifest_sha256', None)
    if digest(unsigned) != manifest['manifest_sha256']:
        raise ValueError('manifest hash mismatch')
    pairs = manifest['pairs']
    seeds = [p['seed'] for p in pairs]
    if (len(seeds) != len(set(seeds)) or set(seeds) != set(manifest['seeds'])
            or any(type(s) is not int for s in seeds)
            or any(sorted(p['order']) != ['H', 'U'] for p in pairs)):
        raise ValueError('ambiguous registered arm inventory')
    steps = plan['steps']
    index = state['next_index']
    if type(index) is not int or not 0 <= index <= len(steps):
        raise ValueError('invalid next operation index')
    step = steps[index] if index < len(steps) else {}
    live, owners, process_view_complete = processes(state_path, step, proc)
    issues, failures, updates = [], [], []
    result_hashes = {}
    complete, markers, claimed, running = 0, 0, 0, 0
    for pair in pairs:
        for arm in pair['order']:
            seed = pair['seed']
            arm_dir = cohort / 'runs' / str(seed) / arm
            marker = arm_dir / 'completed.json'
            try:
                if (arm_dir / 'active.lock').exists():
                    claimed += 1
                if marker.exists():
                    markers += 1
                    completion, _ = read(marker)
                    attempt = completion['attempt']
                    if type(attempt) is not int or attempt < 1:
                        raise ValueError('invalid completed attempt')
                    result_path = arm_dir / f'attempt-{attempt:03d}' / 'result.json'
                    result, _ = read(result_path)
                    expected = {'status': 'completed', 'validated': True, 'seed': seed,
                                'arm': arm, 'attempt': attempt, 'cohort': manifest['cohort'],
                                'manifest_sha256': manifest['manifest_sha256'],
                                'protocol_sha256': manifest['protocol_sha256'],
                                'code_sha256': manifest['code']['sha256'],
                                'input_sha256': pair['inputs'][arm]}
                    if (digest(result) != completion['result_sha256']
                            or any(result.get(k) != v for k, v in expected.items())):
                        raise ValueError('completion metadata mismatch')
                    complete += 1
                    result_hashes[(seed, arm)] = completion['result_sha256']
                    updates.append(marker.stat().st_mtime)
                    continue
                attempts = sorted(arm_dir.glob('attempt-*'))
                if not attempts:
                    continue
                latest = attempts[-1]
                result_path = latest / 'result.json'
                if result_path.exists():
                    result, _ = read(result_path)
                    if result.get('status') == 'failed':
                        failures.append({'seed': seed, 'arm': arm, 'attempt': latest.name,
                                         'error': clipped(result.get('error')),
                                         'result': str(result_path)})
                        issues.append('native_arm_failed')
                progress = latest / 'progress.json'
                if progress.exists():
                    updated = progress.stat().st_mtime
                    updates.append(updated)
                    if (seed, arm) in live and now - updated > stale_seconds:
                        issues.append('progress_stale')
                elif (seed, arm) in live and now - latest.stat().st_mtime > stale_seconds:
                    issues.append('progress_missing')
                if (seed, arm) in live:
                    running += 1
            except READ_ERRORS as exc:
                issues.append('arm_metadata_unreadable')
                failures.append({'seed': seed, 'arm': arm, 'error': clipped(exc),
                                 'path': str(arm_dir)})
    capacity_path = step.get('capacity_path', str(cohort))
    free = shutil.disk_usage(capacity_path).free
    required = step.get('minimum_free_bytes')
    if required is not None and free < required:
        issues.append('capacity_below_plan_bound')
    phase = state['phase']
    if phase == 'in_flight' and not owners:
        issues.append('ownership_unconfirmed')
    if not process_view_complete:
        issues.append('process_view_incomplete')
    if phase in ('failed', 'blocked'):
        issues.append('workflow_' + phase)
    last = state.get('events', [])[-1:]
    if last and last[0]['event'] == 'capacity_hold':
        issues.append('controller_capacity_hold')
    # Receipts here are small, cached certifications, not re-validation of streams.
    certified = 0
    for binding in state.get('completed', {}).values():
        try:
            receipt, raw = read(binding['path'])
            if hashlib.sha256(raw).hexdigest() != binding['sha256']:
                raise ValueError('receipt changed')
            records = receipt.get('records', [])
            identities = {(r['seed'], r['arm']) for r in records}
            expected = {(p['seed'], a) for p in pairs for a in p['order']}
            if (receipt.get('status') == 'C08_native_40_arms_completed'
                    and identities == expected and len(records) == len(expected)
                    and all(result_hashes.get((r['seed'], r['arm'])) == r.get('sha256')
                            and r.get('sha256') is not None for r in records)):
                certified = len(records)
            elif receipt.get('status') == 'C08_native_40_arms_completed':
                issues.append('native_receipt_result_mismatch')
        except READ_ERRORS:
            issues.append('verified_receipt_unreadable_or_changed')
    gates = [s['id'] for s in steps[index:] if s['kind'] == 'gate']
    report = {
        'observed_utc': datetime.fromtimestamp(now, timezone.utc).isoformat(),
        'cohort': clipped(manifest['cohort'], 80), 'phase': clipped(phase, 40),
        'operation': clipped(state.get('in_flight') or step.get('id', 'none'), 80),
        'arms': {'metadata_checked': complete, 'total': len(pairs) * 2,
                 'native_verified_cached': certified, 'markers': markers},
        'workers': {'observed': running if process_view_complete else None,
                    'matched': running, 'claims': claimed,
                    'ownership': 'observed' if owners else 'unconfirmed'},
        'latest_progress_age_s': round(max(0, now - max(updates)), 1) if updates else None,
        'free_GiB': round(free / 1024**3, 2),
        'minimum_GiB': round(required / 1024**3, 2) if required is not None else None,
        'issues': sorted(set(issues)),
        'next_gate': clipped(gates[0] if gates else plan['stop_after'], 80),
        'verification': 'cached native receipts; metadata checks only; no stream revalidation',
    }
    return report, {'state': state, 'plan': plan, 'step': step, 'failures': failures,
                    'state_path': str(state_path), 'cohort_path': str(cohort)}


def observe(*args, **kwargs):
    try:
        return snapshot(*args, **kwargs)
    except READ_ERRORS as exc:
        return {'phase': 'unknown', 'issues': ['observation_unavailable'],
                'error': clipped(exc)}, {}


def tail(path, limit=1200):
    try:
        with Path(path).open('rb') as stream:
            stream.seek(0, 2)
            stream.seek(max(0, stream.tell() - 8192))
            text = stream.read(8192).decode(errors='replace')
        return clipped('\n'.join(text.splitlines()[-12:]), limit)
    except OSError:
        return None


def diagnose(report, detail):
    state_path = detail.get('state_path')
    step = detail.get('step', {})
    packet = {'status': report, 'operation': step.get('id'),
              'side_effects': 'possible; reconcile through admitted native verifier before retry',
              'state_path': state_path, 'receipt_path': step.get('receipt_path')}
    if detail.get('failures'):
        packet['arm_failure'] = detail['failures'][0]
    if state_path and step.get('id'):
        logs = Path(state_path + '.logs')
        paths = [logs / (step['id'] + suffix) for suffix in ('.run.log', '.verify.log')]
        packet['log_paths'] = [str(p) for p in paths]
        # Prefer a failing arm's native log, then verification, then action log.
        failure = packet.get('arm_failure', {})
        if 'seed' in failure:
            paths.append(Path(state_path).parent / f"arm-{failure['seed']}-{failure['arm']}.log")
        for path in reversed(paths):
            excerpt = tail(path)
            if excerpt:
                packet['excerpt_path'], packet['excerpt'] = str(path), excerpt
                break
    return packet


def event_key(report):
    """Routine arm/progress/free-space changes deliberately do not wake the model."""
    arms = report.get('arms', {})
    return (report.get('phase'), report.get('operation'), tuple(report.get('issues', [])),
            arms.get('metadata_checked', -1) == arms.get('total', -2),
            arms.get('native_verified_cached', 0))


def watch(sample, *, interval=30, timeout=900, emit=print, clock=time.monotonic,
          sleep=time.sleep):
    start = clock()
    initial, detail = sample()
    emit(json.dumps({'event': 'initial', 'status': initial}, separators=(',', ':')))
    baseline = event_key(initial)
    pending = None
    while True:
        left = timeout - (clock() - start)
        if left <= 0:
            emit(json.dumps({'event': 'timeout', 'status': initial}, separators=(',', ':')))
            return 2 if initial.get('issues') else 0
        sleep(min(interval, left))
        current, detail = sample()
        key = event_key(current)
        material = (key != baseline or bool(current.get('issues'))
                    or current.get('phase') in ('stopped', 'blocked', 'failed')
                    or key[3])
        # Two consecutive observations prevent transient atomic-file/process races
        # from producing a false blocker. The initial sample counts as one.
        if material and (pending == key or key == baseline):
            event = 'attention' if current.get('issues') else 'boundary'
            emit(json.dumps({'event': event, 'packet': diagnose(current, detail)},
                            separators=(',', ':')))
            return 2 if event == 'attention' else 0
        pending = key if material else None
        initial = current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['status', 'watch', 'diagnose'])
    parser.add_argument('--state', type=Path, required=True)
    parser.add_argument('--cohort', type=Path, required=True)
    parser.add_argument('--interval', type=float, default=30)
    parser.add_argument('--timeout', type=float, default=900)
    parser.add_argument('--stale-seconds', type=float, default=600)
    args = parser.parse_args()
    if not (0 < args.interval <= 60 and 0 < args.timeout <= 86400
            and 0 < args.stale_seconds < float('inf')):
        parser.error('interval: (0,60]; timeout: (0,86400]; stale-seconds: positive finite')
    sample = lambda: observe(args.state, args.cohort, stale_seconds=args.stale_seconds)
    if args.action == 'watch':
        return watch(sample, interval=args.interval, timeout=args.timeout,
                     emit=lambda text: print(text, flush=True))
    report, detail = sample()
    print(json.dumps(diagnose(report, detail) if args.action == 'diagnose' else report,
                     separators=(',', ':')))
    return 2 if report.get('issues') else 0


if __name__ == '__main__':
    sys.exit(main())
