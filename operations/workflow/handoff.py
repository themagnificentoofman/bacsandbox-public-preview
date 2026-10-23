"""Admit the completed C01-C08 storage workflow into the coordinator.

Only operational state is written. Native receipt semantics are reused after an
exact source closure, using the original stable-file verifier. There are no
simulation, upload, deletion, lock-recovery, selection or resume-marker actions.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import time
import zipfile

import controller as c


def need(value, message):
    if not value:
        raise c.WorkflowError(message)


def pin(path):
    path = Path(path)
    c.check_binding({'path': str(path), 'sha256': c.digest(path)})
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': c.digest(path)}


def bound(binding):
    c.check_binding(binding)
    if 'bytes' in binding:
        need(Path(binding['path']).stat().st_size == binding['bytes'], 'Bound size changed')
    return c.read(binding['path'])


def load_native(binding, name):
    c.check_binding(binding)
    path = Path(binding['path'])
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def add(expected, binding):
    p = binding['path']
    row = {k: binding[k] for k in ('path', 'bytes', 'sha256')}
    need(p not in expected or expected[p] == row, 'Conflicting native source: ' + p)
    expected[p] = row


def shape(boundary):
    need(boundary['status'] == 'step1_complete_stopped_before_C08_recovery', 'Step 1 incomplete')
    need(boundary['pending_storage_transaction'] is None, 'Unsettled storage transaction')
    need(boundary['source']['scientific_sha'] == 'f8f196fccd879ce2d87daf1078a374ed20a86d64',
         'Wrong scientific checkpoint')
    need(boundary['C08']['state'] == 'incomplete_not_resumed', 'C08 state changed')


def collect(spec):
    """Build the exact union from already sealed native final records."""
    boundary = bound(spec['boundary'])
    shape(boundary)
    completion = bound(spec['completion'])
    need(completion['status'] == 'step1_complete_with_verified_retention_and_milestone_sync',
         'Missing final Step 1 completion')
    need(completion['final_boundary'] == boundary, 'Completion binds another boundary')
    expected, memberships, receipts, pairs = {}, {}, {}, []
    for binding in [spec['boundary'], spec['completion'], spec['adapter_tests'], *spec['native_helpers'],
                    spec['lifecycle_helper'], spec['prior_lifecycle'], spec['supervisor_journal'],
                    spec['scientific_python'], spec['identity_adapter']]:
        add(expected, binding)
    for candidate, sets, objects in [('C05', 10, 80), ('C06', 14, 90), ('C07', 10, 70)]:
        row = boundary['retention'][candidate]
        terminal = bound(row['terminal'])
        need((terminal['candidate'], terminal['sets'], terminal['arms'], terminal['uploaded_objects'])
             == (candidate, sets, 40, objects), 'Native retention denominator changed')
        need(terminal['status'] == 'complete_candidate_serial_retention_verified'
             and terminal['master'] == row['master'], 'Invalid native retention closure')
        for binding in [row['terminal'], row['master'], *terminal['public_provenance'],
                        *terminal['helper_sources']]:
            add(expected, binding)
        master = bound(row['master'])
        need(master['candidate_retention_complete'] is True and master['main_frozen'] is False,
             'Unexpected scientific status in retention record')
        coverage_path = master['exact_coverage_path']
        coverage_binding = pin(coverage_path)
        need(coverage_binding['sha256'] == master['exact_coverage_sha256'], 'Coverage binding differs')
        add(expected, coverage_binding)
        coverage = bound(coverage_binding)
        for path, record in coverage['files'].items():
            add(expected, {'path': path, **record})
        plan = bound(master['serial_plan'])
        report = bound(plan['source_verification'])
        memberships[report['cohort_path']] = sorted(report['cohort_files'])
        receipts[candidate + '_retention'] = row['terminal']
        review = bound(row['consolidation'])
        need(review['status'] == 'completed_verified' and review['pairs'] == 60,
             'Consolidation review incomplete')
        dedup = bound(review['completion'])
        need(dedup['status'] == 'completed' and dedup['candidate'] == candidate
             and len(dedup['linked']) == 60 and dedup['freeze_allowed_by_this_record'] is False,
             'Consolidation closure differs')
        journal = Path(review['completion']['path']).parent / 'journal.jsonl'
        journal_binding = pin(journal)
        need(journal_binding['sha256'] == dedup['journal_sha256'], 'Consolidation journal changed')
        need(c.read(dedup['plan_path'])['candidate'] == candidate, 'Consolidation plan mismatch')
        for binding in [row['consolidation'], review['completion'], journal_binding]:
            add(expected, binding)
        for path, sha in dedup['source_sha256'].items():
            binding = pin(path)
            need(binding['sha256'] == sha, 'Native consolidation source changed')
            add(expected, binding)
        pairs.extend(dedup['linked'])
        receipts[candidate + '_consolidation'] = row['consolidation']
    snapshot = boundary['C08']['snapshot_verification']
    need(snapshot['matched'] is True and snapshot['readback_sha256'] == snapshot['source']['sha256'],
         'C08 remote verification missing')
    add(expected, snapshot['source'])
    c.check_binding(snapshot['source'])
    with zipfile.ZipFile(snapshot['source']['path']) as z:
        manifest = json.loads(z.read('SNAPSHOT_MANIFEST.json'))
    need(manifest['cohort_files'] == 82 and manifest['resume_performed'] is False,
         'C08 snapshot scope differs')
    cohort = Path(boundary['source']['executing_root']) / 'evidence/grid/C08'
    members = []
    for row in manifest['files']:
        add(expected, {'path': row['source_path'], **row})
        if Path(row['source_path']).is_relative_to(cohort):
            members.append(str(Path(row['source_path']).relative_to(cohort)))
    memberships[str(cohort)] = sorted(members)
    receipts['C08_snapshot'] = spec['c08_verification']
    add(expected, spec['c08_verification'])
    need(bound(spec['c08_verification']) == snapshot, 'C08 proof differs')
    return boundary, expected, memberships, receipts, pairs


def fingerprint(path):
    p = Path(path)
    need(p.is_absolute() and p.resolve(strict=True) == p, 'Canonical path changed')
    s = p.lstat()
    need(stat.S_ISREG(s.st_mode), 'Nonregular source')
    return {k: getattr(s, k) for k in ('st_dev', 'st_ino', 'st_mode', 'st_nlink',
                                      'st_size', 'st_mtime_ns', 'st_ctime_ns')}


def check_members(memberships):
    for directory, expected in memberships.items():
        paths = list(Path(directory).rglob('*'))
        need(not any(p.is_symlink() for p in paths), 'Unexpected source symlink')
        actual = sorted(str(p.relative_to(directory)) for p in paths if p.is_file())
        need(actual == expected, 'Candidate source membership changed: ' + directory)


def verify_pairs(pairs):
    for pair in pairs:
        left, right = fingerprint(pair['U']), fingerprint(pair['H'])
        need((left['st_dev'], left['st_ino']) == (right['st_dev'], right['st_ino'])
             and left['st_size'] == right['st_size'] == pair['bytes'],
             'Completed replay consolidation changed')
        need(left['st_nlink'] == right['st_nlink'] == 2, 'Unexpected shared-inode links')
        need(not os.path.lexists(pair['backup']), 'Historical backup reappeared; do not remove it')


def verify_union(expected, native):
    """One native byte check per unique unchanged inode; preserve every path."""
    result, seen = {}, {}
    for name, row in sorted(expected.items()):
        before = fingerprint(name)
        key = (before['st_dev'], before['st_ino'])
        if key in seen:
            prior = seen[key]
            need(prior['fingerprint'] == before and prior['bytes'] == row['bytes']
                 and prior['sha256'] == row['sha256'], 'Shared inode has conflicting source proof')
            result[name] = {**row, 'fingerprint': before}
        else:
            verified = native.verify_file(name, row)
            result[name] = {**row, 'fingerprint': verified['fingerprint']}
            seen[key] = result[name]
    check_fingerprints(result)
    return result, len(seen)


def check_fingerprints(records):
    for path, row in records.items():
        need(fingerprint(path) == row['fingerprint'], 'Verified source changed: ' + path)


def check_process_snapshot(snap, prior, old_pids, forbidden):
    """Require the prior outer proc mount and persistent PID-1 anchor as well as absence."""
    need(snap['boot_id'] == prior['boot_id'], 'Boot changed; re-establish ownership')
    need(not snap['read_errors'] and not snap['old_pid_paths_present'], 'Unresolved legacy process')
    def mount(s):
        rows = [r for r in s['selected_mountinfo'] if r.split()[4] == '/proc']
        need(len(rows) == 1 and ' - proc ' in rows[0] and 'hidepid' not in rows[0]
             and 'subset=' not in rows[0], 'Incomplete process view')
        fields = rows[0].split()
        return fields[2:5]
    need(mount(snap) == mount(prior), 'Different outer proc filesystem')
    keys = ('outer_pid', 'name', 'outer_ppid', 'start_ticks', 'namespace_pid_chain')
    def anchor(s):
        rows = [r for r in s['processes'] if r['outer_pid'] == 1]
        need(len(rows) == 1, 'Missing outer process anchor')
        return {k: rows[0][k] for k in keys}
    need(anchor(snap) == anchor(prior), 'Outer PID-1 lifetime changed')
    own = snap['self']['outer_pid']
    for row in snap['processes']:
        need(row['outer_pid'] not in old_pids, 'Legacy owner still exists')
        need(row['outer_group'] != 46579 and row['outer_session'] != 46579,
             'Legacy worker group remains')
        need(row['outer_pid'] == own or not set(row.get('argv', [])).intersection(forbidden),
             'Another workflow owner is visible')


def ownership(spec):
    lifecycle = load_native(spec['lifecycle_helper'], 'native_lifecycle_reader')
    previous = bound(spec['prior_lifecycle'])['snapshots'][-1]
    journal = Path(spec['supervisor_journal']['path'])
    c.check_binding(spec['supervisor_journal'])
    events = [json.loads(line) for line in journal.read_text().splitlines()]
    need(events[-1]['event'] == 'paused', 'Legacy supervisor journal changed')
    old_pids = sorted({r['outer_pid'] for r in events[-1]['group_snapshot']}
                      | set(spec['additional_legacy_pids']))
    snap = lifecycle.snapshot(old_pids)
    check_process_snapshot(snap, previous, old_pids, set(spec['legacy_entrypoints']))
    need(not os.path.lexists(spec['resume_marker']), 'Resume marker exists')
    # Do not preserve arbitrary process command lines in a shareable checkpoint.
    for row in snap['processes']:
        row.pop('argv', None)
    snap['self'].pop('argv', None)
    return snap


def identity(spec):
    for binding in (spec['scientific_python'], spec['identity_adapter']):
        c.check_binding(binding)
    need(str(Path(spec['scientific_executable']).resolve(strict=True)) == spec['scientific_python']['path'],
         'Scientific environment executable resolution changed')
    result = subprocess.run([spec['scientific_executable'], '-B',
                             spec['identity_adapter']['path'], spec['boundary']['path']],
                            capture_output=True, text=True, check=False)
    need(result.returncode == 0, 'Scientific identity adapter failed: ' + result.stderr)
    return json.loads(result.stdout)


def admit(spec_path, destination):
    spec_path, destination = Path(spec_path).resolve(), Path(destination).resolve()
    spec_binding = pin(spec_path)
    spec = bound(spec_binding)
    need(spec['scope'] == 'step2_operational_handoff_only', 'Wrong admission scope')
    tests = bound(spec['adapter_tests'])
    need(tests['status'] == 'passed' and tests['failures'] == 0 and tests['skipped'] == 0,
         'Bound native adapter tests have not passed')
    for binding in tests['sources']:
        c.check_binding(binding)
    state = destination / 'STATE.json'
    destination.mkdir(parents=True, exist_ok=True)
    with c.owner(state):
        need(not (destination / 'ADMISSION.json').exists() and not state.exists(),
             'Admission exists; verify/resume instead')
        need(shutil.disk_usage(destination).free >= spec['minimum_free_bytes'], 'Insufficient space')
        native = load_native(spec['native_helpers'][0], 'native_stable_file_verifier')
        native.write_new(destination / 'ADMISSION_INTENT.json', {
            'status': 'admission_started_not_transferred', 'utc': c.now(), 'spec': spec_binding,
            'owner_state_path': str(state), 'C08_resumed': False,
            'failure_policy': 'Preserve this attempt; no automatic admission retry'})
        before = ownership(spec)
        scientific = identity(spec)
        boundary, expected, memberships, receipts, pairs = collect(spec)
        check_members(memberships)
        verify_pairs(pairs)
        records, unique = verify_union(expected, native)
        check_members(memberships)
        verify_pairs(pairs)
        after = ownership(spec)
        c.check_binding(spec_binding)
        check_fingerprints(records)
        admission = {
            'status': 'admitted_operational_handoff_C08_paused', 'utc': c.now(),
            'spec': spec_binding, 'source': pin(__file__), 'controller': pin(c.__file__),
            'scope': spec['scope'], 'owner_state_path': str(state),
            'previous_owner': spec['previous_owner'],
            'successor_owner': spec['workflow_id'],
            'ownership_observations': [before, after], 'scientific_identity': scientific,
            'verified_sources': records, 'unique_inodes_hashed': unique,
            'memberships': memberships, 'receipts': receipts, 'consolidated_pairs': pairs,
            'capacity_free_bytes': shutil.disk_usage(destination).free,
            'native_semantic_validation': 'Reuse completed native validations after exact source-hash closure',
            'C08_resumed': False, 'scientific_lock_ownership_transferred': False,
            'legacy_process_exit_reason': 'Unknown; repeated current outer-proc absence verified',
            'future_gate': 'Separate C08 recovery authorization, fresh identities/capacity and native lock recovery',
        }
        native.write_new(destination / 'ADMISSION.json', admission)
        return {'status': admission['status'], 'verified_paths': len(records),
                'unique_inodes_hashed': unique, 'admission': pin(destination / 'ADMISSION.json')}


def verify_admission(path, sha):
    admission = bound({'path': str(Path(path).resolve()), 'sha256': sha})
    need(admission['status'] == 'admitted_operational_handoff_C08_paused', 'Wrong admission')
    c.check_binding(admission['source'])
    c.check_binding(admission['controller'])
    spec = bound(admission['spec'])
    check_fingerprints(admission['verified_sources'])
    check_members(admission['memberships'])
    verify_pairs(admission['consolidated_pairs'])
    ownership(spec)
    need(not admission['C08_resumed'] and not admission['scientific_lock_ownership_transferred'],
         'Admission crosses scientific boundary')
    return admission


def make_plan(admission_path, output):
    admission_binding = pin(admission_path)
    admission = verify_admission(admission_path, admission_binding['sha256'])
    spec = bound(admission['spec'])
    base = [str(Path(sys.executable).resolve()), '-B', str(Path(__file__).resolve())]
    steps = []
    for name, receipt in admission['receipts'].items():
        steps.append({
            'id': name, 'kind': 'command', 'scope': 'import-completed-native-work',
            'argv': base + ['refuse-reexecution'],
            'verify_argv': base + ['verify-receipt', '--admission', str(admission_path),
                '--sha256', admission_binding['sha256'], '--receipt-id', name],
            'cwd': str(Path(output).parent), 'receipt_path': receipt['path'],
            'receipt_require': {'status': bound(receipt)['status']},
            'capacity_path': str(Path(output).parent), 'minimum_free_bytes': spec['minimum_free_bytes'],
        })
    steps.append({'id': 'C08_recovery_authorization', 'kind': 'gate',
                  'scope': 'scientific-stop', 'reason': admission['future_gate']})
    plan = {'schema_version': 1, 'workflow_id': spec['workflow_id'],
            'authorization_basis': spec['authorization_basis'],
            'allowed_scopes': ['import-completed-native-work', 'scientific-stop'],
            'bindings': [admission_binding, admission['spec'], admission['source'], admission['controller']],
            'stop_after': 'C08_recovery_authorization', 'steps': steps}
    c.validate(plan)
    native = load_native(spec['native_helpers'][0], 'native_plan_publisher')
    native.write_new(output, plan)
    return pin(output)


def verify_receipt(admission_path, sha, name):
    admission = verify_admission(admission_path, sha)
    need(name in admission['receipts'], 'Unknown native receipt')
    binding = admission['receipts'][name]
    bound(binding)
    return {'status': 'native_receipt_reused_after_verified_source_closure', 'receipt': binding,
            'native_semantics_recomputed': False, 'new_remote_readback': False}


def advance_imports(state_path, limit):
    """Bounded automatic import only; never retry a failed native action."""
    need(type(limit) is int and 1 <= limit <= 16, 'Import bound must be 1..16')
    for _ in range(limit):
        status = c.status(state_path)
        if status['phase'] in ('blocked', 'stopped'):
            break
        if status['phase'] in ('failed', 'in_flight'):
            need(status['last_event']['event'] == 'existing_receipt_requires_reconciliation',
                 'Interrupted verification/action requires explicit reconcile')
            c.reconcile(state_path)
        else:
            result = c.advance(state_path, 1)
            if result['phase'] == 'in_flight':
                c.reconcile(state_path)
            elif result['phase'] == 'ready':
                break  # Capacity hold; do not poll it in a tight loop.
    return c.status(state_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['admit', 'make-plan', 'verify-receipt',
                                         'verify-admission', 'refuse-reexecution', 'advance-imports'])
    parser.add_argument('--spec', type=Path)
    parser.add_argument('--destination', type=Path)
    parser.add_argument('--admission', type=Path)
    parser.add_argument('--sha256')
    parser.add_argument('--receipt-id')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--state', type=Path)
    parser.add_argument('--max-actions', type=int, default=8)
    args = parser.parse_args()
    if args.action == 'refuse-reexecution':
        parser.exit(2, 'Completed native work cannot be re-executed by this successor plan\n')
    if args.action == 'admit':
        result = admit(args.spec, args.destination)
    elif args.action == 'make-plan':
        result = make_plan(args.admission, args.output)
    elif args.action == 'verify-receipt':
        result = verify_receipt(args.admission, args.sha256, args.receipt_id)
    elif args.action == 'verify-admission':
        value = verify_admission(args.admission, args.sha256)
        result = {'status': value['status']}
    else:
        result = advance_imports(args.state, args.max_actions)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
