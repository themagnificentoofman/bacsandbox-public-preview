"""Claim one shared operational owner and activate a tested import-only plan."""
import argparse
from pathlib import Path

import controller as c
import handoff as h


def claim(registry, value, publish):
    """Exclusive, append-only fence; repeated identical claims are read-only."""
    registry = Path(registry)
    if registry.exists():
        h.need(c.read(registry) == value, 'Another owner is recorded; explicit migration required')
    else:
        publish(registry, value)
    return h.pin(registry)


def activate(admission_path, implementation_path, tests_path):
    admission_path = Path(admission_path).resolve()
    implementation_path, tests_path = Path(implementation_path).resolve(), Path(tests_path).resolve()
    admission_binding = h.pin(admission_path)
    implementation = c.read(implementation_path)
    tests = c.read(tests_path)
    h.need(implementation['status'] == 'remote_commit_verified'
           and len(implementation['commit']) == 40, 'Verified repository implementation required')
    h.need(tests['status'] == 'passed' and tests['failures'] == 0 and tests['skipped'] == 0,
           'Activation qualification missing')
    for row in implementation['sources'] + tests['sources']:
        c.check_binding(row)
    h.need(h.pin(__file__) in implementation['sources'] and h.pin(__file__) in tests['sources'],
           'Activation source not qualified and repository-bound')
    admission = h.verify_admission(admission_path, admission_binding['sha256'])
    registry = Path(admission['previous_owner']) / 'C01_C08_SUCCESSOR_OWNER.json'
    with c.owner(registry):
        admission = h.verify_admission(admission_path, admission_binding['sha256'])
        spec = h.bound(admission['spec'])
        native = h.load_native(spec['native_helpers'][0], 'native_activation_publisher')
        state = Path(admission['owner_state_path'])
        directory = state.parent
        owner = {'status': 'operational_owner_claimed_C08_paused',
                 'admission': admission_binding, 'state_path': str(state),
                 'implementation': h.pin(implementation_path), 'tests': h.pin(tests_path),
                 'scientific_lock_ownership_transferred': False}
        owner_binding = claim(registry, owner, native.write_new)
        base, plan_path = directory / 'PLAN_BASE.json', directory / 'PLAN.json'
        if not base.exists():
            h.make_plan(admission_path, base)
        plan = c.read(base)
        h.need(plan['workflow_id'] == spec['workflow_id']
               and plan['bindings'][0] == admission_binding, 'Base plan belongs to another admission')
        plan['bindings'] += [owner_binding, h.pin(implementation_path), h.pin(tests_path), h.pin(__file__)]
        c.validate(plan)
        if plan_path.exists():
            h.need(c.read(plan_path) == plan, 'Prepared plan changed; preserve and reconcile')
        else:
            native.write_new(plan_path, plan)
        if state.exists():
            previous_plan, _ = c.loaded(state)
            h.need(previous_plan == plan, 'Existing state binds another plan')
        else:
            c.initialize(plan_path, state)
        result = h.advance_imports(state, 8)
        h.need(result['phase'] == 'blocked' and result['next_index'] == 7
               and len(result['completed']) == 7
               and result['last_event']['step'] == 'C08_recovery_authorization',
               'Import did not reach its scientific stop boundary')
        return {'status': 'step2_complete_coordinator_at_C08_gate',
                'state': h.pin(state), 'plan': h.pin(plan_path), 'registry': owner_binding,
                'implementation_commit': implementation['commit'],
                'receipts_imported': list(result['completed']), 'C08_resumed': False,
                'pending_operation': result['in_flight']}


if __name__ == '__main__':
    import json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--admission', type=Path, required=True)
    p.add_argument('--implementation', type=Path, required=True)
    p.add_argument('--tests', type=Path, required=True)
    args = p.parse_args()
    print(json.dumps(activate(args.admission, args.implementation, args.tests), indent=2))
