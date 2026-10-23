"""Reporting must never certify incomplete evidence or drive execution."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import report


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


class ReportingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cohort = self.root / 'C08'
        self.state = self.root / 'STATE.json'
        self.proc = self.root / 'proc'
        self.proc.mkdir()
        self.manifest = {'cohort': 'diagnostic_v03_C08', 'seeds': [4000],
                         'protocol_sha256': 'p', 'code': {'sha256': 'c'},
                         'pairs': [{'seed': 4000, 'order': ['U', 'H'],
                                    'inputs': {'U': 'u', 'H': 'h'}}]}
        self.manifest['manifest_sha256'] = report.digest(self.manifest)
        write(self.cohort / 'manifest.json', self.manifest)
        self.plan = {'workflow_id': 'test', 'stop_after': 'gate', 'steps': [
            {'id': 'execute', 'kind': 'command', 'capacity_path': str(self.root),
             'minimum_free_bytes': 0, 'receipt_path': str(self.root / 'EXECUTION.json'),
             'argv': ['/python', '-B', '/adapter.py', 'execute', '--config', '/config.json']},
            {'id': 'gate', 'kind': 'gate'}]}
        write(self.root / 'PLAN.json', self.plan)
        self.saved = {'workflow_id': 'test', 'phase': 'in_flight', 'next_index': 0,
                      'in_flight': 'execute', 'completed': {}, 'events': [],
                      'plan': {'path': str(self.root / 'PLAN.json'),
                               'sha256': hashlib.sha256((self.root / 'PLAN.json').read_bytes()).hexdigest()}}
        write(self.state, self.saved)
        self.add_process(100, ['/python', '/adapter.py', 'execute', '--config', '/config.json'])

    def add_process(self, pid, args, state='S'):
        directory = self.proc / str(pid)
        directory.mkdir(exist_ok=True)
        (directory / 'cmdline').write_bytes(b'\0'.join(a.encode() for a in args))
        (directory / 'stat').write_text(f'{pid} (python) {state} 1')

    def arm_path(self, arm='U'):
        return self.cohort / 'runs' / '4000' / arm

    def complete(self, arm='U'):
        result = {'seed': 4000, 'arm': arm, 'attempt': 1, 'status': 'completed',
                  'validated': True, 'cohort': self.manifest['cohort'],
                  'manifest_sha256': self.manifest['manifest_sha256'],
                  'input_sha256': arm.lower(), 'protocol_sha256': 'p', 'code_sha256': 'c'}
        path = self.arm_path(arm)
        write(path / 'attempt-001' / 'result.json', result)
        write(path / 'completed.json', {'attempt': 1, 'result_sha256': report.digest(result)})
        return result

    def sample(self, **kwargs):
        return report.observe(self.state, self.cohort, proc=self.proc, **kwargs)

    def test_marker_alone_is_not_completion(self):
        write(self.arm_path() / 'completed.json', {'attempt': 1, 'result_sha256': 'fake'})
        status, _ = self.sample()
        self.assertEqual(status['arms']['metadata_checked'], 0)
        self.assertIn('arm_metadata_unreadable', status['issues'])

    def test_metadata_is_not_native_verification(self):
        self.complete()
        status, _ = self.sample()
        self.assertEqual(status['arms']['metadata_checked'], 1)
        self.assertEqual(status['arms']['native_verified_cached'], 0)
        self.assertEqual(status['next_gate'], 'gate')

    def test_result_or_input_tampering_is_detected(self):
        result = self.complete()
        result['input_sha256'] = 'other'
        write(self.arm_path() / 'attempt-001' / 'result.json', result)
        write(self.arm_path() / 'completed.json', {'attempt': 1, 'result_sha256': report.digest(result)})
        self.assertEqual(self.sample()[0]['arms']['metadata_checked'], 0)

    def test_retained_lock_and_progress_do_not_claim_live_worker(self):
        write(self.arm_path() / 'active.lock', {'pid': 100})
        write(self.arm_path() / 'attempt-001' / 'progress.json', {'reached_time_s': 5})
        status, _ = self.sample()
        self.assertEqual(status['workers']['observed'], 0)
        self.assertEqual(status['workers']['claims'], 1)

    def test_process_identity_ignores_shell_text_and_status_command(self):
        self.add_process(100, ['/python', '/controller.py', 'status', '--state', str(self.state)])
        self.add_process(101, ['/bin/sh', '-c', '/adapter.py execute --config /config.json'])
        status, _ = self.sample()
        self.assertIn('ownership_unconfirmed', status['issues'])

    def test_exact_worker_match_and_stale_progress(self):
        write(self.arm_path() / 'active.lock', {'pid': 1})
        progress = self.arm_path() / 'attempt-001' / 'progress.json'
        write(progress, {'reached_time_s': 5})
        self.add_process(101, ['/python', '/adapter.py', 'arm', '--config', '/config.json',
                               '--seed', '4000', '--arm', 'U'])
        status, _ = self.sample(clock=lambda: progress.stat().st_mtime + 700)
        self.assertEqual(status['workers']['observed'], 1)
        self.assertIn('progress_stale', status['issues'])
        self.add_process(101, ['/python', '/adapter.py', 'arm', '--config', '/other.json',
                               '--seed', '4000', '--arm', 'U'])
        self.assertEqual(self.sample()[0]['workers']['observed'], 0)

    def test_missing_process_view_is_unknown(self):
        status, _ = report.observe(self.state, self.cohort, proc=self.root / 'absent')
        self.assertIsNone(status['workers']['observed'])
        self.assertIn('process_view_incomplete', status['issues'])

    def test_failed_arm_keeps_bounded_diagnostics(self):
        write(self.arm_path() / 'attempt-001' / 'result.json',
              {'status': 'failed', 'error': 'bad' * 10000})
        (self.root / 'arm-4000-U.log').write_text('irrelevant\n' * 10000 + 'native exception\n')
        status, detail = self.sample()
        packet = report.diagnose(status, detail)
        self.assertIn('native_arm_failed', status['issues'])
        self.assertIn('native exception', packet['excerpt'])
        self.assertLessEqual(len(packet['excerpt']), 1200)
        self.assertLessEqual(len(packet['arm_failure']['error']), 180)
        self.assertIn('reconcile', packet['side_effects'])

    def test_capacity_bound(self):
        class Space:
            free = -1
        with patch.object(report.shutil, 'disk_usage', return_value=Space()):
            self.assertIn('capacity_below_plan_bound', self.sample()[0]['issues'])

    def test_caches_only_bound_native_receipt(self):
        for arm in ('U', 'H'):
            self.complete(arm)
        path = self.root / 'EXECUTION.json'
        write(path, {'status': 'C08_native_40_arms_completed',
                     'records': [{'seed': 4000, 'arm': a,
                                  'sha256': report.read(self.arm_path(a) / 'completed.json')[0]['result_sha256']}
                                 for a in ('U', 'H')]})
        self.saved['completed']['execute'] = {'path': str(path),
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
        write(self.state, self.saved)
        self.assertEqual(self.sample()[0]['arms']['native_verified_cached'], 2)
        path.write_text('{}')
        status, _ = self.sample()
        self.assertEqual(status['arms']['native_verified_cached'], 0)
        self.assertIn('verified_receipt_unreadable_or_changed', status['issues'])

    def test_receipt_from_other_cohort_cannot_certify_current_results(self):
        for arm in ('U', 'H'):
            self.complete(arm)
        path = self.root / 'EXECUTION.json'
        write(path, {'status': 'C08_native_40_arms_completed',
                     'records': [{'seed': 4000, 'arm': a, 'sha256': 'other-cohort'} for a in ('U', 'H')]})
        self.saved['completed']['execute'] = {'path': str(path),
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
        write(self.state, self.saved)
        status, _ = self.sample()
        self.assertEqual(status['arms']['native_verified_cached'], 0)
        self.assertIn('native_receipt_result_mismatch', status['issues'])

    def test_malformed_or_changed_metadata_fails_closed(self):
        (self.root / 'PLAN.json').write_text('{}')
        status, _ = self.sample()
        self.assertEqual(status['phase'], 'unknown')
        self.assertEqual(status['issues'], ['observation_unavailable'])

    def test_no_workflow_writes_or_large_stream_reads(self):
        self.complete()
        stream = self.arm_path() / 'attempt-001' / 'trajectory.jsonl.gz'
        stream.write_bytes(b'not a gzip file')
        before = {str(p): p.read_bytes() for p in self.root.rglob('*') if p.is_file()}
        self.sample()
        after = {str(p): p.read_bytes() for p in self.root.rglob('*') if p.is_file()}
        self.assertEqual(before, after)
        self.assertLess(len(json.dumps(self.sample()[0])), 1400)


class WatchTests(unittest.TestCase):
    def run_watch(self, reports, timeout=5):
        now, outputs = [0], []
        sequence = iter(reports)
        latest = [reports[-1]]
        def sample():
            latest[0] = next(sequence, latest[0])
            return latest[0], {}
        def sleep(seconds):
            now[0] += seconds
        code = report.watch(sample, interval=1, timeout=timeout, clock=lambda: now[0],
                            sleep=sleep, emit=lambda text: outputs.append(json.loads(text)))
        return code, outputs

    def status(self, count=0, issues=None, phase='in_flight'):
        return {'phase': phase, 'operation': 'execute', 'issues': issues or [],
                'arms': {'metadata_checked': count, 'total': 40, 'native_verified_cached': 0}}

    def test_routine_arm_updates_do_not_emit(self):
        code, rows = self.run_watch([self.status(i) for i in range(6)])
        self.assertEqual(code, 0)
        self.assertEqual([r['event'] for r in rows], ['initial', 'timeout'])
        self.assertEqual(rows[-1]['status']['arms']['metadata_checked'], 5)

    def test_failure_returns_bounded_packet_after_confirmation(self):
        good, bad = self.status(), self.status(issues=['native_arm_failed'])
        code, rows = self.run_watch([good, bad, bad])
        self.assertEqual(code, 2)
        self.assertEqual(rows[-1]['event'], 'attention')

    def test_single_transient_error_does_not_trigger_attention(self):
        good, bad = self.status(), self.status(issues=['observation_unavailable'])
        code, rows = self.run_watch([good, bad, good, good])
        self.assertEqual(code, 0)
        self.assertEqual([r['event'] for r in rows], ['initial', 'timeout'])

    def test_completion_is_boundary_not_scientific_approval(self):
        code, rows = self.run_watch([self.status(), self.status(40), self.status(40)])
        self.assertEqual(code, 0)
        self.assertEqual(rows[-1]['event'], 'boundary')
        self.assertEqual(rows[-1]['packet']['status']['arms']['native_verified_cached'], 0)

    def test_initial_blocker_and_final_timeout_blocker_are_reported(self):
        bad = self.status(issues=['ownership_unconfirmed'])
        code, rows = self.run_watch([bad, bad])
        self.assertEqual(code, 2)
        self.assertEqual(rows[-1]['event'], 'attention')
        code, rows = self.run_watch([self.status(), bad], timeout=1)
        self.assertEqual(code, 2)
        self.assertIn('ownership_unconfirmed', rows[-1]['status']['issues'])


if __name__ == '__main__':
    unittest.main()
