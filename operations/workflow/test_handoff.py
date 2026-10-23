"""Disposable recovery/admission tests; never mutate production evidence."""
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import controller as c
import handoff as h
import test_controller


class HandoffTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.calls = []

    def native_verify(self, path, row):
        self.calls.append(path)
        self.assertEqual(c.digest(path), row['sha256'])
        self.assertEqual(Path(path).stat().st_size, row['bytes'])
        return {'fingerprint': h.fingerprint(path)}

    def pair(self):
        u, v = self.root / 'U', self.root / 'H'
        u.write_bytes(b'preserved native bytes')
        os.link(u, v)
        return u, v

    def test_hash_once_per_inode_with_both_paths_preserved(self):
        u, v = self.pair()
        expected = {str(p): h.pin(p) for p in (u, v)}
        records, count = h.verify_union(expected, SimpleNamespace(verify_file=self.native_verify))
        self.assertEqual((len(records), count, len(self.calls)), (2, 1, 1))

    def test_conflicting_hardlink_hash_rejected(self):
        u, v = self.pair()
        expected = {str(p): h.pin(p) for p in (u, v)}
        expected[str(u)]['sha256'] = '0' * 64
        with self.assertRaises((c.WorkflowError, AssertionError)):
            h.verify_union(expected, SimpleNamespace(verify_file=self.native_verify))

    def test_modified_bytes_invalidate_verified_cache(self):
        u, v = self.pair()
        records, _ = h.verify_union({str(p): h.pin(p) for p in (u, v)},
                                   SimpleNamespace(verify_file=self.native_verify))
        u.write_bytes(b'x' * u.stat().st_size)
        with self.assertRaises(c.WorkflowError):
            h.check_fingerprints(records)

    def test_replaced_file_invalidates_cache(self):
        p = self.root / 'data'
        p.write_bytes(b'same bytes')
        records, _ = h.verify_union({str(p): h.pin(p)}, SimpleNamespace(verify_file=self.native_verify))
        replacement = self.root / 'new'
        replacement.write_bytes(p.read_bytes())
        os.replace(replacement, p)
        with self.assertRaises(c.WorkflowError):
            h.check_fingerprints(records)

    def test_symlink_and_new_member_rejected(self):
        p = self.root / 'data'
        p.write_bytes(b'1')
        (self.root / 'link').symlink_to(p)
        with self.assertRaises(c.WorkflowError):
            h.fingerprint(self.root / 'link')
        with self.assertRaises(c.WorkflowError):
            h.check_members({str(self.root): ['data']})

    def test_reappearing_backup_is_preserved_and_refused(self):
        u, v = self.pair()
        backup = self.root / 'backup'
        pair = {'U': str(u), 'H': str(v), 'bytes': u.stat().st_size, 'backup': str(backup)}
        h.verify_pairs([pair])
        backup.write_bytes(b'reappeared')
        with self.assertRaises(c.WorkflowError):
            h.verify_pairs([pair])
        self.assertEqual(backup.read_bytes(), b'reappeared')

    def test_identity_invokes_virtual_environment_path(self):
        executable = self.root / 'python'
        executable.symlink_to(Path(os.sys.executable).resolve())
        adapter = self.root / 'adapter.py'
        adapter.write_text('pass')
        spec = {'scientific_executable': str(executable),
                'scientific_python': h.pin(executable.resolve()), 'identity_adapter': h.pin(adapter),
                'boundary': {'path': 'boundary.json'}}
        with patch.object(h.subprocess, 'run', return_value=SimpleNamespace(returncode=0, stdout='{}', stderr='')) as run:
            h.identity(spec)
        self.assertEqual(run.call_args.args[0][0], str(executable))

    def snapshot(self):
        anchor = {'outer_pid': 1, 'name': 'supervisord', 'outer_ppid': 0, 'start_ticks': 39,
                  'namespace_pid_chain': [1], 'outer_group': 1, 'outer_session': 1}
        return {'boot_id': 'same', 'read_errors': [], 'old_pid_paths_present': [],
                'selected_mountinfo': ['94 93 0:30 / /proc rw - proc proc rw'],
                'processes': [anchor], 'self': {'outer_pid': 100}}

    def test_outer_namespace_proof_accepts_complete_view(self):
        s = self.snapshot()
        h.check_process_snapshot(s, copy.deepcopy(s), [40], {'legacy.py'})

    def test_no_processes_in_isolated_view_is_not_closure(self):
        s, old = self.snapshot(), self.snapshot()
        s['selected_mountinfo'] = ['94 93 0:99 / /proc rw - proc proc rw']
        with self.assertRaises(c.WorkflowError):
            h.check_process_snapshot(s, old, [40], {'legacy.py'})

    def test_live_legacy_child_and_visibility_error_refuse(self):
        for change in ('child', 'error', 'boot', 'anchor'):
            s, old = self.snapshot(), self.snapshot()
            if change == 'child':
                s['processes'].append({'outer_pid': 99, 'outer_group': 46579, 'outer_session': 46579})
            elif change == 'error':
                s['read_errors'] = ['could not inspect process']
            elif change == 'boot':
                s['boot_id'] = 'different'
            else:
                s['processes'][0]['start_ticks'] = 40
            with self.assertRaises(c.WorkflowError):
                h.check_process_snapshot(s, old, [40], {'legacy.py'})

    def test_controller_import_recovery_uses_real_native_verifier(self):
        source = os.environ.get('BACSANDBOX_NATIVE_VERIFIER')
        if not source:
            self.skipTest('Set BACSANDBOX_NATIVE_VERIFIER for production adapter qualification')
        native = h.load_native(h.pin(source), 'fixture_native_verifier')
        p = self.root / 'payload'
        p.write_bytes(b'disposable fixture')
        binding = h.pin(p)
        records, count = h.verify_union({str(p): binding}, native)
        self.assertEqual(count, 1)
        p.write_bytes(b'corrupted fixture!')
        with self.assertRaises(ValueError):
            native.verify_file(p, binding)
        with self.assertRaises(c.WorkflowError):
            h.check_fingerprints(records)


class BoundedImportTests(unittest.TestCase):
    def setUp(self):
        self.fixture = test_controller.WorkflowTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)

    def test_existing_receipts_import_once_then_gate(self):
        f = self.fixture
        plan = f.configure(n=2)
        # Prepare two genuine fixture side effects outside the coordinator.
        import subprocess
        for step in plan['steps']:
            subprocess.run(step['argv'], check=True)
        self.assertEqual(h.advance_imports(f.state, 1)['next_index'], 1)
        self.assertEqual(h.advance_imports(f.state, 2)['phase'], 'stopped')
        self.assertEqual((f.root / 'count0').read_text(), '1')
        self.assertEqual((f.root / 'count1').read_text(), '1')

    def test_failed_verifier_is_not_automatically_retried(self):
        f = self.fixture
        f.configure(n=1)
        with patch.object(c, 'invoke', side_effect=c.WorkflowError('failure')):
            with self.assertRaises(c.WorkflowError):
                c.advance(f.state)
        before = f.state.read_bytes()
        with self.assertRaises(c.WorkflowError):
            h.advance_imports(f.state, 2)
        self.assertEqual(f.state.read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
