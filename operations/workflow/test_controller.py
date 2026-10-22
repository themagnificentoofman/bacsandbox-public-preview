"""Recovery and ownership contracts; no production commands or evidence."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('workflow_controller', Path(__file__).with_name('controller.py'))
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.plan = self.root / 'plan.json'
        self.state = self.root / 'state.json'
        self.helper = self.root / 'native.py'
        self.helper.write_text('''import json, pathlib, sys
action, receipt, counter = sys.argv[1:]
p, count = pathlib.Path(receipt), pathlib.Path(counter)
if action == 'run':
    count.write_text(str(int(count.read_text()) + 1) if count.exists() else '1')
    p.write_text(json.dumps({'status': 'complete', 'operation': p.stem}))
elif action == 'verify':
    assert count.read_text() == '1'
    assert json.loads(p.read_text()) == {'status': 'complete', 'operation': p.stem}
else:
    raise RuntimeError('bad action')
''')

    def configure(self, n=2, **updates):
        steps = []
        for i in range(n):
            receipt = self.root / f'step{i}.json'
            common = [str(receipt), str(self.root / f'count{i}')]
            steps.append(dict(id=f'step{i}', scope='retention', kind='command',
                              argv=[sys.executable, str(self.helper), 'run', *common],
                              verify_argv=[sys.executable, str(self.helper), 'verify', *common],
                              cwd=str(self.root), receipt_path=str(receipt),
                              receipt_require={'status': 'complete', 'operation': f'step{i}'},
                              capacity_path=str(self.root), minimum_free_bytes=0))
        plan = dict(schema_version=1, workflow_id='test', authorization_basis='disposable test',
                    allowed_scopes=['retention'], stop_after=f'step{n-1}', steps=steps,
                    bindings=[{'path': str(self.helper), 'sha256': c.digest(self.helper)}])
        plan.update(updates)
        self.plan.write_text(json.dumps(plan))
        c.initialize(self.plan, self.state)
        return plan

    def test_batch_boundary_and_no_duplicate_effect(self):
        self.configure()
        self.assertEqual(c.advance(self.state, 1)['next_index'], 1)
        self.assertFalse((self.root / 'count1').exists())
        self.assertEqual(c.advance(self.state)['phase'], 'stopped')
        c.advance(self.state)
        self.assertEqual((self.root / 'count0').read_text(), '1')

    def test_stop_before_later_step(self):
        self.configure(stop_after='step0')
        self.assertEqual(c.advance(self.state)['phase'], 'stopped')
        self.assertFalse((self.root / 'count1').exists())

    def test_crash_after_side_effect_reconciles_without_rerun(self):
        self.configure(n=1)
        with patch.object(c, 'verify', side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt):
                c.advance(self.state)
        self.assertEqual(c.status(self.state)['phase'], 'in_flight')
        with self.assertRaises(c.WorkflowError):
            c.advance(self.state)
        self.assertEqual(c.reconcile(self.state)['phase'], 'stopped')
        self.assertEqual((self.root / 'count0').read_text(), '1')

    def test_existing_receipt_is_verified_not_reexecuted(self):
        p = self.configure(n=1)
        subprocess.run(p['steps'][0]['argv'], check=True)
        self.assertEqual(c.advance(self.state)['phase'], 'in_flight')
        self.assertEqual(c.reconcile(self.state)['phase'], 'stopped')
        self.assertEqual((self.root / 'count0').read_text(), '1')

    def test_receipt_without_actual_effect_cannot_reconcile(self):
        self.configure(n=1)
        (self.root / 'step0.json').write_text('{"status":"complete","operation":"step0"}')
        c.advance(self.state)
        with self.assertRaises(c.WorkflowError):
            c.reconcile(self.state)
        self.assertEqual(c.status(self.state)['next_index'], 0)

    def test_binding_change_blocks_before_mutation(self):
        self.configure()
        self.helper.write_text(self.helper.read_text() + '\n# changed\n')
        with self.assertRaises(c.WorkflowError):
            c.advance(self.state)
        self.assertFalse((self.root / 'count0').exists())

    def test_changed_completed_receipt_blocks_next_action(self):
        self.configure()
        c.advance(self.state, 1)
        (self.root / 'step0.json').write_text('{}')
        with self.assertRaises(c.WorkflowError):
            c.advance(self.state)
        self.assertFalse((self.root / 'count1').exists())

    def test_plan_edit_is_not_an_implicit_authorization(self):
        self.configure()
        self.plan.write_text(self.plan.read_text() + ' ')
        with self.assertRaises(c.WorkflowError):
            c.advance(self.state)

    def test_capacity_hold_has_no_side_effect(self):
        self.configure()
        with patch.object(c.shutil, 'disk_usage', return_value=type('D', (), {'free': -1})()):
            self.assertEqual(c.advance(self.state)['next_index'], 0)
        self.assertFalse((self.root / 'count0').exists())
        self.assertEqual(c.advance(self.state)['phase'], 'stopped')

    def test_activation_block_and_status_read_only(self):
        self.configure(activation_blocker='finish legacy transaction first')
        self.assertEqual(c.advance(self.state)['phase'], 'blocked')
        before = self.state.read_bytes()
        c.status(self.state)
        self.assertEqual(before, self.state.read_bytes())
        self.assertFalse((self.root / 'count0').exists())

    def test_single_owner(self):
        self.configure()
        with c.owner(self.state):
            p = subprocess.run([sys.executable, str(Path(c.__file__)), 'advance', '--state', str(self.state)],
                               capture_output=True, text=True)
        self.assertEqual(p.returncode, 2)
        self.assertIn('holds the lock', p.stderr)

    def test_child_retains_lock_when_parent_descriptor_closes(self):
        self.configure()
        read_fd, write_fd = os.pipe()
        child = None
        try:
            with c.owner(self.state) as lock_fd:
                child = subprocess.Popen([sys.executable, '-c', 'import os,sys; os.read(int(sys.argv[1]),1)', str(read_fd)],
                                         pass_fds=(lock_fd, read_fd))
            with self.assertRaises(c.WorkflowError):
                with c.owner(self.state):
                    pass
        finally:
            os.write(write_fd, b'x')
            os.close(write_fd)
            os.close(read_fd)
            if child:
                child.wait(timeout=10)
        with c.owner(self.state):
            pass

    def test_failed_action_has_no_automatic_retry(self):
        self.configure()
        with patch.object(c, 'invoke', side_effect=c.WorkflowError('native failure')):
            with self.assertRaises(c.WorkflowError):
                c.advance(self.state)
        self.assertEqual(c.status(self.state)['phase'], 'failed')
        with self.assertRaises(c.WorkflowError):
            c.advance(self.state)
        self.assertFalse((self.root / 'count0').exists())

    def test_state_cannot_be_silently_reinitialized(self):
        self.configure()
        with self.assertRaises(c.WorkflowError):
            c.initialize(self.plan, self.state)


if __name__ == '__main__':
    unittest.main()
