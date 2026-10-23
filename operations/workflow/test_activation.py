"""A shared operational owner cannot be silently replaced after interruption."""
import json
from pathlib import Path
import tempfile
import unittest

import activate
import controller as c


class ActivationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.registry = Path(self.tmp.name) / 'OWNER.json'
        self.writes = 0

    def publish(self, path, value):
        with Path(path).open('x') as f:
            json.dump(value, f)
        self.writes += 1

    def test_claim_recovery_reuses_identical_owner_without_second_write(self):
        value = {'state_path': '/one/state', 'admission': 'exact admission'}
        a = activate.claim(self.registry, value, self.publish)
        b = activate.claim(self.registry, value, self.publish)
        self.assertEqual(a, b)
        self.assertEqual(self.writes, 1)

    def test_different_state_cannot_claim_same_workflow(self):
        activate.claim(self.registry, {'state_path': '/one'}, self.publish)
        with self.assertRaises(c.WorkflowError):
            activate.claim(self.registry, {'state_path': '/two'}, self.publish)
        self.assertEqual(c.read(self.registry), {'state_path': '/one'})

    def test_changed_implementation_requires_explicit_migration(self):
        activate.claim(self.registry, {'commit': 'first'}, self.publish)
        with self.assertRaises(c.WorkflowError):
            activate.claim(self.registry, {'commit': 'second'}, self.publish)
        self.assertEqual(self.writes, 1)
