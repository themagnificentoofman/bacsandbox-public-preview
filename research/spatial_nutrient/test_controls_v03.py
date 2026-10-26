"""Production-path controls, retained evidence integrity, and exact denominator."""
import copy
import gzip
import json
from pathlib import Path
import tempfile
import unittest

from . import controls_v03 as controls
from .inputs import scoped_python_random
from .runner import atomic_json, digest


class IsolatedControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan = json.loads(controls.PLAN_PATH.read_text())
        cls.base = json.loads(controls.BASE_PATH.read_text())

    def fixture(self, condition, role, lifespan=.4):
        spec = next(s for s in controls.fixture_specs() if s['condition'] == condition and
                    s['role'] == role and s['lifespan_multiplier'] == lifespan)
        return controls._make_input(self.plan, self.base, spec)

    def test_prespecified_denominator_and_no_inference_seeds(self):
        specs = controls.fixture_specs()
        self.assertEqual(len(specs), 240)
        self.assertEqual(len({s['fixture_id'] for s in specs}), 240)
        self.assertEqual({s['seed'] for s in specs}, set(range(6000, 6020)))

    def test_all_lifespan_extremes_starve_without_feed_with_closed_ledgers(self):
        for role in controls.ROLES:
            for lifespan in controls.LIFESPANS:
                with self.subTest(role=role, lifespan=lifespan):
                    inputs = self.fixture('no_feed', role, lifespan)
                    events = []
                    with scoped_python_random(inputs['streams']['ecology']['python_state']):
                        sim = controls.build_control(inputs, events.append)
                        for _ in range(1000):
                            sim.update(800, 600, .1)
                            if not sim.species[0].members:
                                break
                    self.assertFalse(sim.species[0].members)
                    self.assertEqual(sim.experiment.cumulative['births'], 0)
                    self.assertEqual(sim.experiment.cumulative['feed_pg'], 0)
                    self.assertEqual([e['cause'] for e in events if e['kind'] == 'death'], ['starvation'])
                    self.assertAlmostEqual(sim.experiment.cumulative['maintenance_pg'], 2.5)
                    self.assertAlmostEqual(sim.experiment.last_tick['mass_residual_pg'], 0)

    def test_actual_rates_stochastic_division_and_accounted_node_replenishment(self):
        for role, uptake in [('focal', .5), ('resident', 1)]:
            with self.subTest(role=role):
                inputs = self.fixture('abundant_local_feed', role)
                with scoped_python_random(inputs['streams']['ecology']['python_state']):
                    sim = controls.build_control(inputs)
                    nodes = {}
                    self.assertEqual(sim.species[0].consumption_speed, uptake)
                    first_node = None
                    for _ in range(6000):
                        controls.replenish_abundance(sim, nodes)
                        if first_node is None:
                            first_node = nodes[0]
                        self.assertIs(nodes[0], first_node)
                        self.assertEqual(nodes[0].amount, 100)
                        sim.update(800, 600, .1)
                        if sim.experiment.cumulative['births']:
                            # Exercise the newly born cell's own node next tick.
                            controls.replenish_abundance(sim, nodes)
                            sim.update(800, 600, .1)
                            break
                self.assertEqual(sim.experiment.cumulative['births'], 1)
                self.assertEqual(sim.experiment.cumulative['deaths'], 0)
                self.assertEqual(len(nodes), 2)
                self.assertEqual(len(sim.food_nodes), 2)
                self.assertGreater(sim.experiment.cumulative['feed_pg'], 200)
                self.assertLess(sim.step_idx * .1, 600)
                self.assertAlmostEqual(sim.experiment.last_tick['mass_residual_pg'], 0)
                for cell in sim.species[0].members:
                    self.assertEqual(cell.base_size, 1)
                    self.assertEqual(cell.speed, 0)
                    self.assertEqual((cell.x, cell.y), (0, 0))

    def test_replenishment_never_hides_motion_or_overfull_node(self):
        inputs = self.fixture('abundant_local_feed', 'focal')
        with scoped_python_random(inputs['streams']['ecology']['python_state']):
            sim = controls.build_control(inputs)
            nodes = {}
            controls.replenish_abundance(sim, nodes)
            nodes[0].amount += 1
            with self.assertRaisesRegex(ValueError, '100 pg'):
                controls.replenish_abundance(sim, nodes)
            nodes[0].amount = 100
            sim.species[0].members[0].x += 1
            with self.assertRaisesRegex(ValueError, 'moved'):
                controls.replenish_abundance(sim, nodes)

    def test_denominator_and_input_tampering_fail_closed(self):
        with tempfile.TemporaryDirectory() as parent:
            output = Path(parent) / 'controls'
            manifest = controls.register_controls(output)
            row = manifest['fixtures'][0]
            path = output / 'inputs' / (row['fixture_id'] + '.json')
            inputs = json.loads(path.read_text())
            inputs['founder_state']['food_level'] = 100
            atomic_json(path, inputs)
            with self.assertRaisesRegex(ValueError, 'input hash'):
                controls._read_input(output, manifest, row)
            # Rehashing an altered founder must not defeat the fixed input rule.
            changed_row = {**row, 'input_sha256': digest(inputs)}
            with self.assertRaisesRegex(ValueError, 'founder or random streams'):
                controls._read_input(output, manifest, changed_row)
            shortened = copy.deepcopy(manifest)
            shortened['fixtures'].pop()
            shortened.pop('manifest_sha256')
            shortened['manifest_sha256'] = digest(shortened)
            atomic_json(output / 'manifest.json', shortened)
            with self.assertRaisesRegex(ValueError, 'entire240'):
                controls._read_manifest(output, check_runtime=False)

    def test_full_horizon_retained_source_hashes_detect_changed_events(self):
        with tempfile.TemporaryDirectory() as parent:
            output = Path(parent) / 'controls'
            manifest = controls.register_controls(output)
            row = manifest['fixtures'][0]
            result = controls._run_fixture(output, row['fixture_id'])
            self.assertTrue(result['passed'])
            self.assertEqual(result['completed_ticks'], 6000)
            self.assertEqual(controls._read_result(output, manifest, row), result)
            changed = copy.deepcopy(result)
            changed['first_starvation_s'] += .1
            changed.pop('result_sha256')
            changed['result_sha256'] = digest(changed)
            result_path = output / 'fixtures' / row['fixture_id'] / 'result.json'
            atomic_json(result_path, changed)
            with self.assertRaisesRegex(ValueError, 'summary differs from retained sources'):
                controls._read_result(output, manifest, row)
            atomic_json(result_path, result)
            path = output / 'fixtures' / row['fixture_id'] / 'events.jsonl.gz'
            with path.open('ab') as handle:
                handle.write(b'altered evidence')
            with self.assertRaisesRegex(ValueError, 'artifact hash mismatch'):
                controls._read_result(output, manifest, row)

    def test_failed_attempt_keeps_readable_complete_compressed_evidence(self):
        with tempfile.TemporaryDirectory() as parent:
            path = Path(parent) / 'partial.jsonl.gz'
            with self.assertRaisesRegex(RuntimeError, 'deliberate fixture failure'):
                with controls._gzip_jsonlines(path) as write:
                    write({'tick': 0, 'kind': 'maintenance', 'amount_pg': .01})
                    raise RuntimeError('deliberate fixture failure')
            with gzip.open(path, 'rt') as handle:
                rows = [json.loads(row) for row in handle]
            self.assertEqual(rows, [{'tick': 0, 'kind': 'maintenance', 'amount_pg': .01}])


if __name__ == '__main__':
    unittest.main()
