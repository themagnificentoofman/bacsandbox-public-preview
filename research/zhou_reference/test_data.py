"""Source roles, counterfactual initialization and endpoint identity contracts."""
from copy import deepcopy
import unittest
from unittest.mock import patch

import numpy as np

from research.zhou_reference import data


class SourceRoleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = data.load_baselines()
        cls.calibration = data.load_calibration()
        cls.evaluation = data.load_evaluation()

    def test_roles_counts_identity_and_endpoint_mask(self):
        summary = data.role_manifest(self.inputs, self.calibration, self.evaluation)
        self.assertEqual(summary["endpoint_counts"], dict(zip(data.PHASES[1:], (39, 38))))
        self.assertEqual(len(self.inputs["animals"]), 39)
        for animal in self.inputs["animals"]:
            self.assertNotIn("followups", animal)
            self.assertNotIn("raw_log10", animal)
        self.assertEqual(len(self.calibration), 20)
        self.assertEqual(len(self.evaluation), 57)
        treated14 = [r for r in self.evaluation if r["experiment"] == "Iteration 2"
                     and r["diet"] == "Mirpur-18 + Glucomannan" and r["sex"] == "Male"
                     and r["mouse_id"] == "14"]
        control14 = [r for r in self.evaluation if r["experiment"] == "Iteration 2"
                     and r["diet"] == "Mirpur-18" and r["sex"] == "Male"
                     and r["mouse_id"] == "14"]
        self.assertEqual(len(treated14), 1)
        self.assertEqual(len(control14), 2)
        self.assertNotEqual(treated14[0]["animal_id"], control14[0]["animal_id"])

    def test_evaluation_rows_cannot_be_mislabeled_as_calibration(self):
        for evaluation_row in self.evaluation:
            bad = dict(evaluation_row, role="calibration")
            rows = list(self.calibration)
            rows[0] = bad
            with self.assertRaises(ValueError):
                data.validate_calibration(rows, self.inputs)
        duplicate = list(self.calibration)
        duplicate[-1] = duplicate[0]
        with self.assertRaises(ValueError):
            data.validate_calibration(duplicate, self.inputs)

    def test_control_loader_never_exports_evaluation_followup_values(self):
        mapped = data._mapping(data.DEFAULT_SOURCE)
        for row in mapped["target_observations.csv"]:
            if row["phase_label"] != data.PHASES[0] and not data._is_calibration(row):
                # This object cannot be serialized or interpreted as a measurement.
                row["raw_log10"] = object()
        with patch.object(data, "_mapping", return_value=mapped):
            actual = data.load_calibration()
        self.assertEqual(actual, self.calibration)

    def test_explicit_zero_scenarios_preserve_source_values(self):
        animal = next(a for a in self.inputs["animals"] if any(a["zero_mask"]))
        before = data.fingerprint(animal)
        arrays = [data.initialize(animal, {"zero_policy": z, "biomass_scale": 1,
                                          "abundance_scale": 1e10})
                  for z in ("absent", "low", "high")]
        mask = np.array(animal["zero_mask"])
        for array, value in zip(arrays, (0, 1e-4, 1e-2)):
            np.testing.assert_array_equal(array[mask], value)
        np.testing.assert_array_equal(arrays[0][~mask], arrays[2][~mask])
        self.assertEqual(before, data.fingerprint(animal))
        doubled = data.initialize(animal, {"zero_policy": "absent", "biomass_scale": 2})
        np.testing.assert_array_equal(doubled, 2 * arrays[0])

    def test_scales_and_zero_policy_are_required_and_validated(self):
        animal = self.inputs["animals"][0]
        for bad in (-1, 0, float("nan"), float("inf"), True):
            with self.assertRaises(ValueError):
                data.initialize(animal, {"zero_policy": "absent", "biomass_scale": bad})
        with self.assertRaises((KeyError, ValueError)):
            data.initialize(animal, {})
        with self.assertRaises(ValueError):
            data.initialize(animal, {"zero_policy": "assume_detection_limit"})

    def test_chronology_is_explicit_and_keeps_iteration_difference(self):
        self.assertEqual(data.get_times("Iteration 1", "age"), (11, 29))
        self.assertEqual(data.get_times("Iteration 2", "age"), (12, 30))
        self.assertEqual(data.get_times("Iteration 1", "label"), (10, 28))
        with self.assertRaises(ValueError):
            data.get_times("Iteration 1", "choose_best_fit")
        with self.assertRaises(ValueError):
            data.phase_index("Weaning")


if __name__ == "__main__":
    unittest.main()
