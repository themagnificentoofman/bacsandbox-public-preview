import copy
import unittest

from research.spatial_nutrient.qualification import (
    REQUIRED_CHECKS, bounded_discrete_mean_interval, canonical_sha256,
    compare_timesteps, qualification_gate,
)
from research.spatial_nutrient.test_analysis import records_for


class NumericalQualificationTests(unittest.TestCase):
    def compare(self, coarse, fine, n):
        return compare_timesteps(records_for(coarse, seed_start=2000, cohort="numerical"),
                                 records_for(fine, seed_start=2000, cohort="numerical"),
                                 range(2000, 2000 + n), horizon_s=1.,
                                 bootstrap_resamples=200)

    def test_zero_changes_at_50_is_inconclusive_not_converged(self):
        result = self.compare([(1, 1)] * 50, [(1, 1)] * 50, 50)
        self.assertEqual(result["Q_mean"], 0.)
        self.assertEqual(result["status"], "inconclusive")
        self.assertGreater(result["Q_interval"][1], .05)
        self.assertEqual(result["next_cumulative_pair_count"], 200)
        self.assertEqual(result["next_registered_seeds"], list(range(2000, 2200)))

    def test_sufficient_registered_evidence_can_pass(self):
        result = self.compare([(1, 1)] * 1000, [(1, 1)] * 1000, 1000)
        self.assertEqual(result["status"], "passed")
        self.assertLess(result["Q_interval"][1], .05)
        self.assertTrue(result["registered_look"])

    def test_common_arm_error_cannot_hide_behind_zero_difference(self):
        result = self.compare([(0, 0)] * 1000, [(1, 1)] * 1000, 1000)
        self.assertEqual(result["Q_mean"], 0.)
        self.assertEqual(result["Q_status"], "passed")
        self.assertEqual(result["arms"]["U"]["status"], "failed")
        self.assertEqual(result["status"], "failed")

    def test_effect_reversal_is_visible(self):
        result = self.compare([(1, 0)] * 50, [(0, 1)] * 50, 50)
        self.assertEqual(result["Q_mean"], 2.)
        self.assertEqual(result["status"], "failed")

    def test_unknown_sample_size_cannot_qualify(self):
        result = self.compare([(1, 1)] * 500, [(1, 1)] * 500, 500)
        self.assertFalse(result["registered_look"])
        self.assertEqual(result["status"], "inconclusive")

    def test_discrete_interval_remembers_unseen_categories(self):
        lower, upper = bounded_discrete_mean_interval([0] * 50, [-2, -1, 0, 1, 2])
        self.assertLess(lower, 0)
        self.assertGreater(upper, 0)
        with self.assertRaises(ValueError):
            bounded_discrete_mean_interval([.3], [-1, 0, 1])

    def test_required_diagnostics_are_not_silently_omitted(self):
        records = records_for([(1, 1)] * 50, seed_start=2000, cohort="numerical")
        changed = copy.deepcopy(records)
        del changed[0]["diagnostics"]["births"]
        with self.assertRaisesRegex(ValueError, "births"):
            compare_timesteps(records, changed, range(2000, 2050), horizon_s=1.)


class FreezeGateTests(unittest.TestCase):
    def evidence(self):
        checks = {name: {"status": "passed", "evidence_sha256": "a" * 64} for name in REQUIRED_CHECKS}
        checks["pilot"].update(seeds=list(range(20)), completed_pairs=20)
        checks["numerical"].update(n_pairs=1000, registered_seeds=list(range(2000, 3000)),
                                    registered_look=True, Q_status="passed",
                                    arms={arm: {"status": "passed"} for arm in ("U", "H")})
        criterion = {"metric": "mean_paired_time_integrated_half_food_contrast", "threshold": .05,
                     "operator": ">="}
        checks["exposure"].update(criterion=criterion, criterion_sha256=canonical_sha256(criterion),
                                  measured=.3)
        return {"protocol_sha256": "b" * 64, "code_sha256": "c" * 64, "checks": checks}

    def gate(self, evidence):
        return qualification_gate(evidence, protocol_sha256="b" * 64, code_sha256="c" * 64)

    def test_complete_version_matched_evidence_allows_freeze(self):
        self.assertTrue(self.gate(self.evidence())["freeze_allowed"])

    def test_missing_evidence_or_uncertainty_blocks_freeze(self):
        evidence = self.evidence()
        del evidence["checks"]["lifecycle"]
        result = self.gate(evidence)
        self.assertEqual(result["status"], "inconclusive")
        self.assertFalse(result["freeze_allowed"])
        evidence = self.evidence()
        evidence["checks"]["numerical"]["status"] = "inconclusive"
        self.assertFalse(self.gate(evidence)["freeze_allowed"])

    def test_changed_code_or_protocol_invalidates_freeze(self):
        for field in ["protocol_sha256", "code_sha256"]:
            evidence = self.evidence()
            evidence[field] = "d" * 64
            self.assertEqual(self.gate(evidence)["status"], "failed")

    def test_incomplete_pilot_and_outcome_tuned_exposure_are_rejected(self):
        evidence = self.evidence()
        evidence["checks"]["pilot"]["completed_pairs"] = 19
        self.assertFalse(self.gate(evidence)["freeze_allowed"])
        evidence = self.evidence()
        criterion = evidence["checks"]["exposure"]["criterion"]
        criterion["metric"] = "delta_p"
        evidence["checks"]["exposure"]["criterion_sha256"] = canonical_sha256(criterion)
        self.assertFalse(self.gate(evidence)["freeze_allowed"])


if __name__ == "__main__":
    unittest.main()
