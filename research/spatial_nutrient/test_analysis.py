import copy
import unittest

from research.spatial_nutrient.analysis import (
    IncompleteExperimentError, analyze_pairs, validate_records, wilson_interval,
)


def records_for(outcomes, *, seed_start=0, cohort="pilot", horizon=1.):
    records = []
    for seed, (u, h) in enumerate(outcomes, seed_start):
        for arm, survives in (("U", u), ("H", h)):
            records.append({"seed": seed, "arm": arm, "cohort": cohort,
                            "status": "completed", "validated": True,
                            "reached_time_s": horizon, "initial_focal_count": 5,
                            "focal_count": 5 if survives else 0,
                            "extinction_time_s": None if survives else horizon,
                            "diagnostics": {"births": 0, "deaths": 0 if survives else 5,
                                            "food_balance_residual_pg": 0., "exposure_contrast": .5}})
    return records


class PairedAnalysisTests(unittest.TestCase):
    def analyze(self, outcomes):
        return analyze_pairs(records_for(outcomes), range(len(outcomes)),
                             cohort="pilot", horizon_s=1., dt_s=.1,
                             bootstrap_resamples=1000)

    def test_hand_counted_table_preserves_pairing_and_zero_outcomes(self):
        result = self.analyze([(1, 1), (0, 1), (0, 1), (1, 0), (0, 0)])
        self.assertEqual(result["delta_p"], .2)
        self.assertEqual(result["contingency"], {"both_survive": 1, "H_only": 2,
                                                  "U_only": 1, "neither": 1})
        self.assertEqual(result["n_pairs"], 5)
        self.assertEqual(result["arms"]["H"]["probability"], .6)
        self.assertEqual(result["absolute_focal_counts"][-1], {"seed": 4, "U": 0, "H": 0})

    def test_degenerate_samples_use_nonzero_width_finite_sample_fallback(self):
        for outcomes, expected in [([(1, 1)] * 4, 0), ([(0, 0)] * 4, 0),
                                   ([(0, 1)] * 4, 1), ([(1, 0)] * 4, -1)]:
            with self.subTest(expected=expected, outcomes=outcomes[0]):
                result = self.analyze(outcomes)
                self.assertEqual(result["delta_p"], expected)
                self.assertTrue(result["bootstrap"]["degenerate"])
                self.assertEqual(result["reported_interval_method"], "hoeffding_degenerate_fallback")
                self.assertGreater(result["reported_interval"][1], result["reported_interval"][0])

    def test_last_tick_extinction_is_observed_survival_is_censored(self):
        result = self.analyze([(0, 1)])
        uniform, heterogeneous = result["extinction_records"]
        self.assertEqual(uniform["time_s"], 1.)
        self.assertTrue(uniform["event_observed"])
        self.assertFalse(uniform["right_censored"])
        self.assertFalse(heterogeneous["event_observed"])
        self.assertTrue(heterogeneous["right_censored"])

    def test_fixed_rng_and_input_order_replay_exactly(self):
        records = records_for([(0, 1), (1, 0), (1, 1), (0, 1)])
        kwargs = dict(registered_seeds=range(4), cohort="pilot", horizon_s=1., dt_s=.1)
        first = analyze_pairs(records, **kwargs)
        second = analyze_pairs(list(reversed(records)), **kwargs)
        self.assertEqual(first, second)

    def test_missing_duplicate_failed_and_extra_arms_are_not_dropped(self):
        original = records_for([(0, 1), (1, 0)])
        failed = copy.deepcopy(original)
        failed[-1]["status"] = "failed"
        not_valid = copy.deepcopy(original)
        not_valid[0]["validated"] = False
        extra = original + records_for([(1, 1)], seed_start=20)
        for records in [original[:-1], original + [original[0]], failed, not_valid, extra]:
            with self.subTest(records=records):
                with self.assertRaises(IncompleteExperimentError):
                    analyze_pairs(records, range(2), cohort="pilot", horizon_s=1., dt_s=.1)

    def test_mixed_cohort_partial_horizon_bad_count_and_resurrection_rejected(self):
        for field, value in [("cohort", "main"), ("reached_time_s", .9),
                             ("focal_count", -1), ("focal_count", True),
                             ("focal_count", 1), ("extinction_time_s", .25),
                             ("initial_focal_count", 0)]:
            records = records_for([(0, 1)])
            records[0][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(IncompleteExperimentError):
                validate_records(records, [0], cohort="pilot", horizon_s=1., dt_s=.1)

    def test_wilson_limits_include_boundary_and_are_nontrivial(self):
        zero = wilson_interval(0, 10)
        full = wilson_interval(10, 10)
        self.assertAlmostEqual(zero[0], 0.)
        self.assertGreater(zero[1], 0.)
        self.assertAlmostEqual(full[1], 1.)
        self.assertLess(full[0], 1.)

    def test_establishment_uses_window_minimum_not_endpoint(self):
        records = records_for([(1, 1)])
        records[0].update(focal_count=10, tail_min_focal_count=2)
        records[1].update(focal_count=10, tail_min_focal_count=5)
        result = analyze_pairs(records, [0], cohort="pilot", horizon_s=1., dt_s=.1)
        self.assertEqual(result["secondary_establishment"]["thresholds"]["5"], {"U": 0, "H": 1})
        self.assertEqual(result["secondary_establishment"]["thresholds"]["10"], {"U": 0, "H": 0})
        del records[1]["tail_min_focal_count"]
        with self.assertRaises(IncompleteExperimentError):
            analyze_pairs(records, [0], cohort="pilot", horizon_s=1., dt_s=.1)


if __name__ == "__main__":
    unittest.main()
