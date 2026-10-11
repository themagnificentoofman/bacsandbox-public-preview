"""Analytical and numerical invariants for the scaled research reference."""

from dataclasses import replace
import math
import unittest
from unittest.mock import patch

import numpy as np

from research.zhou_reference.model import (
    ModelNumericalError, ModelParameters, StageDomainError, check_convergence, simulate,
)
from research.zhou_reference import model


class ReferenceModelTests(unittest.TestCase):
    def setUp(self):
        self.initial = np.linspace(0.01, 0.18, 18)
        self.times = [0.0, 1.0, 10.0, 28.0]
        self.p = ModelParameters()

    def test_pure_shared_loss_matches_analytical_solution(self):
        p = replace(self.p, qmax_R_target=0, qmax_R_background=0,
                    qmax_G_target=0, qmax_G_background=0, loss=0.3)
        result = simulate(self.initial, 5, self.times, p, True)
        expected = self.initial[None, :] * np.exp(-0.3 * np.array(self.times)[:, None])
        np.testing.assert_allclose(result.biomass, expected, rtol=2e-13, atol=0)
        np.testing.assert_allclose(result.target_log_change,
                                   -0.3 * np.array(self.times) / math.log(10), atol=2e-14)

    def test_resource_input_and_washout_match_analytical_solution(self):
        p = replace(self.p, qmax_R_target=0, qmax_R_background=0,
                    qmax_G_target=0, qmax_G_background=0,
                    initial_R=0.3, initial_G=0.2, input_R=0.8, input_G=0.4,
                    washout_R=0.2, washout_G=0.7)
        result = simulate(self.initial, 0, self.times, p, True)
        for col, resource, flow, decay in [(0, 0.3, 0.8, 0.2), (1, 0.2, 0.4, 0.7)]:
            expected = flow / decay + (resource - flow / decay) * np.exp(-decay * np.array(self.times))
            np.testing.assert_allclose(result.resources[:, col], expected, rtol=1e-9, atol=1e-11)

    def test_closed_unit_yield_system_conserves_resource_plus_biomass(self):
        p = replace(self.p, loss=0, input_R=0, input_G=0,
                    initial_R=2, initial_G=1, washout_R=0, washout_G=0)
        result = simulate(self.initial, 7, [0, 0.25, 1, 2], p, True)
        totals = result.biomass.sum(axis=1) + result.resources.sum(axis=1)
        np.testing.assert_allclose(totals, self.initial.sum() + 3, rtol=1e-8, atol=1e-10)
        np.testing.assert_allclose(result.accounting["growth"].sum(axis=1),
                                   result.accounting["uptake"].sum(axis=1), atol=1e-12)

    def test_resource_and_biomass_accounts_close_with_unequal_yields(self):
        p = replace(self.p, Y_R_target=0.4, Y_G_target=0.7,
                    Y_R_background=0.6, Y_G_background=0.2)
        result = simulate(self.initial, 4, self.times, p, True)
        expected_resource = (np.array([p.initial_R, p.initial_G])
                             + result.accounting["input"] - result.accounting["outflow"]
                             - result.accounting["uptake"])
        np.testing.assert_allclose(result.resources, expected_resource, atol=1e-12)
        expected_total_biomass = (self.initial.sum() + result.accounting["growth"].sum(axis=1)
                                  - result.accounting["loss"].sum(axis=1))
        np.testing.assert_allclose(result.biomass.sum(axis=1), expected_total_biomass, atol=2e-8)

    def test_identical_inputs_are_exact_under_both_diet_labels(self):
        p = replace(self.p, input_G=0, initial_G=0.4)
        control = simulate(self.initial, 2, self.times, p, False)
        treated = simulate(self.initial, 2, self.times, p, True)
        np.testing.assert_array_equal(control.biomass, treated.biomass)
        np.testing.assert_array_equal(control.resources, treated.resources)

    def test_no_glucomannan_ignores_glucomannan_capability(self):
        control = simulate(self.initial, 2, self.times, self.p, False)
        changed = simulate(self.initial, 2, self.times,
                           replace(self.p, qmax_G_target=4, qmax_G_background=8), False)
        np.testing.assert_array_equal(control.target_log_change, changed.target_log_change)
        np.testing.assert_array_equal(control.resources[:, 1], 0)

    def test_tied_background_ratios_stay_constant(self):
        result = simulate(self.initial, 3, self.times, self.p, True)
        background = np.arange(18) != 3
        ratios = result.biomass[:, background] / result.biomass[:, [0]]
        expected = self.initial[background] / self.initial[0]
        np.testing.assert_allclose(ratios, np.broadcast_to(expected, ratios.shape), rtol=3e-14)

    def test_zero_labels_and_rare_positive_labels_are_preserved(self):
        initial = self.initial.copy()
        initial[1] = 0
        initial[2] = 1e-250
        result = simulate(initial, 0, self.times, self.p, True)
        np.testing.assert_array_equal(result.biomass[:, 1], 0)
        self.assertTrue(np.all(result.biomass[:, 2] > 0))
        ratio = result.biomass[:, 2] / result.biomass[:, 3]
        np.testing.assert_allclose(ratio, initial[2] / initial[3], rtol=8e-14, atol=0)

    def test_empty_background_remains_empty(self):
        initial = np.zeros(18)
        initial[9] = 0.3
        result = simulate(initial, 9, self.times, self.p, True)
        np.testing.assert_array_equal(result.biomass[:, np.arange(18) != 9], 0)
        self.assertTrue(np.all(result.biomass[:, 9] > 0))

    def test_baseline_and_optional_times_do_not_change_integration_origin(self):
        full = simulate(self.initial, 0, self.times, self.p, True)
        later = simulate(self.initial, 0, self.times[1:], self.p, True)
        np.testing.assert_array_equal(full.biomass[1:], later.biomass)
        only_zero = simulate(self.initial, 0, [0], self.p, True)
        np.testing.assert_array_equal(only_zero.biomass[0], self.initial)
        np.testing.assert_array_equal(only_zero.target_log_change, 0)

    def test_delayed_input_start_is_an_exact_boundary(self):
        p = replace(self.p, qmax_R_target=0, qmax_R_background=0,
                    qmax_G_target=0, qmax_G_background=0,
                    input_G=0.4, washout_G=0.7, glucomannan_start_day=2)
        times = np.array([0, 1, 2, 2.5, 4])
        result = simulate(self.initial, 0, times, p, True)
        expected = 0.4 / 0.7 * (1 - np.exp(-0.7 * np.maximum(0, times - 2)))
        np.testing.assert_allclose(result.resources[:, 1], expected, rtol=1e-9, atol=1e-11)
        np.testing.assert_allclose(result.accounting["input"][:, 1],
                                   0.4 * np.maximum(0, times - 2), atol=1e-14)

    def test_prespecified_tighter_solver_converges(self):
        standard = simulate(self.initial, 5, self.times, self.p, True)
        tighter = simulate(self.initial, 5, self.times, self.p, True, tighter=True)
        self.assertTrue(check_convergence(standard, tighter)["passed"])
        np.testing.assert_allclose(standard.resources, tighter.resources, atol=2e-8)

    def test_positive_biomass_underflow_is_failure_not_extinction(self):
        p = replace(self.p, qmax_R_target=0, qmax_R_background=0,
                    qmax_G_target=0, qmax_G_background=0, loss=100)
        with self.assertRaisesRegex(ModelNumericalError, "underflow"):
            simulate(self.initial, 0, [0, 30], p, False)

    def test_historical_control_stage_fixture_converges_with_valid_attempts(self):
        # Published baseline fixture; the failure was found in control-only
        # calibration. No observed follow-up abundance enters this regression.
        initial = [23.224617781999974, 0.40414732296326944, 1.579523102444941,
                   4.733712236125903, 0.21973815586002415, 5.3743847457332725,
                   0.2576906721311042, 4.158946992046722, 0.0004963430432462628,
                   0.00031907767065831184, 0.19835995192591682, 0.002162637545573002,
                   0.0015599352787739695, 0.11245715236979617, 0.23289124650605042,
                   3.5090035034530036, 4.464853845521766, 0.1]
        p = replace(self.p, loss=0.022505466356404744)
        standard = simulate(initial, 0, [10, 28], p, False)
        tight = simulate(initial, 0, [10, 28], p, False, tighter=True)
        self.assertTrue(check_convergence(standard, tight)["passed"])
        # Adaptive trial steps can differ across supported numerical runtimes.
        # Require convergence and valid diagnostics, not a particular overshoot.
        for result, limits in [(standard, [1.0, 0.5, 0.25]),
                               (tight, [0.5, 0.25, 0.125])]:
            with self.subTest(limits=limits):
                attempts = result.diagnostics["integration_attempts"]
                self.assertGreaterEqual(len(attempts), 1)
                self.assertLessEqual(len(attempts), len(limits))
                self.assertEqual([a["max_step_days"] for a in attempts],
                                 limits[:len(attempts)])
                self.assertEqual([a["success"] for a in attempts],
                                 [False] * (len(attempts) - 1) + [True])
                for failed in attempts[:-1]:
                    self.assertEqual(failed["failure_type"], "StageDomainError")
                    self.assertLessEqual(min(failed["stage"]["denominators"]), 0)
                self.assertEqual(result.diagnostics["nfev"],
                                 sum(a["rhs_calls"] for a in attempts))

    def test_injected_stage_failure_recovers_on_identical_segment(self):
        actual_solver = model.solve_ivp
        calls = []

        def fail_first(fun, interval, initial, **kwargs):
            calls.append((interval, initial.copy()))
            if len(calls) == 1:
                bad = initial.copy()
                bad[2] = -2
                fun(interval[0], bad)
            return actual_solver(fun, interval, initial, **kwargs)

        with patch.object(model, "solve_ivp", side_effect=fail_first):
            result = simulate(self.initial, 0, [0, 1], self.p, False)
        attempts = result.diagnostics["integration_attempts"]
        self.assertEqual([a["max_step_days"] for a in attempts], [1.0, 0.5])
        self.assertEqual([a["success"] for a in attempts], [False, True])
        self.assertEqual(attempts[0]["failure_type"], "StageDomainError")
        self.assertEqual(calls[0][0], calls[1][0])
        np.testing.assert_array_equal(calls[0][1], calls[1][1])
        self.assertEqual(result.diagnostics["nfev"], sum(a["rhs_calls"] for a in attempts))
        tight = simulate(self.initial, 0, [0, 1], self.p, False, tighter=True)
        self.assertTrue(check_convergence(result, tight)["passed"])

    def test_stage_retries_are_bounded_and_restart_identical_segment_state(self):
        bad_stage = np.zeros(14)
        bad_stage[2] = -2
        for tighter, limits in [(False, [1, 0.5, 0.25]), (True, [0.5, 0.25, 0.125])]:
            with self.subTest(tighter=tighter):
                with patch.object(model, "solve_ivp", side_effect=lambda *a, **k:
                                  (_ for _ in ()).throw(StageDomainError(0.5, bad_stage, np.ones(2)))) as solver:
                    with self.assertRaisesRegex(ModelNumericalError, "retries exhausted") as caught:
                        simulate(self.initial, 0, [0, 1], self.p, False, tighter=tighter)
                self.assertEqual(solver.call_count, 3)
                self.assertIsInstance(caught.exception.__cause__, StageDomainError)
                attempts = caught.exception.integration_attempts
                self.assertEqual([a["max_step_days"] for a in attempts], limits)
                self.assertTrue(all(not a["success"] and "stage" in a for a in attempts))
                for call in solver.call_args_list:
                    np.testing.assert_array_equal(call.args[2], solver.call_args_list[0].args[2])
                    self.assertEqual(call.kwargs["method"], "DOP853")
                    self.assertEqual(call.kwargs["rtol"], 1e-10 if tighter else 1e-8)

    def test_nonfinite_stage_is_not_retried(self):
        def nonfinite_trial(fun, interval, initial, **kwargs):
            bad = initial.copy()
            bad[2] = float("nan")
            fun(interval[0], bad)
        with patch.object(model, "solve_ivp", side_effect=nonfinite_trial) as solver:
            with self.assertRaisesRegex(ModelNumericalError, "Nonfinite solver stage") as caught:
                simulate(self.initial, 0, [0, 1], self.p, False)
        self.assertEqual(solver.call_count, 1)
        self.assertEqual(len(caught.exception.integration_attempts), 1)
        self.assertEqual(caught.exception.integration_attempts[0]["rhs_calls"], 1)

    def test_accepted_state_and_accounting_failures_are_not_retried(self):
        actual_solver = model.solve_ivp
        for index, value, message in [(2, -1, "Negative accepted"), (8, 10, "Accounting residual")]:
            def corrupted_accepted(*args, **kwargs):
                solution = actual_solver(*args, **kwargs)
                solution.y[index, -1] = value
                return solution
            with self.subTest(index=index):
                with patch.object(model, "solve_ivp", side_effect=corrupted_accepted) as solver:
                    with self.assertRaisesRegex(ModelNumericalError, message) as caught:
                        simulate(self.initial, 0, [0, 1], self.p, False)
                self.assertEqual(solver.call_count, 1)
                self.assertEqual(len(caught.exception.integration_attempts), 1)

    def test_invalid_inputs_fail_before_integration(self):
        for name, value in [("K_R", 0), ("input_G", -1), ("loss", float("nan")),
                            ("qmax_G_target", float("inf")), ("input_R", True)]:
            with self.subTest(name=name), self.assertRaises(ValueError):
                replace(self.p, **{name: value})
        for times in [[], [1, 0], [0, 0], [-1], [float("inf")]]:
            with self.subTest(times=times), self.assertRaises(ValueError):
                simulate(self.initial, 0, times, self.p, False)
        for initial in [np.ones(17), -self.initial, np.full(18, np.nan), np.zeros(18)]:
            with self.subTest(initial=str(initial)), self.assertRaises(ValueError):
                simulate(initial, 0, [0, 1], self.p, False)
        with self.assertRaises(ValueError):
            simulate(self.initial, True, [0, 1], self.p, False)
        with self.assertRaises(ValueError):
            simulate(self.initial, 0, [0, 1], self.p, "treated")


if __name__ == "__main__":
    unittest.main()
