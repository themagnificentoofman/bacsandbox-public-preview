"""Meaningful numerical-gate, leakage, mask and immutable-freeze regression tests."""
from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from research.zhou_reference import runner


SETTINGS = {"lower":0.0,"upper":2.0,"starts":[0.2,1.0,1.8],"method":"trf",
    "loss":"linear","jac":"3-point","tr_solver":"exact","x_scale":1.0,
    "ftol":1e-8,"xtol":1e-8,"gtol":1e-8,"max_nfev":2000,
    "required_successful_starts":3,"utility_margin_multiplier":2.0,
    "control_equivalence_floor":1e-12,"rank_fd_steps":[1e-4,5e-5],
    "rank_relative_singular_threshold":1e-6}


def rows(value=7.2):
    return [{"animal_id":str(i//2),"raw_log10":value} for i in range(20)]


def _small_ode_task(loss):
    from research.zhou_reference.model import ModelParameters,simulate
    result=simulate(np.full(18,0.01),0,[1.0,2.0],ModelParameters(loss=loss),supplemented=True)
    return {"loss":loss,"target":result.target_log_change.tolist(),"resources":result.resources.tolist(),
            "accounting":result.diagnostics}


class FitTests(unittest.TestCase):
    def test_recovers_identifiable_loss_and_preserves_all_starts(self):
        report = runner.fit_control(rows(),np.full(20,8.0),
            lambda loss,tighter:np.full(20,8.0-loss),SETTINGS,1e-6)
        self.assertTrue(report["passed"])
        self.assertEqual(report["successful_starts"],3)
        self.assertAlmostEqual(report["best_loss"],.8,places=6)
        for start in report["starts"]:
            self.assertAlmostEqual(2*start["optimizer_cost"],start["objective"],places=12)
            self.assertGreaterEqual(start["residual_calls"],start["nfev"])
            self.assertEqual(start["rank_diagnostic"]["rank"],1)

    def test_equal_animal_then_equal_endpoint_objective(self):
        unequal = [{"animal_id":"a"},{"animal_id":"a"},{"animal_id":"b"}]
        np.testing.assert_allclose(runner.residual_weights(unequal)**2,[.25,.25,.5])

    def test_persistence_cannot_pass_strict_improvement(self):
        report = runner.fit_control(rows(8.0),np.full(20,8.0),
            lambda loss,tighter:np.full(20,8.0-loss),SETTINGS,1e-6)
        self.assertFalse(report["passed"])
        self.assertIn("no_strict_control_improvement_beyond_numerical_margin",report["failures"])

    def test_rank_zero_does_not_become_identifiable_from_good_fit(self):
        report = runner.fit_control(rows(),np.full(20,8.0),
            lambda loss,tighter:np.full(20,7.2),SETTINGS,1e-6)
        self.assertFalse(report["passed"])
        self.assertIn("unstable_or_deficient_control_jacobian_rank",report["failures"])
        self.assertEqual(len(report["equivalent_losses"]),3)

    def test_numerical_failures_are_retained_and_do_not_use_penalty_residuals(self):
        def predictor(loss,tighter):
            if loss > 1.5:
                raise RuntimeError("intentional numerical failure")
            return np.full(20,8.0-loss)
        report = runner.fit_control(rows(),np.full(20,8),predictor,SETTINGS,1e-6)
        self.assertFalse(report["passed"])
        self.assertEqual(len(report["starts"]),3)
        self.assertIn("intentional",report["starts"][2]["failure"])
        self.assertIn("fewer_than_required_successful_starts",report["failures"])

    def test_tighter_prediction_disagreement_blocks_freeze(self):
        report = runner.fit_control(rows(),np.full(20,8),
            lambda loss,tighter:np.full(20,8-loss+(1e-4 if tighter else 0)),SETTINGS,1e-6)
        self.assertFalse(report["passed"])
        self.assertIn("control_prediction_convergence_failed",report["failures"])

    def test_active_upper_bound_is_reported(self):
        report = runner.fit_control(rows(5),np.full(20,8),
            lambda loss,tighter:np.full(20,8-loss),SETTINGS,1e-6)
        self.assertTrue(report["passed"])
        self.assertAlmostEqual(report["best_loss"],2,places=6)
        self.assertTrue(any(s["active_bounds"]==[1] for s in report["starts"]))


class FreezeAndRoleTests(unittest.TestCase):
    def test_parallel_actual_ode_jobs_match_sequential_and_preserve_order(self):
        inputs=[1.8,0.2,1.0]
        sequential=list(runner._ordered_work(_small_ode_task,inputs,1))
        parallel=list(runner._ordered_work(_small_ode_task,inputs,2))
        self.assertEqual(sequential,parallel)
        self.assertEqual([r["loss"] for r in parallel],inputs)
        with self.assertRaisesRegex(ValueError,"1 through 8"):
            list(runner._ordered_work(_small_ode_task,inputs,9))

    def test_stale_engineering_evidence_cannot_lock(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/"protocol.json").write_text(json.dumps({"status":"specification_frozen"}))
            (root/"checks.json").write_text(json.dumps({"passed":True,"code_hashes":{}}))
            with patch.object(runner,"validate_protocol",return_value=True):
                with self.assertRaisesRegex(ValueError,"exact executable"):
                    runner.freeze_specification(root/"protocol.json",root/"checks.json",root/"lock.json")
            self.assertFalse((root/"lock.json").exists())

    def test_code_or_protocol_changes_reject_existing_lock(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            p=root/"protocol.json"
            p.write_text(json.dumps({"status":"specification_frozen"}))
            c=root/"checks.json"
            c.write_text(json.dumps({"passed":True,"code_hashes":runner._code_hashes()}))
            lock=root/"lock.json"
            with patch.object(runner,"validate_protocol",return_value=True):
                runner.freeze_specification(p,c,lock)
            runner._check_lock(p,lock)
            with patch.object(runner,"_code_hashes",return_value={"altered":"bytes"}):
                with self.assertRaisesRegex(ValueError,"source changed"):
                    runner._check_lock(p,lock)
            p.write_text(p.read_text()+"\n")
            with self.assertRaisesRegex(ValueError,"Protocol changed"):
                runner._check_lock(p,lock)

    def test_failed_calibration_cannot_open_evaluation_or_create_freeze(self):
        from research.zhou_reference import data
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            protocol=root/"p.json"
            protocol.write_text("{}")
            lock=root/"lock.json"
            lock.write_text("{}")
            baseline={"source_hashes":{},"source_fingerprint":"x"}
            report={"passed":False,"evaluation_opened":False}
            with patch.object(runner,"_check_lock",return_value={}), \
                    patch.object(data,"load_baselines",return_value=baseline), \
                    patch.object(data,"load_calibration",return_value=()), \
                    patch.object(data,"load_evaluation",side_effect=AssertionError("evaluation opened")) as evaluation, \
                    patch.object(data,"role_manifest",return_value={}), \
                    patch.object(runner,"calibrate_scenarios",return_value=report):
                status=runner.main(["calibrate","--sources",tmp,"--protocol",str(protocol),
                    "--lock",str(lock),"--output",str(root/"output")])
            self.assertEqual(status,2)
            evaluation.assert_not_called()
            self.assertTrue((root/"output/source_input_lock.json").exists())
            self.assertFalse((root/"output/executable_freeze.json").exists())

    def test_evaluation_requires_both_control_and_full_cohort_gates(self):
        with self.assertRaisesRegex(ValueError,"Calibration gate"):
            runner.evaluate_scenarios({},[],[],{}, {"passed":False},{"passed":True})
        with self.assertRaisesRegex(ValueError,"Full-cohort engineering"):
            runner.evaluate_scenarios({},[],[],{}, {"passed":True},{"passed":False})

    def test_control_cache_preserves_only_mathematical_equivalence(self):
        from research.zhou_reference.model import ModelParameters
        from dataclasses import asdict
        s={"params":asdict(ModelParameters()),"chronology":"age_derived","zero_policy":"absent",
           "zero_relative":{"absent":0,"low":1e-4,"high":1e-2},"biomass_scale":1,"abundance_scale":1e10}
        other=copy.deepcopy(s)
        other["params"]["qmax_G_target"]=4
        other["params"]["input_G"]=1
        self.assertEqual(runner._control_key(s),runner._control_key(other))
        other["params"]["qmax_R_target"]=2
        self.assertNotEqual(runner._control_key(s),runner._control_key(other))


class ContrastTests(unittest.TestCase):
    def test_pre_evaluation_convergence_mask_never_scores_missing_final(self):
        from research.data_readiness.map_benchmark import COMMUNITY
        animals=[]
        for iteration in ("Iteration 1","Iteration 2"):
            for supplemented in (False,True):
                animals.append({"animal_id":str(len(animals)),"experiment":iteration,"community":COMMUNITY,
                    "diet":"Mirpur-18 + Glucomannan" if supplemented else "Mirpur-18",
                    "sex":"Female","mouse_id":"1"})
        animals.append({"animal_id":"missing","experiment":"Iteration 2","community":COMMUNITY,
            "diet":"Mirpur-18 + Glucomannan","sex":"Male","mouse_id":"14"})
        first={a["animal_id"]:{"predictions":[8.,8.]} for a in animals}
        second=copy.deepcopy(first)
        second["missing"]["predictions"][1]=100.
        self.assertEqual(runner._outcome_free_contrast_change({"animals":animals},first,second),0.0)
        second["missing"]["predictions"][0]=9.
        self.assertAlmostEqual(runner._outcome_free_contrast_change({"animals":animals},first,second),.5)

    def test_protocol_numerical_config_cannot_silently_disagree_with_solver(self):
        protocol=json.loads((Path(runner.__file__).parent/"protocol_v05.json").read_text())
        self.assertTrue(runner.validate_protocol(protocol))
        protocol["numerical"]["rtol"]=1e-3
        with self.assertRaisesRegex(ValueError,"implemented solver"):
            runner.validate_protocol(protocol)

    def test_masks_and_calibration_assisted_labels(self):
        individual=[]
        for iteration in ("Iteration 1","Iteration 2"):
            for phase in ("P35","P53"):
                for treated in (False,True):
                    # Only iteration2 treated P53 loses one endpoint.
                    count=1 if iteration=="Iteration 2" and phase=="P53" and treated else 2
                    for _ in range(count):
                        individual.append({"experiment":iteration,"phase_label":phase,
                            "supplemented":treated,"baseline_log10":8+treated*.1,
                            "observed_log10":8+treated*.3,"predicted_log10":8+treated*.4,
                            "no_increment_log10":8+treated*.2})
        contrasts=runner.contrast_rows(individual,{})
        self.assertEqual(len(contrasts),16)
        for c in contrasts:
            self.assertEqual(c["n_control"],2)
            self.assertEqual(c["result_role"],"calibration_assisted" if c["experiment"]=="Iteration 1" else "retrospective_transfer")
            if c["comparator"]=="baseline_persistence":
                self.assertAlmostEqual(c["predicted_baseline_adjusted"],0)
            if c["comparator"]=="model":
                self.assertAlmostEqual(c["signed_log10_error"],.1)
                self.assertAlmostEqual(c["predicted_fold"],10**.4)
        final=[c for c in contrasts if c["experiment"]=="Iteration 2" and c["phase_label"]=="P53"]
        self.assertTrue(all(c["n_treated"]==1 for c in final))


if __name__ == "__main__":
    unittest.main()
