"""Control-only calibration and separately gated retrospective reporting.

Own-model numerical rules are prescribed by the versioned protocol, as described
in DOCUMENTATION.md. Evaluation loading occurs only in the evaluation entrypoint.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import multiprocessing
import platform
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import openpyxl
import scipy
from scipy.optimize import least_squares


def _stamp():
    return datetime.now(timezone.utc).isoformat()


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _environment():
    return {"python": sys.version, "executable": sys.executable,
            "numpy": np.__version__, "scipy": scipy.__version__, "openpyxl": openpyxl.__version__,
            "platform": platform.platform()}


def _code_hashes():
    package = Path(__file__).parent
    files = [package/name for name in ("model.py", "data.py", "runner.py")]
    readiness = package.parent/"data_readiness"
    files.extend(sorted(readiness.glob("*.py")))
    files.extend(sorted((readiness/"provenance").glob("*.json")))
    files.append(package.parent.parent/"species_data.csv")
    files.append(package.parent.parent/"audit_benchmark_table.py")
    return {str(path.relative_to(package.parent.parent)): _hash(path) for path in files}


def validate_protocol(protocol):
    """Reject configuration/code disagreement before any optimization."""
    from .model import ModelParameters, ACCOUNTING_TOLERANCE, LOG_CONVERGENCE_TOLERANCE
    expected = {"method":"DOP853","rtol":1e-8,"atol":1e-10,"max_step":1.0,
                "tight_rtol":1e-10,"tight_atol":1e-12,"tight_max_step":0.5,
                "max_log10_change":LOG_CONVERGENCE_TOLERANCE,
                "mass_balance_relative_tolerance":ACCOUNTING_TOLERANCE,
                "stage_domain_retry_max_steps":[1.0,0.5,0.25],
                "tight_stage_domain_retry_max_steps":[0.5,0.25,0.125],
                "stage_domain_retry_trigger":"finite_nonpositive_monod_denominator_in_internal_stage"}
    if protocol.get("protocol_version") != "0.5.1":
        raise ValueError("This executable requires the v0.5.1 numerical amendment; use historical code for v0.5")
    if protocol["numerical"] != expected:
        raise ValueError("Protocol numerical settings differ from implemented solver contract")
    c = protocol["calibration"]
    if (c["parameter"] != "loss" or c["lower"] != 0 or c["upper"] != 2
            or c["starts"] != [0.2,1.0,1.8] or c["required_successful_starts"] != 3):
        raise ValueError("Loss calibration differs from the one-parameter v0.5 contract")
    scenarios = protocol["scenarios"]
    if len(scenarios) != 96 or len({s["scenario_id"] for s in scenarios}) != 96:
        raise ValueError("The full locked 96-condition design is required")
    for scenario in scenarios:
        ModelParameters(**scenario["params"])
        if scenario["chronology"] not in ("age_derived","phase_labels"):
            raise ValueError("Unknown scenario chronology")
        if scenario["zero_policy"] not in ("absent","low","high"):
            raise ValueError("Unknown source-zero counterfactual")
        if scenario["zero_relative"] != {"absent":0.0,"low":1e-4,"high":1e-2}:
            raise ValueError("Source-zero scenario differs from the specified alternatives")
        if not all(math.isfinite(scenario[k]) and scenario[k]>0 for k in ("biomass_scale","abundance_scale")):
            raise ValueError("Invalid initialization scale")
    for observation in protocol["observation_maps"]:
        if (not math.isfinite(observation["treated_endpoint_factor"])
                or observation["treated_endpoint_factor"] <= 0 or observation["ramp_days"] != 1.0):
            raise ValueError("Invalid observation map")
    return True


def freeze_specification(protocol_path, engineering_path, lock_path):
    """Persist an exact specification/code lock without opening outcome inputs."""
    protocol = json.loads(Path(protocol_path).read_text())
    validate_protocol(protocol)
    if protocol.get("status") not in ("frozen", "specification_frozen"):
        raise ValueError("The scientific specification must be marked frozen before locking")
    report = json.loads(Path(engineering_path).read_text())
    if report.get("passed") is not True:
        raise ValueError("Engineering report does not establish a passing gate")
    if report.get("code_hashes") != _code_hashes():
        raise ValueError("Engineering report must bind passing tests to the exact executable sources")
    lock = {"locked_at": _stamp(), "protocol_sha256": _hash(protocol_path),
            "code_hashes": _code_hashes(), "environment": _environment(),
            "engineering_report_sha256": _hash(engineering_path),
            "engineering_report": report, "evaluation_opened": False}
    if Path(lock_path).exists():
        raise FileExistsError("Do not overwrite an existing scientific specification lock")
    _write(lock_path, lock)
    return lock


def _check_lock(protocol_path, lock_path):
    lock = json.loads(Path(lock_path).read_text())
    if lock["protocol_sha256"] != _hash(protocol_path):
        raise ValueError("Protocol changed after its specification lock")
    if lock["code_hashes"] != _code_hashes():
        raise ValueError("Executable source changed after its specification lock")
    if lock["environment"] != _environment():
        raise ValueError("Execution environment changed after its specification lock")
    if lock["engineering_report"].get("passed") is not True:
        raise ValueError("Engineering gate is not passing")
    return lock


def residual_weights(rows):
    """Equal animals, then equal available follow-ups within each animal."""
    counts = Counter(r["animal_id"] for r in rows)
    if not counts:
        raise ValueError("Empty calibration data")
    return np.array([1 / math.sqrt(len(counts) * counts[r["animal_id"]]) for r in rows])


def finite_difference_rank(residual_function, normalized, steps, threshold):
    results = []
    for step in steps:
        lo, hi = max(0.0, normalized-step), min(1.0, normalized+step)
        derivative = (residual_function(hi)-residual_function(lo)) / (hi-lo)
        if not np.isfinite(derivative).all():
            raise ValueError("Nonfinite finite-difference derivative")
        singular = np.linalg.svd(derivative[:, None], compute_uv=False)
        largest = float(singular[0])
        rank = int(np.count_nonzero(singular > threshold*largest)) if largest > 0 else 0
        results.append({"step": step, "rank": rank, "singular_values": singular.tolist()})
    return {"steps": results, "stable": len({r["rank"] for r in results}) == 1,
            "rank": results[0]["rank"]}


def fit_control(rows, baseline_predictions, predictor, settings, tolerance):
    """Fit a single loss coefficient; predictor(loss, tighter) cannot see evaluation rows.

    A failed integration fails that optimizer start; no artificial residual is
    substituted. Every start is retained, including failures and active bounds.
    """
    weights = residual_weights(rows)
    observed = np.array([r["raw_log10"] for r in rows], dtype=float)
    baseline = np.asarray(baseline_predictions, dtype=float)
    if baseline.shape != observed.shape or not np.isfinite(observed).all():
        raise ValueError("Malformed control observations/baseline predictions")
    lower, upper = settings["lower"], settings["upper"]
    if not (math.isfinite(lower) and math.isfinite(upper) and lower < upper):
        raise ValueError("Loss bounds must be finite and ordered")
    cache = {}

    def predictions(z, tighter=False):
        key = (float(z), bool(tighter))
        if key not in cache:
            value = np.asarray(predictor(lower+(upper-lower)*z, tighter), dtype=float)
            if value.shape != observed.shape or not np.isfinite(value).all():
                raise ValueError("Undefined/nonfinite control prediction")
            cache[key] = value
        return cache[key]

    def residual(z, tighter=False):
        return (predictions(float(z), tighter)-observed)*weights

    def objective(z, tighter=False):
        r = residual(z, tighter)
        return float(r @ r)

    starts = []
    for start in settings["starts"]:
        entry = {"initial_loss": start, "residual_calls": 0, "success": False}

        def fun(z):
            entry["residual_calls"] += 1
            return residual(float(z[0]))

        try:
            fit = least_squares(fun, [(start-lower)/(upper-lower)], bounds=([0.0], [1.0]),
                method=settings["method"], loss=settings["loss"], jac=settings["jac"],
                tr_solver=settings["tr_solver"], x_scale=settings["x_scale"],
                ftol=settings["ftol"], xtol=settings["xtol"], gtol=settings["gtol"],
                max_nfev=settings["max_nfev"])
            z = float(fit.x[0])
            entry.update(success=bool(fit.success), normalized_loss=z,
                loss=lower+(upper-lower)*z, objective=objective(z),
                optimizer_cost=float(fit.cost), nfev=int(fit.nfev),
                njev=int(fit.njev) if fit.njev is not None else None,
                status=int(fit.status), message=fit.message,
                active_bounds=fit.active_mask.tolist(), optimality=float(fit.optimality))
            if fit.success:
                entry["tight_objective"] = objective(z, True)
                entry["max_prediction_change"] = float(np.max(np.abs(predictions(z)-predictions(z, True))))
                entry["rank_diagnostic"] = finite_difference_rank(
                    lambda v: residual(v, True), z, settings["rank_fd_steps"],
                    settings["rank_relative_singular_threshold"])
                entry["control_predictions"] = predictions(z).tolist()
                entry["control_residuals"] = (predictions(z)-observed).tolist()
        except Exception as exc:
            entry.update(success=False, failure_type=type(exc).__name__, failure=str(exc))
            if getattr(exc,"integration_attempts",None):
                entry["integration_attempts"]=exc.integration_attempts
        starts.append(entry)
    successful = [s for s in starts if s["success"]]
    persistence = float(np.sum(((baseline-observed)*weights)**2))
    report = {"starts": starts, "persistence_objective": persistence,
              "successful_starts": len(successful), "passed": False,
              "evaluation_opened": False, "failures": []}
    if len(successful) < settings["required_successful_starts"]:
        report["failures"].append("fewer_than_required_successful_starts")
    if not successful:
        return report
    best = min(successful, key=lambda s: s["objective"])
    margin = max(abs(s["objective"]-s["tight_objective"]) for s in successful)
    equivalence = max(settings["control_equivalence_floor"], settings["utility_margin_multiplier"]*margin)
    equivalent = [s for s in successful if s["objective"]-best["objective"] <= equivalence]
    improvement = persistence-best["objective"]
    if not improvement > settings["utility_margin_multiplier"]*margin:
        report["failures"].append("no_strict_control_improvement_beyond_numerical_margin")
    if any(s["max_prediction_change"] > tolerance for s in successful):
        report["failures"].append("control_prediction_convergence_failed")
    if any(not s["rank_diagnostic"]["stable"] or s["rank_diagnostic"]["rank"] != 1 for s in equivalent):
        report["failures"].append("unstable_or_deficient_control_jacobian_rank")
    # No equivalent solution is selected using evaluation agreement.
    report.update(best_loss=best["loss"], best_objective=best["objective"],
        objective_numerical_margin=margin, improvement_over_persistence=improvement,
        equivalence_tolerance=equivalence,
        equivalent_losses=sorted({s["loss"] for s in equivalent}),
        equivalent_loss_range=[min(s["loss"] for s in equivalent), max(s["loss"] for s in equivalent)],
        passed=not report["failures"])
    return report


def _control_key(scenario):
    # G=0 and zero G input in every control: only its R kinetics affect controls.
    # This equivalence is exact, not a similarity selected from fitted results.
    if scenario["params"]["initial_G"] != 0:
        return _digest(scenario)
    p = scenario["params"]
    keys = ("qmax_R_target", "qmax_R_background", "K_R", "Y_R_target",
            "Y_R_background", "input_R", "washout_R", "initial_R")
    return _digest({"parameters": {k:p[k] for k in keys},
        "chronology": scenario["chronology"], "zero_policy": scenario["zero_policy"],
        "zero_relative": scenario["zero_relative"], "biomass_scale": scenario["biomass_scale"],
        "abundance_scale": scenario["abundance_scale"]})


def _model_predictions(baselines, scenario, loss, animals, tighter=False, no_increment=False,integration_audit=None):
    from .data import get_times, initialize
    from .model import simulate
    result = {}
    parameters = dict(scenario["params"], loss=loss)
    for a in animals:
        times = get_times(a["experiment"], scenario["chronology"])
        supplemented = a["diet"] != "Mirpur-18" and not no_increment
        context={"animal_id":a["animal_id"],"scenario_id":scenario["scenario_id"],
                 "loss":loss,"tighter":tighter,"no_increment":no_increment}
        try:
            simulation = simulate(initialize(a, scenario), baselines["target_index"], times,
                                  parameters, supplemented=supplemented, tighter=tighter)
        except Exception as exc:
            # A failed full-cohort prediction must remain attributable even
            # outside calibration, where no shared integration audit is passed.
            exc.integration_context=context
            if integration_audit is not None:
                attempts=getattr(exc,"integration_attempts",[])
                _audit_integration(integration_audit,context,attempts,False,
                                   {"failure_type":type(exc).__name__,"failure":str(exc)})
            raise
        if integration_audit is not None:
            _audit_integration(integration_audit,context,simulation.diagnostics.get("integration_attempts",[]),True)
        result[a["animal_id"]] = {
            "predictions": (a["target_log10"] + simulation.target_log_change).tolist(),
            "times": list(times), "diagnostics": simulation.diagnostics}
    return result


def _audit_integration(audit,context,attempts,success,failure=None):
    """Count every solve and retain all failed/retried stages, including recoveries."""
    audit["simulation_calls"]=audit.get("simulation_calls",0)+1
    audit["solver_attempts"]=audit.get("solver_attempts",0)+len(attempts)
    failed_attempts=sum(not attempt["success"] for attempt in attempts)
    audit["failed_solver_attempts"]=audit.get("failed_solver_attempts",0)+failed_attempts
    if success:
        category="recovered_retry_simulations" if failed_attempts else "successful_without_retry_simulations"
    else:
        category="failed_simulations"
    audit[category]=audit.get(category,0)+1
    if failed_attempts or not success:
        audit.setdefault("retried_or_failed_simulations",[]).append(
            dict(context,success=success,integration_attempts=attempts,**(failure or {})))


def _row_predictions(rows, predictions):
    from .data import phase_index
    return np.array([predictions[r["animal_id"]]["predictions"][phase_index(r["phase_label"])] for r in rows])


def _fit_control_work(payload):
    key,baselines,rows,scenario,settings,tolerance = payload
    lookup = {a["animal_id"]:a for a in baselines["animals"]}
    controls = [a for a in baselines["animals"] if a["animal_id"] in {r["animal_id"] for r in rows}]
    y0 = [lookup[r["animal_id"]]["target_log10"] for r in rows]
    integration_audit={"retried_or_failed_simulations":[]}
    def predict(loss,tighter):
        return _row_predictions(rows,_model_predictions(baselines,scenario,loss,controls,tighter,
                                                       integration_audit=integration_audit))
    try:
        report=fit_control(rows,y0,predict,settings,tolerance)
    except Exception as exc:
        report={"passed":False,"starts":[],"successful_starts":0,"evaluation_opened":False,
            "failures":["control_worker_failure"],"failure_type":type(exc).__name__,"failure":str(exc)}
        if getattr(exc,"integration_attempts",None):
            report["integration_attempts"]=exc.integration_attempts
    report["integration_audit"]=integration_audit
    return key,report


def _ordered_work(function,payloads,workers):
    if isinstance(workers,bool) or not isinstance(workers,int) or not 1 <= workers <= 8:
        raise ValueError("workers must be an integer from 1 through 8")
    if workers == 1:
        yield from map(function,payloads)
    else:
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context("spawn")) as pool:
            yield from pool.map(function,payloads)


def calibrate_scenarios(baselines, rows, protocol, progress=None, workers=1):
    from .data import validate_calibration
    validate_calibration(rows, baselines)
    fits, scenarios = {}, []
    unique={}
    for s in protocol["scenarios"]:
        key = _control_key(s)
        unique.setdefault(key,s)
    payloads=[(key,baselines,rows,s,protocol["calibration"],protocol["numerical"]["max_log10_change"])
              for key,s in unique.items()]
    for key,fit in _ordered_work(_fit_control_work,payloads,workers):
        fits[key]=fit
        if progress:
            progress({"control_fit_count":len(fits),"scenario_id":unique[key]["scenario_id"],
                      "passed":fit["passed"],"failures":fit["failures"]})
    for s in protocol["scenarios"]:
        key=_control_key(s)
        scenarios.append({"scenario_id": s["scenario_id"], "control_fit_key": key, "passed": fits[key]["passed"]})
    return {"calibration_started_without_evaluation": True, "evaluation_opened": False,
            "source_fingerprint": baselines["source_fingerprint"],
            "calibration_count": len(rows), "calibration_rows": list(rows),
            "fits": fits, "scenarios": scenarios,"workers":workers,
            "passed": all(s["passed"] for s in scenarios)}


def _verify_scenario_work(payload):
    baselines,protocol,calibration_report,scenario = payload
    try:
        return verify_prediction_set(baselines,dict(protocol,scenarios=[scenario]),calibration_report)
    except Exception as exc:
        failure={"scenario_id":scenario["scenario_id"],"passed":False,
                 "failure_type":type(exc).__name__,"failure":str(exc)}
        if getattr(exc,"integration_attempts",None):
            failure["integration_attempts"]=exc.integration_attempts
        if getattr(exc,"integration_context",None):
            failure["integration_context"]=exc.integration_context
        return {"passed":False,"runs":[],"evaluation_opened":False,
            "failures":[failure]}


def verify_prediction_set(baselines, protocol, calibration_report, progress=None, workers=1):
    """Outcome-free, full-cohort numerical gate before opening evaluation values."""
    if not calibration_report.get("passed"):
        raise ValueError("Cannot advance failed control calibration")
    if workers != 1:
        runs,failures=[],[]
        payloads=[(baselines,protocol,calibration_report,s) for s in protocol["scenarios"]]
        for index,result in enumerate(_ordered_work(_verify_scenario_work,payloads,workers)):
            runs.extend(result["runs"])
            failures.extend(result["failures"])
            if progress:
                progress({"prediction_scenario":protocol["scenarios"][index]["scenario_id"],
                          "numerical_failures":len(failures)})
        return {"passed":not failures,"runs":runs,"failures":failures,
                "evaluation_opened":False,"workers":workers}
    fit_keys = {s["scenario_id"]:s["control_fit_key"] for s in calibration_report["scenarios"]}
    tolerance = protocol["numerical"]["max_log10_change"]
    runs, failures = [], []
    no_g_cache = {}
    for s in protocol["scenarios"]:
        fit = calibration_report["fits"][fit_keys[s["scenario_id"]]]
        for member, loss in enumerate(fit["equivalent_losses"]):
            run = {"scenario_id":s["scenario_id"],"equivalent_member":member,"loss":loss}
            try:
                pred = _model_predictions(baselines,s,loss,baselines["animals"])
                tight = _model_predictions(baselines,s,loss,baselines["animals"],True)
                key = (_control_key(s),loss)
                if key not in no_g_cache:
                    no_g_cache[key] = (
                        _model_predictions(baselines,s,loss,baselines["animals"],no_increment=True),
                        _model_predictions(baselines,s,loss,baselines["animals"],True,True))
                no_g,no_g_tight = no_g_cache[key]
                neutral = dict(s,params=dict(s["params"],input_G=0.0))
                label_check = _model_predictions(baselines,neutral,loss,baselines["animals"])
                if any(no_g[k]["predictions"] != label_check[k]["predictions"] for k in no_g):
                    raise ValueError("Identical inputs differ under diet labels")
                def max_change(a,b):
                    return max(float(np.max(np.abs(np.array(a[k]["predictions"])-b[k]["predictions"]))) for k in a)
                change,no_change = max_change(pred,tight),max_change(no_g,no_g_tight)
                contrast_change = max(_outcome_free_contrast_change(baselines,pred,tight),
                                      _outcome_free_contrast_change(baselines,no_g,no_g_tight))
                if max(change,no_change,contrast_change) > tolerance:
                    raise ValueError("Prediction/contrast numerical convergence failed")
                run.update(predictions=pred,tight_predictions=tight,no_increment=no_g,
                    tight_no_increment=no_g_tight,max_prediction_change=change,
                    max_no_increment_change=no_change,max_contrast_change=contrast_change,
                    exact_no_effect_check=True,passed=True)
                runs.append(run)
            except Exception as exc:
                failure=dict(run,passed=False,failure_type=type(exc).__name__,failure=str(exc))
                if getattr(exc,"integration_attempts",None):
                    failure["integration_attempts"]=exc.integration_attempts
                if getattr(exc,"integration_context",None):
                    failure["integration_context"]=exc.integration_context
                failures.append(failure)
        if progress:
            progress({"prediction_scenario":s["scenario_id"],"numerical_failures":len(failures)})
    return {"passed":not failures,"runs":runs,"failures":failures,"evaluation_opened":False,"workers":workers}


def _outcome_free_contrast_change(baselines,first,second):
    """Known endpoint eligibility is distinct from its unopened measured outcome."""
    from .data import META
    from research.data_readiness.map_benchmark import MISSING_KEY
    differences=[]
    for experiment in ("Iteration 1","Iteration 2"):
        for endpoint in (0,1):
            by_arm={True:[],False:[]}
            for animal in baselines["animals"]:
                if animal["experiment"] != experiment:
                    continue
                if endpoint == 1 and tuple(animal[k] for k in META) == MISSING_KEY:
                    continue
                key=animal["animal_id"]
                by_arm[animal["diet"] != "Mirpur-18"].append(
                    first[key]["predictions"][endpoint]-second[key]["predictions"][endpoint])
            if not all(by_arm.values()):
                raise ValueError("Incomplete predeclared contrast cohort")
            differences.append(abs(float(np.mean(by_arm[True])-np.mean(by_arm[False]))))
    return max(differences)


def summarize_sensitivity(contrasts,individual):
    """Conditional ranges and matched contrasts, never probabilities or selection."""
    model=[r for r in contrasts if r["comparator"]=="model"]
    summaries=[]
    for experiment,phase in sorted({(r["experiment"],r["phase_label"]) for r in model}):
        group=[r for r in model if (r["experiment"],r["phase_label"])==(experiment,phase)]
        predictions=[r["predicted_log10_contrast"] for r in group]
        errors=[r["absolute_log10_error"] for r in group]
        summaries.append({"experiment":experiment,"phase_label":phase,
            "prediction_range":[min(predictions),max(predictions)],"absolute_error_range":[min(errors),max(errors)],
            "positive_prediction_count":sum(v>0 for v in predictions),
            "negative_prediction_count":sum(v<0 for v in predictions),
            "zero_prediction_count":sum(v==0 for v in predictions),"reported_members":len(group),
            "interpretation":"Counts describe this finite analyst design, not probabilities or robustness outside its range"})
    evaluation=[r for r in individual if r["role"]=="evaluation"]
    individual_summaries=[]
    groups={tuple(r[k] for k in ("scenario_id","observation_map","equivalent_member")) for r in evaluation}
    for key in sorted(groups):
        subset=[r for r in evaluation if tuple(r[k] for k in ("scenario_id","observation_map","equivalent_member"))==key]
        for experiment in ("Iteration 1","Iteration 2"):
            values=[r["residual_log10"] for r in subset if r["experiment"]==experiment]
            individual_summaries.append({"scenario_id":key[0],"observation_map":key[1],"equivalent_member":key[2],
                "experiment":experiment,"evaluation_count":len(values),
                "evaluation_MSE":float(np.mean(np.square(values))),"evaluation_MAE":float(np.mean(np.abs(values)))})
    # Compare finite panels within the same chronology, zero and observation
    # condition; intervals retain every numerically equivalent fitted loss.
    intervals={}
    for row in model:
        key=tuple(row[k] for k in ("chronology","zero_policy","observation_map","experiment","phase_label","panel"))
        intervals.setdefault(key,[]).append(row["predicted_log10_contrast"])
    ranges={k:(min(v),max(v)) for k,v in intervals.items()}
    effects,interactions=[],[]
    for key,(low,high) in sorted(ranges.items()):
        base=key[:-1]+("center",)
        if key[-1] != "center" and base in ranges:
            c_low,c_high=ranges[base]
            effects.append({"chronology":key[0],"zero_policy":key[1],"observation_map":key[2],
                "experiment":key[3],"phase_label":key[4],"panel":key[5],
                "change_from_center_range":[low-c_high,high-c_low]})
    for prefix in sorted({k[:-1] for k in ranges}):
        for background in ("none","stronger"):
            names=("center","exposure_high","background_"+background,"exposure_high_background_"+background)
            if all(prefix+(n,) in ranges for n in names):
                a,b,c,d=[ranges[prefix+(n,)] for n in names]
                interactions.append({"chronology":prefix[0],"zero_policy":prefix[1],"observation_map":prefix[2],
                    "experiment":prefix[3],"phase_label":prefix[4],"background":background,
                    "finite_interaction_range":[d[0]-b[1]-c[1]+a[0],d[1]-b[0]-c[0]+a[1]],
                    "definition":"joint - exposure_high - background_alternative + center"})
    return {"endpoint_ranges":summaries,"individual_evaluation_errors":individual_summaries,
            "matched_panel_changes":effects,"prespecified_interactions":interactions}


def _csv(path, rows):
    if not rows:
        raise ValueError("No reporting rows")
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def contrast_rows(individual, meta):
    """Four endpoint-by-iteration contrasts; all comparators use exact same masks."""
    results = []
    groups = sorted({(r["experiment"],r["phase_label"]) for r in individual})
    for experiment, phase in groups:
        group = [r for r in individual if (r["experiment"],r["phase_label"]) == (experiment, phase)]
        treated = [r for r in group if r["supplemented"]]
        controls = [r for r in group if not r["supplemented"]]
        if not treated or not controls:
            raise ValueError("Missing diet arm in endpoint mask")
        def delta(field):
            return float(np.mean([r[field] for r in treated])-np.mean([r[field] for r in controls]))
        observed, baseline = delta("observed_log10"), delta("baseline_log10")
        for comparator, prediction in (("model",delta("predicted_log10")),
            ("baseline_persistence",baseline), ("mechanistic_no_increment",delta("no_increment_log10")),
            ("zero_between_diet_contrast",0.0)):
            results.append(dict(meta, experiment=experiment, phase_label=phase,
                result_role="calibration_assisted" if experiment == "Iteration 1" else "retrospective_transfer",
                comparator=comparator, n_treated=len(treated), n_control=len(controls),
                observed_log10_contrast=observed, predicted_log10_contrast=prediction,
                observed_fold=10**observed, predicted_fold=10**prediction,
                signed_log10_error=prediction-observed, absolute_log10_error=abs(prediction-observed),
                observed_baseline_adjusted=observed-baseline,
                predicted_baseline_adjusted=prediction-baseline))
    return results


def evaluate_scenarios(baselines, calibration_rows, evaluation_rows, protocol, calibration_report, prediction_set):
    """Called only after immutable executable gate. Outcomes are never used to select."""
    from .data import phase_index
    if calibration_report.get("passed") is not True:
        raise ValueError("Calibration gate failed; evaluation is prohibited")
    if prediction_set.get("passed") is not True:
        raise ValueError("Full-cohort engineering gate failed; evaluation is prohibited")
    all_rows = list(calibration_rows) + list(evaluation_rows)
    lookup = {a["animal_id"]: a for a in baselines["animals"]}
    fit_keys = {s["scenario_id"]:s["control_fit_key"] for s in calibration_report["scenarios"]}
    individual, contrasts, failures, convergence = [], [], [], []
    tolerance = protocol["numerical"]["max_log10_change"]
    frozen_predictions = {(r["scenario_id"],r["equivalent_member"]):r for r in prediction_set["runs"]}
    for s in protocol["scenarios"]:
        fit = calibration_report["fits"][fit_keys[s["scenario_id"]]]
        for member, loss in enumerate(fit["equivalent_losses"]):
            try:
                frozen = frozen_predictions[(s["scenario_id"],member)]
                if frozen["loss"] != loss or frozen.get("passed") is not True:
                    raise ValueError("Frozen prediction member does not match calibration")
                pred,tight = frozen["predictions"],frozen["tight_predictions"]
                no_g,no_g_tight = frozen["no_increment"],frozen["tight_no_increment"]
                max_change = max(float(np.max(np.abs(np.array(pred[a]["predictions"])-tight[a]["predictions"]))) for a in pred)
                max_no_g_change = max(float(np.max(np.abs(np.array(no_g[a]["predictions"])-no_g_tight[a]["predictions"]))) for a in pred)
                if max(max_change,max_no_g_change) > tolerance:
                    raise ValueError("Evaluation/no-increment log prediction convergence failed")
                convergence.append({"scenario_id":s["scenario_id"], "equivalent_member":member,
                    "max_prediction_change":max_change, "max_no_increment_change":max_no_g_change})
                for observation in protocol["observation_maps"]:
                    meta = {"scenario_id":s["scenario_id"], "panel":s["panel"], "chronology":s["chronology"],
                        "zero_policy":s["zero_policy"], "observation_map":observation["id"],
                        "equivalent_member":member, "loss":loss}
                    current, tighter_current = [], []
                    for row in all_rows:
                        a = lookup[row["animal_id"]]
                        idx = phase_index(row["phase_label"])
                        supplemented = a["diet"] != "Mirpur-18"
                        t = pred[a["animal_id"]]["times"][idx]
                        shift = math.log10(observation["treated_endpoint_factor"])*min(t/observation["ramp_days"],1.0) if supplemented else 0.0
                        value = pred[a["animal_id"]]["predictions"][idx] + shift
                        # Observation sidecar describes the assigned arm, including
                        # its no-increment counterfactual, and never changes controls.
                        record = dict(meta,animal_id=a["animal_id"],experiment=a["experiment"],diet=a["diet"],
                            phase_label=row["phase_label"], role=row["role"], supplemented=supplemented,
                            baseline_log10=a["target_log10"],observed_log10=row["raw_log10"],
                            predicted_log10=value, residual_log10=value-row["raw_log10"],
                            no_increment_log10=no_g[a["animal_id"]]["predictions"][idx]+shift)
                        current.append(record)
                        tighter_current.append(dict(record,
                            predicted_log10=tight[a["animal_id"]]["predictions"][idx]+shift,
                            no_increment_log10=no_g_tight[a["animal_id"]]["predictions"][idx]+shift))
                    c = contrast_rows(current,meta)
                    ct = contrast_rows(tighter_current,meta)
                    contrast_change = max(abs(a["predicted_log10_contrast"]-b["predicted_log10_contrast"]) for a,b in zip(c,ct))
                    if contrast_change > tolerance:
                        raise ValueError("Endpoint contrast convergence failed")
                    individual.extend(current)
                    contrasts.extend(c)
            except Exception as exc:
                failures.append({"scenario_id":s["scenario_id"],"equivalent_member":member,"loss":loss,
                    "failure_type":type(exc).__name__,"failure":str(exc)})
    return {"individual":individual,"contrasts":contrasts,"failures":failures,
            "convergence":convergence,"passed":not failures,"evaluation_opened":True}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    lock = sub.add_parser("lock")
    lock.add_argument("--protocol",required=True)
    lock.add_argument("--engineering-report",required=True)
    lock.add_argument("--output",required=True)
    for command in ("calibrate","evaluate"):
        child = sub.add_parser(command)
        child.add_argument("--sources",required=True)
        child.add_argument("--protocol",required=True)
        child.add_argument("--lock",required=True)
        child.add_argument("--output",required=True)
        if command == "calibrate":
            child.add_argument("--workers",type=int,default=1,choices=range(1,9))
        if command == "evaluate":
            child.add_argument("--calibration",required=True)
    args = parser.parse_args(argv)
    if args.command == "lock":
        freeze_specification(args.protocol,args.engineering_report,args.output)
        return 0
    lock_record = _check_lock(args.protocol,args.lock)
    protocol = json.loads(Path(args.protocol).read_text())
    from .data import load_baselines, load_calibration
    baseline = load_baselines(args.sources)
    calibration = load_calibration(args.sources)
    output = Path(args.output)
    output.mkdir(parents=True,exist_ok=False)
    if args.command == "calibrate":
        from .data import role_manifest
        input_lock = dict(role_manifest(baseline,calibration),locked_at=_stamp(),lock_sha256=_hash(args.lock))
        _write(output/"source_input_lock.json",input_lock)
        report = calibrate_scenarios(baseline,calibration,protocol,lambda v: print(json.dumps(v),flush=True),workers=args.workers)
        report.update(completed_at=_stamp(),lock_sha256=_hash(args.lock),source_hashes=baseline["source_hashes"],
                      source_input_lock_sha256=_hash(output/"source_input_lock.json"))
        _write(output/"calibration.json",report)
        if report["passed"]:
            predictions = verify_prediction_set(baseline,protocol,report,lambda v: print(json.dumps(v),flush=True),workers=args.workers)
            _write(output/"frozen_predictions.json",predictions)
            if not predictions["passed"]:
                print(json.dumps({"full_cohort_engineering_passed":False,"evaluation_opened":False}),flush=True)
                return 2
            freeze = {"frozen_at":_stamp(),"lock_sha256":_hash(args.lock),
                "calibration_sha256":_hash(output/"calibration.json"),
                "source_input_lock_sha256":_hash(output/"source_input_lock.json"),
                "frozen_predictions_sha256":_hash(output/"frozen_predictions.json"),
                "source_fingerprint":baseline["source_fingerprint"],
                "evaluation_enabled":True,"environment":_environment(),"code_hashes":_code_hashes(),
                "calibration_gate_passed":True,"engineering_gate_passed":True,
                "protocol_sha256":lock_record["protocol_sha256"]}
            _write(output/"executable_freeze.json",freeze)
        print(json.dumps({"calibration_passed":report["passed"],"evaluation_opened":False}),flush=True)
        return 0 if report["passed"] else 2
    frozen_dir = Path(args.calibration)
    freeze = json.loads((frozen_dir/"executable_freeze.json").read_text())
    if not freeze.get("evaluation_enabled") or freeze["source_fingerprint"] != baseline["source_fingerprint"]:
        raise ValueError("Evaluation freeze missing/disabled or sources changed")
    if freeze["lock_sha256"] != _hash(args.lock) or freeze["calibration_sha256"] != _hash(frozen_dir/"calibration.json"):
        raise ValueError("Calibration or lock changed after executable freeze")
    report = json.loads((frozen_dir/"calibration.json").read_text())
    from .data import role_manifest
    input_lock = json.loads((frozen_dir/"source_input_lock.json").read_text())
    if freeze["source_input_lock_sha256"] != _hash(frozen_dir/"source_input_lock.json"):
        raise ValueError("Source input lock changed after freeze")
    current_roles = role_manifest(baseline,calibration)
    if any(input_lock[k] != v for k,v in current_roles.items()):
        raise ValueError("Derived baseline/calibration inputs changed after freeze")
    if freeze["frozen_predictions_sha256"] != _hash(frozen_dir/"frozen_predictions.json"):
        raise ValueError("Outcome-free frozen predictions changed after freeze")
    predictions = json.loads((frozen_dir/"frozen_predictions.json").read_text())
    # The evaluation loader is deliberately imported only beyond every gate.
    from .data import load_evaluation
    evaluation = load_evaluation(args.sources)
    role_manifest(baseline,calibration,evaluation)
    result = evaluate_scenarios(baseline,calibration,evaluation,protocol,report,predictions)
    if result["individual"]:
        _write(output/"sensitivity_summary.json",summarize_sensitivity(result["contrasts"],result["individual"]))
        _csv(output/"individual_predictions.csv",result.pop("individual"))
        _csv(output/"endpoint_contrasts.csv",result.pop("contrasts"))
    result.update(completed_at=_stamp(),executable_freeze_sha256=_hash(frozen_dir/"executable_freeze.json"),
                  source_fingerprint=baseline["source_fingerprint"],environment=_environment())
    _write(output/"evaluation.json",result)
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
