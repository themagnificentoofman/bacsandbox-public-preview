"""Retained v0.3 qualification and final freeze, with unresolved anomalies blocked.

This extends the prospective study admission layer; it does not weaken the
50/200/1000-pair conservative numerical equivalence rule. All choices are Own.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

from .protocol import digest
from .qualification import REQUIRED_CHECKS, compare_timesteps
from .qualification_run import _bind, _source, _test_identity, junit_checks
from .study_v03 import (PLAN_SHA256, _current_identity, _verify_sources, audit_cohort,
                        load_plan, read_sealed, seal, selected_protocol, verify_admission,
                        verify_grid, verify_selection, write_new)

V03_TESTS = (
    "test_v03_plan_rejects_rehashed_weakening",
    "test_v03_grid_registration_precedes_execution",
    "test_v03_eligibility_counts_only_u_both_populations",
    "test_v03_no_eligible_candidate_stays_blocked",
    "test_v03_fresh_pilot_does_not_require_h_births",
    "test_v03_freeze_rejects_unresolved_anomaly",
)


def fresh_pilot_decision(audit):
    u = [r for r in audit["by_arm"] if r["arm"] == "U"]
    h = [r for r in audit["by_arm"] if r["arm"] == "H"]
    seeds = list(range(5000, 5020))
    complete = (audit["seeds"] == seeds and audit["completed_pairs"] == 20
                and sorted(r["seed"] for r in u) == seeds and sorted(r["seed"] for r in h) == seeds)
    birth_runs = sum(r["births"]["focal"] > 0 and r["births"]["resident"] > 0 for r in u)
    culls = sum(r["capacity_culling_count"] for r in audit["by_arm"])
    return {"status": "passed" if complete and audit["status"] == "passed" and birth_runs >= 10 and culls == 0 else "failed",
            "seeds": seeds, "completed_pairs": audit["completed_pairs"], "U_denominator": len(u),
            "U_runs_with_both_population_births": birth_runs, "capacity_culling_count": culls,
            "require_H_births": False, "survival_count_gate": None}


def _anomalies_from_audit(label, audit):
    result = copy.deepcopy(audit["retained_attempts"])
    if audit["status"] != "passed":
        result.append({"id": label + ":accounting_or_replay", "kind": "failed_integrity_audit"})
    if any(r["capacity_culling_count"] for r in audit["by_arm"]):
        result.append({"id": label + ":capacity_culling", "kind": "capacity_culling"})
    for arm in ("U", "H"):
        rows = [r for r in audit["by_arm"] if r["arm"] == arm]
        survivors = sum(r["focal_count"] > 0 for r in rows)
        if survivors in (0, len(rows)):
            result.append({"id": f"{label}:{arm}:degenerate_survival", "kind": "degenerate_survival",
                           "surviving_runs": survivors, "denominator": len(rows)})
        for population in ("focal", "resident"):
            births = sum(r["births"][population] for r in rows)
            if births == 0:
                result.append({"id": f"{label}:{arm}:{population}:zero_births", "kind": "zero_births",
                               "births": 0, "population": population, "denominator": len(rows)})
    return result


def anomaly_inventory(selection, pilot=None, *, coarse=None, fine=None):
    selection = verify_selection(selection)
    sources = copy.deepcopy(selection["sources"])
    items = []
    for candidate in selection["candidates"]:
        items.extend(_anomalies_from_audit(candidate["id"], candidate["audit"]))
    paths = {"pilot": pilot, "coarse": coarse, "fine": fine}
    numerical_records = {}
    for label, path in paths.items():
        if path is not None:
            _, records, audit = audit_cohort(path, sources, prefix=label)
            items.extend(_anomalies_from_audit(label, audit))
            if label in ("coarse", "fine"):
                numerical_records[label] = records
    if numerical_records:
        if set(numerical_records) != {"coarse", "fine"}:
            raise ValueError("numerical anomaly review requires both resolutions")
        summary = {}
        for arm in ("U", "H"):
            indexed = {label: {r["seed"]: r for r in rows if r["arm"] == arm}
                       for label, rows in numerical_records.items()}
            if set(indexed["coarse"]) != set(indexed["fine"]):
                raise ValueError("numerical anomaly review has unpaired seeds")
            seeds = sorted(indexed["coarse"])
            summary[arm] = {"pairs": len(seeds), "pathwise_survival_flips": sum(
                (indexed["coarse"][s]["focal_count"] > 0) != (indexed["fine"][s]["focal_count"] > 0) for s in seeds),
                "mean_changes": {key: math.fsum(indexed["fine"][s]["diagnostics"][key] -
                    indexed["coarse"][s]["diagnostics"][key] for s in seeds) / len(seeds)
                    for key in ("births", "deaths", "uptake_pg", "exposure_contrast", "early_exposure_contrast")}}
        items.append({"id": "numerical:descriptive_diagnostics_review", "kind": "mandatory_numerical_diagnostic_review",
                      "summary": summary, "acceptance_threshold": None})
    if len({item["id"] for item in items}) != len(items):
        raise ValueError("duplicate anomaly identity")
    return sorted(items, key=lambda item: item["id"]), sources


def make_anomaly_audit(selection_path, pilot_path, output, *, dispositions=None, coarse=None, fine=None):
    from .runner import code_identity, runtime_identity
    items, sources = anomaly_inventory(selection_path, pilot_path, coarse=coarse, fine=fine)
    dispositions = dispositions or {}
    unknown = set(dispositions) - {item["id"] for item in items}
    if unknown:
        raise ValueError(f"dispositions do not match observed anomalies: {sorted(unknown)}")
    entries = []
    for item in items:
        supplied = dispositions.get(item["id"], {})
        disposition = supplied.get("disposition", "unresolved")
        explanation = supplied.get("explanation", "")
        paths = supplied.get("sources", [])
        if disposition not in ("unresolved", "explained", "resolved"):
            raise ValueError("unsupported anomaly disposition")
        if disposition != "unresolved" and (not isinstance(explanation, str) or not explanation.strip() or not paths):
            raise ValueError("resolved/explained anomalies need an explanation and retained trace sources")
        source_keys = []
        for n, path in enumerate(paths):
            key = f"anomaly:{item['id']}:{n}"
            _source(sources, key, path)
            source_keys.append(key)
        entries.append({"observed": item, "disposition": disposition, "explanation": explanation,
                        "sources": source_keys})
    _source(sources, "selection", selection_path)
    report = seal({"study_version": "0.3", "study_plan_sha256": PLAN_SHA256,
                   "status": "passed" if all(e["disposition"] != "unresolved" for e in entries) else "unresolved",
                   "code": code_identity(), "runtime": runtime_identity(), "entries": entries, "sources": sources,
                   "inputs": {"selection": str(Path(selection_path).resolve()),
                              "pilot": str(Path(pilot_path).resolve()) if pilot_path else None,
                              "coarse": str(Path(coarse).resolve()) if coarse else None,
                              "fine": str(Path(fine).resolve()) if fine else None}})
    write_new(output, report)
    return report


def verify_anomaly_audit(path, *, selection, pilot, coarse=None, fine=None):
    report = read_sealed(path)
    _current_identity(report)
    expected_inputs = {"selection": str(Path(selection).resolve()), "pilot": str(Path(pilot).resolve()) if pilot else None,
                       "coarse": str(Path(coarse).resolve()) if coarse else None,
                       "fine": str(Path(fine).resolve()) if fine else None}
    if report["inputs"] != expected_inputs or report["study_plan_sha256"] != PLAN_SHA256:
        raise ValueError("anomaly audit does not cover the qualification evidence")
    _verify_sources(report["sources"])
    items, sources = anomaly_inventory(selection, pilot, coarse=coarse, fine=fine)
    if [e["observed"] for e in report["entries"]] != items:
        raise ValueError("anomaly inventory changed or an observed anomaly was omitted")
    if any(report["sources"].get(k) != v for k, v in sources.items()):
        raise ValueError("anomaly audit lacks a retained execution source")
    for entry in report["entries"]:
        if entry["disposition"] not in ("unresolved", "explained", "resolved"):
            raise ValueError("unsupported anomaly disposition")
        if entry["disposition"] != "unresolved" and (not entry["explanation"].strip() or not entry["sources"]
                or any(key not in report["sources"] for key in entry["sources"])):
            raise ValueError("anomaly disposition lacks a retained explanation/trace")
    expected = "passed" if all(e["disposition"] != "unresolved" for e in report["entries"]) else "unresolved"
    if report["status"] != expected:
        raise ValueError("anomaly status does not match its dispositions")
    return report


def _test_checks(tests, code, runtime, sources):
    _source(sources, "tests", tests)
    metadata_path = Path(str(tests) + ".meta.json")
    _source(sources, "tests_meta", metadata_path)
    metadata = json.loads(metadata_path.read_text())
    if (metadata["code_sha256"] != code["sha256"] or metadata["runtime"] != runtime
            or metadata["report_sha256"] != sources["tests"]["sha256"]
            or metadata.get("test_source_sha256") != _test_identity()
            or metadata.get("unchanged_during_tests") is not True):
        raise ValueError("JUnit evidence does not bind the unchanged code/runtime/test version")
    checks = junit_checks(tests)
    observed = {}
    for case in ET.parse(tests).getroot().iter("testcase"):
        status = "failed" if case.find("failure") is not None or case.find("error") is not None else (
            "not_run" if case.find("skipped") is not None else "passed")
        observed.setdefault(case.get("name", ""), []).append(status)
    required = {name: observed.get(name, []) for name in V03_TESTS}
    checks["study_v03"] = {"status": "passed" if all(v == ["passed"] for v in required.values()) else "failed",
                            "tests": required, "sources": ["tests", "tests_meta"]}
    if metadata.get("returncode") != 0:
        for check in checks.values():
            check["status"] = "failed"
    return checks


def _numerical(protocol, coarse, fine, sources):
    from .runner import read_inputs
    cm, cr, ca = audit_cohort(coarse, sources, prefix="coarse")
    fm, fr, fa = audit_cohort(fine, sources, prefix="fine")
    verify_admission(coarse, cm)
    verify_admission(fine, fm)
    for label, folder in (("coarse", coarse), ("fine", fine)):
        admission_path = Path(folder) / "study_admission.json"
        admission = read_sealed(admission_path)
        _source(sources, label + ":admission", admission_path)
        _source(sources, label + ":pilot_evidence", admission["pilot_evidence_path"])
    expected_fine = copy.deepcopy(protocol)
    expected_fine["clock"]["dt_model_seconds"] = protocol["numerical_gates"]["refined_dt"]
    expected_fine["clock"]["steps"] = int(round(protocol["clock"]["horizon_model_seconds"] / expected_fine["clock"]["dt_model_seconds"]))
    if (cm["cohort"] != "numerical_coarse" or fm["cohort"] != "numerical_fine"
            or cm["protocol"] != protocol or fm["protocol"] != expected_fine or cm["seeds"] != fm["seeds"]):
        raise ValueError("numerical refinement changed the selected configuration or seed pairing")
    for seed in cm["seeds"]:
        for arm in ("U", "H"):
            u, v = [read_inputs(folder, seed, arm) for folder in (coarse, fine)]
            for field in ("pair_id", "geometry", "founders", "feed_schedule", "streams", "orientation_rad", "arm"):
                if u[field] != v[field]:
                    raise ValueError(f"numerical refinement changed physical input {field}")
    gates, clock = protocol["numerical_gates"], protocol["clock"]
    result = compare_timesteps(cr, fr, cm["seeds"], horizon_s=clock["horizon_model_seconds"],
                              coarse_dt_s=clock["dt_model_seconds"], fine_dt_s=gates["refined_dt"],
                              cohort=cm["cohort"], fine_cohort=fm["cohort"], tolerance=gates["effect_tolerance"],
                              confidence=gates["confidence_level"],
                              bootstrap_resamples=protocol["analysis"]["bootstrap_resamples"],
                              bootstrap_seed=protocol["analysis"]["bootstrap_seed"])
    result["integrity_audits"] = {"coarse": ca, "fine": fa}
    if ca["status"] != "passed" or fa["status"] != "passed":
        result["status"] = "failed"
    result["sources"] = [k for k in sources if k.startswith(("coarse:", "fine:"))]
    return result


def qualification_gate(report):
    checks = report.get("checks", {})
    required = (*REQUIRED_CHECKS, "study_v03", "controls", "selection", "anomalies", "focal_exposure")
    failed, pending = [], []
    if report.get("study_version") != "0.3" or report.get("study_plan_sha256") != PLAN_SHA256:
        failed.append("study plan mismatch")
    for name in required:
        check = checks.get(name, {})
        if check.get("status") == "failed":
            failed.append(name)
        elif check.get("status") != "passed":
            pending.append(name)
        if check and check != _bind(check):
            failed.append(name + " evidence digest")
    pilot = checks.get("pilot", {})
    if pilot.get("status") == "passed" and (pilot.get("seeds") != list(range(5000, 5020))
            or pilot.get("completed_pairs") != 20 or pilot.get("U_denominator") != 20
            or pilot.get("U_runs_with_both_population_births", 0) < 10
            or pilot.get("capacity_culling_count") != 0 or pilot.get("require_H_births") is not False
            or pilot.get("survival_count_gate") is not None):
        failed.append("fresh pilot rule mismatch")
    numerical = checks.get("numerical", {})
    if numerical.get("status") == "passed":
        n = numerical.get("n_pairs")
        if type(n) is not int or n not in (50, 200, 1000) or numerical.get("registered_seeds") != list(range(2000, 2000 + n)):
            failed.append("numerical prespecified look mismatch")
        if numerical.get("registered_look") is not True or numerical.get("Q_status") != "passed":
            failed.append("numerical uncertainty criterion")
        metadata = {"coarse_dt_s": .1, "fine_dt_s": .05, "horizon_s": 600., "tolerance": .05,
                    "confidence_family": .95, "confidence_per_look": 1 - (1 - .95) / 3,
                    "interval_method": "Bonferroni_Clopper_Pearson_category_simplex"}
        if any(numerical.get(key) != value for key, value in metadata.items()):
            failed.append("numerical prespecified interval rule mismatch")
        for entry in [numerical] + [numerical.get("arms", {}).get(arm, {}) for arm in ("U", "H")]:
            interval = entry.get("Q_interval", entry.get("survival_change_interval", []))
            if (len(interval) != 2 or not all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in interval)
                    or interval[0] > interval[1] or interval[0] < -.05 or interval[1] > .05):
                failed.append("numerical interval outside tolerance")
        if any(numerical.get("arms", {}).get(arm, {}).get("status") != "passed" for arm in ("U", "H")):
            failed.append("numerical arm criterion")
    exposure = checks.get("exposure", {})
    criterion = {"metric": "mean_paired_time_integrated_normalized_half_food_contrast", "threshold": .05, "operator": ">="}
    if exposure.get("status") == "passed" and (exposure.get("criterion") != criterion
            or exposure.get("criterion_sha256") != digest(criterion)
            or not math.isfinite(exposure.get("measured", math.nan)) or exposure["measured"] < .05):
        failed.append("resource geography criterion")
    status = "failed" if failed else "inconclusive" if pending else "passed"
    return {"status": status, "freeze_allowed": status == "passed", "reasons": failed, "pending": pending}


def _focal_exposure_check(records, sources, protocol):
    """Require focal-specific support/normalization, with no outcome threshold."""
    focal = next(s for s in protocol["species"] if s["role"] == "focal")
    horizon = protocol["clock"]["horizon_model_seconds"]
    failures, rows = [], []
    def require(condition):
        if not condition:
            raise ValueError("invalid focal exposure support or normalization")
    for record in records:
        try:
            exposure = record["exposure"]
            item = exposure["populations"][focal["name"]]
            early = record["diagnostics"]["early_exposure_contrast"]
            require(exposure["early_window_seconds"] == [0, 25])
            require(math.isclose(exposure["early_observed_seconds"], min(25, horizon), abs_tol=1e-8))
            require(math.isfinite(early) and -1 <= early <= 1)
            seconds, uptake, rates = item["cell_seconds"], item["uptake_pg"], item["uptake_pg_per_cell_second"]
            for values in (seconds, uptake):
                require(set(values) == {"rich", "poor", "total"})
                require(all(math.isfinite(v) and v >= 0 for v in values.values()))
                require(math.isclose(values["total"], values["rich"] + values["poor"], rel_tol=1e-9, abs_tol=1e-8))
            require(set(rates) == set(seconds))
            for region, value in seconds.items():
                require(rates[region] is None if value == 0 else math.isclose(rates[region], uptake[region] / value, rel_tol=1e-9, abs_tol=1e-8))
            require(item["rich_cell_time_fraction"] is None if seconds["total"] == 0 else
                    math.isclose(item["rich_cell_time_fraction"], seconds["rich"] / seconds["total"], rel_tol=1e-9, abs_tol=1e-8))
            require(item["living_count"] == record["focal_count"])
            require(item["initial_living_count"] == focal["init"])
            support = item["time_support"]
            require(math.isclose(support["observed_through_seconds"], horizon, abs_tol=1e-7))
            require(support["right_censored_at_observation_end"] is (record["focal_count"] > 0))
            extinction = support["first_extinction_time_seconds"]
            if record["extinction_time_s"] is None:
                require(extinction is None)
            else:
                require(math.isclose(extinction, record["extinction_time_s"], abs_tol=1e-7))
            for field in ("first_occupied_tick_start_seconds", "last_occupied_tick_end_seconds"):
                value = support[field]
                require(value is None if seconds["total"] == 0 else math.isfinite(value) and -1e-7 <= value <= horizon + 1e-7)
            require(item["residence_convention"] == "tick_start_position_and_abundance")
            require(item["uptake_convention"] == "event_location_before_movement")
            weighted = item["cell_time_weighted_field_contrast"]
            require(weighted is None if seconds["total"] == 0 else math.isfinite(weighted) and -1 <= weighted <= 1)
            rows.append({"seed": record["seed"], "arm": record["arm"], "early_exposure_contrast": early,
                         "early_observed_seconds": exposure["early_observed_seconds"], "focal": item})
        except (KeyError, TypeError, ValueError, AssertionError) as exc:
            failures.append({"seed": record["seed"], "arm": record["arm"], "reason": type(exc).__name__})
    return {"status": "passed" if not failures else "failed", "failures": failures,
            "threshold": None, "sources": sources, "per_arm": rows}


def _build_report(selection, pilot, tests, *, anomalies=None, coarse=None, fine=None):
    from .runner import code_identity, runtime_identity
    if (coarse is None) != (fine is None):
        raise ValueError("both numerical resolutions must be provided")
    selected = verify_selection(selection)
    protocol = selected_protocol(selected)
    sources = copy.deepcopy(selected["sources"])
    manifest, records, audited = audit_cohort(pilot, sources, prefix="pilot")
    verify_admission(pilot, manifest)
    if manifest["cohort"] != "pilot" or manifest["protocol"] != protocol:
        raise ValueError("fresh pilot does not use the configuration selected before H outcomes")
    _source(sources, "selection_report", selection)
    _source(sources, "pilot_admission", Path(pilot) / "study_admission.json")
    code, runtime = code_identity(), runtime_identity()
    checks = _test_checks(tests, code, runtime, sources)
    pilot_sources = [k for k in sources if k.startswith("pilot")]
    for name in ("food_ledger", "population_ledger", "disabled_mechanisms"):
        entry = checks.setdefault(name, {"status": "passed", "sources": []})
        entry["audit"] = audited["accounting"][name]
        entry["sources"].extend(pilot_sources)
        if not all(item["passed"] for item in entry["audit"]):
            entry["status"] = "failed"
    checks["food_ledger"]["pruning"] = audited["accounting"]["pruning"]
    if not all(item["passed"] for item in checks["food_ledger"]["pruning"]):
        checks["food_ledger"]["status"] = "failed"
    checks["matched_inputs"]["sources"].extend(pilot_sources)
    checks["pilot"] = {**fresh_pilot_decision(audited), "audit": audited, "sources": pilot_sources}
    indexed = {(r["seed"], r["arm"]): r for r in records}
    contrasts = [indexed[(s, "H")]["diagnostics"]["exposure_contrast"] - indexed[(s, "U")]["diagnostics"]["exposure_contrast"]
                 for s in manifest["seeds"]]
    measured = math.fsum(contrasts) / len(contrasts)
    criterion = {"metric": "mean_paired_time_integrated_normalized_half_food_contrast", "threshold": .05, "operator": ">="}
    checks["exposure"] = {"status": "passed" if math.isfinite(measured) and measured >= .05 else "failed",
                          "criterion": criterion, "criterion_sha256": digest(criterion), "measured": measured,
                          "per_pair": contrasts, "sources": pilot_sources}
    checks["focal_exposure"] = _focal_exposure_check(records, pilot_sources, protocol)
    index = verify_grid(selected["grid_path"])
    checks["controls"] = {"status": "passed", "controls_path": index["controls_path"],
                          "controls_report_sha256": index["controls_report_sha256"], "sources": ["grid_registration"]}
    checks["selection"] = {"status": "passed", "selected_candidate": selected["selected_candidate"],
                           "selection_sha256": selected["report_sha256"], "sources": ["selection_report", "grid_registration"]}
    checks["numerical"] = _numerical(protocol, coarse, fine, sources) if coarse else {
        "status": "not_run", "sources": [], "detail": "Numerical qualification absent; freeze blocked."}
    if anomalies:
        anomaly_report = verify_anomaly_audit(anomalies, selection=selection, pilot=pilot, coarse=coarse, fine=fine)
        _source(sources, "anomaly_report", anomalies)
        sources.update({"anomaly_source:" + k: v for k, v in anomaly_report["sources"].items()})
        checks["anomalies"] = {"status": anomaly_report["status"], "entries": anomaly_report["entries"],
                               "report_sha256": anomaly_report["report_sha256"],
                               "sources": ["anomaly_report"] + [k for k in sources if k.startswith("anomaly_source:")]}
    else:
        checks["anomalies"] = {"status": "not_run", "sources": [], "detail": "Retained anomaly inventory/disposition required."}
    checks = {name: _bind(check) for name, check in checks.items()}
    report = {"schema_version": 1, "study_version": "0.3", "study_plan_sha256": PLAN_SHA256,
              "protocol_sha256": digest(protocol), "protocol": protocol, "code": code,
              "code_sha256": code["sha256"], "runtime": runtime, "test_source_sha256": _test_identity(),
              "checks": checks, "sources": sources,
              "inputs": {"selection": str(Path(selection).resolve()), "pilot": str(Path(pilot).resolve()),
                         "tests": str(Path(tests).resolve()), "anomalies": str(Path(anomalies).resolve()) if anomalies else None,
                         "coarse": str(Path(coarse).resolve()) if coarse else None,
                         "fine": str(Path(fine).resolve()) if fine else None}}
    report["verdict"] = qualification_gate(report)
    return seal(report)


def qualify(selection, pilot, tests, output, *, anomalies=None, coarse=None, fine=None):
    if Path(output).exists():
        raise FileExistsError("qualification reports are append-only")
    report = _build_report(selection, pilot, tests, anomalies=anomalies, coarse=coarse, fine=fine)
    write_new(output, report)
    return report


def verify_report(report, *, rebuild=True):
    if report != seal(report):
        raise ValueError("qualification report hash mismatch")
    _current_identity(report)
    if report["test_source_sha256"] != _test_identity():
        raise ValueError("qualification tests changed")
    _verify_sources(report["sources"])
    if report["verdict"] != qualification_gate(report):
        raise ValueError("qualification verdict is not derived from the gates")
    if rebuild and report != _build_report(**report["inputs"]):
        raise ValueError("qualification evidence no longer reproduces its retained report")
    return report


def freeze(evidence, output):
    if Path(output).exists():
        raise FileExistsError("freezes are append-only")
    report = verify_report(read_sealed(evidence))
    if not qualification_gate(report)["freeze_allowed"]:
        raise ValueError("freeze blocked: " + json.dumps(report["verdict"], sort_keys=True))
    result = seal({"status": "qualified", "study_version": "0.3", "study_plan_sha256": PLAN_SHA256,
                   "protocol_sha256": report["protocol_sha256"], "code_sha256": report["code_sha256"],
                   "runtime": report["runtime"], "evidence": report}, "freeze_sha256")
    write_new(output, result)
    return result


def verify_freeze(path, *, protocol_sha256, code_sha256, runtime):
    record = copy.deepcopy(path) if isinstance(path, dict) else read_sealed(path, "freeze_sha256")
    if record != seal(record, "freeze_sha256"):
        raise ValueError("freeze hash mismatch")
    if (record.get("status") != "qualified" or record.get("study_version") != "0.3"
            or record.get("study_plan_sha256") != PLAN_SHA256 or record.get("protocol_sha256") != protocol_sha256
            or record.get("code_sha256") != code_sha256 or record.get("runtime") != runtime):
        raise ValueError("freeze does not match the requested study/configuration/code/runtime")
    report = verify_report(record["evidence"])
    if (report["protocol_sha256"] != protocol_sha256 or report["code_sha256"] != code_sha256
            or report["runtime"] != runtime or not qualification_gate(report)["freeze_allowed"]):
        raise ValueError("freeze evidence is unqualified or mismatched")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    qual = sub.add_parser("qualify")
    for name in ("selection", "pilot", "tests", "output"):
        qual.add_argument("--" + name, required=True)
    for name in ("anomalies", "coarse", "fine"):
        qual.add_argument("--" + name)
    audit = sub.add_parser("anomalies")
    audit.add_argument("--selection", required=True); audit.add_argument("--pilot"); audit.add_argument("--output", required=True)
    audit.add_argument("--coarse"); audit.add_argument("--fine"); audit.add_argument("--dispositions")
    frozen = sub.add_parser("freeze"); frozen.add_argument("--evidence", required=True); frozen.add_argument("--output", required=True)
    args = vars(parser.parse_args(argv))
    command = args.pop("command")
    if command == "qualify":
        result = qualify(**args)
    elif command == "freeze":
        result = freeze(**args)
    else:
        dispositions = args.pop("dispositions")
        result = make_anomaly_audit(args.pop("selection"), args.pop("pilot"), args.pop("output"),
                                    dispositions=json.loads(Path(dispositions).read_text()) if dispositions else None, **args)
    print(json.dumps(result.get("verdict", {k: result[k] for k in ("status",) if k in result}), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
