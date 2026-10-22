"""Build retained qualification evidence, then freeze only a qualified version.

Usage: python -m research.spatial_nutrient.qualification_run tests --output tests.xml
       python -m research.spatial_nutrient.qualification_run qualify --pilot PILOT \
           --tests tests.xml --output qualification.json [--coarse COARSE --fine FINE]
       python -m research.spatial_nutrient.qualification_run freeze \
           --evidence qualification.json --output freeze.json
"""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

from .analysis import analyze_pairs
from .qualification import (REQUIRED_CHECKS, canonical_sha256, compare_timesteps,
                            qualification_gate)


TEST_REQUIREMENTS = {
    "lifecycle": ["test_original_capacity_fails_and_amendment_passes",
                  "test_fed_isolated_cell_really_divides",
                  "test_variable_daughter_capacity_overflow_returns_to_parent"],
    "starvation": ["test_no_feed_cell_starves_and_population_reconciles"],
    "matched_inputs": ["test_between_arm_tampering_rejected_even_with_recomputed_hash",
                       "test_exact_feed_mass_times_and_physical_invariance"],
    "food_ledger": ["test_uptake_then_cull_remains_in_ledger"],
    "population_ledger": ["test_no_feed_cell_starves_and_population_reconciles",
                          "test_daughter_not_updated_on_birth_tick"],
    "scheduling": ["test_container_order_does_not_change_dynamics",
                   "test_arbitrary_neutral_labels_preserve_total_dynamics",
                   "test_scheduler_gives_both_populations_first_access"],
    "replay": ["test_duplicate_uniform_inputs_replay_exactly"],
    "null_allocator": ["test_zero_contrast_ecological_replay_is_exact"],
    "food_node_limit": ["test_node_limit_never_silently_deletes_food"],
    "failure_retry": ["test_failure_retained_and_retry_uses_same_saved_input"],
    "endpoint": ["test_last_tick_extinction_is_observed_survival_is_censored",
                 "test_completed_results_validate_physical_endpoint_and_analysis"],
    "analysis": ["test_hand_counted_table_preserves_pairing_and_zero_outcomes",
                 "test_degenerate_samples_use_nonzero_width_finite_sample_fallback",
                 "test_missing_duplicate_failed_and_extra_arms_are_not_dropped",
                 "test_zero_changes_at_50_is_inconclusive_not_converged",
                 "test_common_arm_error_cannot_hide_behind_zero_difference"],
}


def file_sha256(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def _test_identity():
    from .runner import ROOT
    files = sorted([*ROOT.joinpath("bacteria").glob("test_*.py"),
                    *ROOT.joinpath("research/spatial_nutrient").glob("test_*.py")])
    return canonical_sha256({str(p.relative_to(ROOT)): file_sha256(p) for p in files})


def run_tests(output):
    """Run the fixed qualification tests and bind the JUnit report to its code."""
    from .runner import ROOT, atomic_json, code_identity, runtime_identity
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    meta_path = Path(str(output) + ".meta.json")
    if output.exists() or meta_path.exists():
        raise FileExistsError("test evidence is append-only; choose a new output path")
    code, runtime, tests = code_identity(), runtime_identity(), _test_identity()
    command = [sys.executable, "-m", "pytest", "bacteria/test_experiment_support.py",
               "research/spatial_nutrient", "-q", f"--junitxml={output}"]
    result = subprocess.run(command, cwd=ROOT, check=False)
    unchanged = code == code_identity() and runtime == runtime_identity() and tests == _test_identity()
    metadata = {"code_sha256": code["sha256"], "runtime": runtime,
                "test_source_sha256": tests, "unchanged_during_tests": unchanged,
                "returncode": result.returncode, "command": command,
                "report_sha256": file_sha256(output) if output.exists() else None}
    atomic_json(meta_path, metadata)
    return metadata


def junit_checks(path):
    """Only observed, unique, passing test cases satisfy a named check."""
    path = Path(path)
    root = ET.parse(path).getroot()
    names = {}
    for case in root.iter("testcase"):
        status = ("failed" if case.find("failure") is not None or case.find("error") is not None
                  else "not_run" if case.find("skipped") is not None else "passed")
        names.setdefault(case.attrib.get("name", ""), []).append(status)
    checks = {}
    for name, requirements in TEST_REQUIREMENTS.items():
        observed = {test: names.get(test, []) for test in requirements}
        status = ("failed" if any("failed" in states or len(states) > 1 for states in observed.values())
                  else "passed" if all(states == ["passed"] for states in observed.values())
                  else "not_run")
        checks[name] = {"status": status, "tests": observed, "sources": ["tests", "tests_meta"]}
    return checks


def _bind(check):
    check = dict(check)
    check.pop("evidence_sha256", None)
    check["evidence_sha256"] = canonical_sha256(check)
    return check


def _source(sources, key, path):
    path = Path(path).resolve()
    sources[key] = {"path": str(path), "sha256": file_sha256(path)}


def _audit_pilot(output, manifest, records, sources):
    """Reconcile retained event streams against final counts and external mass."""
    from .runner import artifact_path, digest, input_path, open_jsonl, read_inputs
    from .inputs import assert_matched_arms
    output = Path(output)
    protocol = manifest["protocol"]
    initial = {sp["name"]: sp["init"] for sp in protocol["species"]}
    tolerance = protocol["execution"]["mass_tolerance"]
    prune_limit = protocol["numerical_gates"]["max_pruned_fraction_of_sources"]
    results = {"food_ledger": [], "population_ledger": [], "disabled_mechanisms": [],
               "matched_inputs": [], "pruning": [], "births_by_arm": {"U": 0, "H": 0}}
    for pair in manifest["pairs"]:
        u = read_inputs(output, pair["seed"], "U")
        h = read_inputs(output, pair["seed"], "H")
        assert_matched_arms(u, h)
        for arm, inputs in (("U", u), ("H", h)):
            if digest(inputs) != pair["inputs"][arm]:
                raise ValueError("pilot input no longer matches registration")
            _source(sources, f"pilot:{pair['seed']}:{arm}:inputs",
                    input_path(output, pair["seed"], arm))
    for record in records:
        prefix = f"pilot:{record['seed']}:{record['arm']}"
        arm_dir = output / "runs" / str(record["seed"]) / record["arm"]
        attempt = arm_dir / f"attempt-{record['attempt']:03d}"
        for name, path in (("result", attempt / "result.json"), ("events", artifact_path(attempt, "events.jsonl")),
                           ("trajectory", artifact_path(attempt, "trajectory.jsonl")),
                           ("completion", arm_dir / "completed.json")):
            _source(sources, prefix + ":" + name, path)
        births, deaths, transfer = Counter(), Counter(), Counter()
        forbidden = []
        with open_jsonl(artifact_path(attempt, "events.jsonl")) as handle:
            for line in handle:
                event = json.loads(line)
                kind = event["kind"]
                if kind in ("birth", "death"):
                    (births if kind == "birth" else deaths)[event["population"]] += 1
                if kind == "death" and event["cause"] not in ("starvation", "capacity"):
                    forbidden.append(event["cause"])
                if "amount_pg" in event:
                    amount = float(event["amount_pg"])
                    if not math.isfinite(amount) or amount < 0:
                        raise ValueError("invalid retained transfer")
                    transfer[kind] += amount
        last = None
        with open_jsonl(artifact_path(attempt, "trajectory.jsonl")) as handle:
            for line in handle:
                last = json.loads(line)
        if last is None or abs(last["time_s"] - protocol["clock"]["horizon_model_seconds"]) > 1e-9:
            raise ValueError("pilot trajectory lacks its physical endpoint")
        ledger = record["ledger"]
        for kind, key in (("feed", "feed_pg"), ("recycling", "recycled_pg"),
                          ("uptake", "uptake_pg"), ("washout", "washout_pg"), ("prune", "pruned_pg")):
            if abs(transfer[kind] - ledger[key]) > tolerance * max(1., transfer[kind]):
                raise ValueError(f"retained event ledger mismatch: {prefix}/{key}")
        mass = last["tick"]["mass_end_pg"]
        predicted = ledger["feed_pg"] + ledger["recycled_pg"] - ledger["uptake_pg"] - ledger["washout_pg"] - ledger["pruned_pg"]
        scale = max(1., ledger["feed_pg"] + ledger["recycled_pg"])
        residual = mass - predicted
        results["food_ledger"].append({"arm": prefix, "residual_pg": residual,
                                      "passed": math.isfinite(residual) and abs(residual) <= tolerance * scale})
        passed = all(record["counts"].get(pop) == count + births[pop] - deaths[pop]
                     for pop, count in initial.items())
        passed = passed and sum(births.values()) == ledger["births"] and sum(deaths.values()) == ledger["deaths"]
        results["population_ledger"].append({"arm": prefix, "passed": passed})
        results["disabled_mechanisms"].append({"arm": prefix, "forbidden_death_causes": forbidden,
                                               "passed": not forbidden})
        fraction = ledger["pruned_pg"] / scale
        results["pruning"].append({"arm": prefix, "fraction_of_sources": fraction,
                                   "passed": fraction <= prune_limit})
        results["births_by_arm"][record["arm"]] += ledger["births"]
    return results


def qualify(pilot, tests, output, *, coarse=None, fine=None):
    from .runner import (atomic_json, code_identity, collect, digest, input_path,
                         read_inputs, runtime_identity)
    if (coarse is None) != (fine is None):
        raise ValueError("provide both coarse and fine cohorts, or neither")
    output = Path(output)
    if output.exists():
        raise FileExistsError("qualification reports are append-only")
    manifest, records = collect(pilot)
    code, runtime = code_identity(), runtime_identity()
    if manifest["code"]["sha256"] != code["sha256"] or manifest["runtime"] != runtime:
        raise ValueError("pilot code/runtime does not match the qualification candidate")
    if manifest["cohort"] != "pilot" or manifest["seeds"] != list(range(20)):
        raise ValueError("qualification requires the complete registered pilot")
    protocol = manifest["protocol"]
    p_clock = protocol["clock"]
    analysis = analyze_pairs(records, manifest["seeds"], cohort="pilot",
                             horizon_s=p_clock["horizon_model_seconds"], dt_s=p_clock["dt_model_seconds"],
                             bootstrap_resamples=protocol["analysis"]["bootstrap_resamples"],
                             bootstrap_seed=protocol["analysis"]["bootstrap_seed"])
    sources = {}
    _source(sources, "pilot_manifest", Path(pilot) / "manifest.json")
    _source(sources, "tests", tests)
    meta_path = Path(str(tests) + ".meta.json")
    _source(sources, "tests_meta", meta_path)
    metadata = json.loads(meta_path.read_text())
    if (metadata["code_sha256"] != code["sha256"] or metadata["runtime"] != runtime
            or metadata["report_sha256"] != sources["tests"]["sha256"]
            or metadata.get("test_source_sha256") != _test_identity()
            or metadata.get("unchanged_during_tests") is not True):
        raise ValueError("JUnit evidence was not produced by this unchanged code/runtime/test version")
    checks = junit_checks(tests)
    if metadata.get("returncode") != 0:
        for entry in checks.values():
            entry["status"] = "failed"
            entry["detail"] = "qualification test process failed"
    audited = _audit_pilot(pilot, manifest, records, sources)
    audited_sources = [key for key in sources if key.startswith("pilot")]
    for name in ("food_ledger", "population_ledger", "disabled_mechanisms"):
        entry = checks.setdefault(name, {"status": "passed", "sources": []})
        entry["pilot_audit"] = audited[name]
        entry["sources"].extend(audited_sources)
        if not all(item["passed"] for item in audited[name]):
            entry["status"] = "failed"
    checks["food_ledger"]["pruning_audit"] = audited["pruning"]
    if not all(item["passed"] for item in audited["pruning"]):
        checks["food_ledger"]["status"] = "failed"
    checks["matched_inputs"]["sources"].extend(audited_sources)
    checks["pilot"] = {"status": "passed" if all(v > 0 for v in audited["births_by_arm"].values()) else "inconclusive",
                        "seeds": manifest["seeds"], "completed_pairs": analysis["n_pairs"],
                        "births_by_arm": audited["births_by_arm"], "sources": audited_sources,
                        "detail": "births must occur in each arm cohort to exercise population recovery"}
    indexed = {(r["seed"], r["arm"]): r for r in records}
    contrasts = [indexed[(s, "H")]["diagnostics"]["exposure_contrast"] -
                 indexed[(s, "U")]["diagnostics"]["exposure_contrast"] for s in manifest["seeds"]]
    measured = math.fsum(contrasts) / len(contrasts)
    criterion = {"metric": "mean_paired_time_integrated_normalized_half_food_contrast",
                 "threshold": protocol["numerical_gates"]["min_realized_food_contrast"], "operator": ">="}
    checks["exposure"] = {"status": "passed" if math.isfinite(measured) and measured >= criterion["threshold"] else "failed",
                          "criterion": criterion, "criterion_sha256": canonical_sha256(criterion),
                          "measured": measured, "per_pair": contrasts, "sources": audited_sources,
                          "pool": "all available food, including recycled material",
                          "normalization": "(rich-poor)/(rich+poor); zero when no food, averaged over physical time"}
    if coarse is None:
        checks["numerical"] = {"status": "not_run", "sources": [],
                               "detail": "coarse and refined cohorts not supplied; freeze blocked"}
    else:
        cm, cr = collect(coarse)
        fm, fr = collect(fine)
        if cm["cohort"] != "numerical_coarse" or fm["cohort"] != "numerical_fine":
            raise ValueError("numerical cohorts must retain their separate coarse/fine labels")
        expected_fine = copy.deepcopy(protocol)
        expected_fine["clock"]["dt_model_seconds"] = protocol["numerical_gates"]["refined_dt"]
        expected_fine["clock"]["steps"] = int(round(p_clock["horizon_model_seconds"] / expected_fine["clock"]["dt_model_seconds"]))
        if cm["protocol"] != protocol or fm["protocol"] != expected_fine:
            raise ValueError("numerical protocols must differ only by the registered dt and step count")
        for label, m, folder in (("coarse", cm, coarse), ("fine", fm, fine)):
            if m["code"]["sha256"] != code["sha256"] or m["runtime"] != runtime:
                raise ValueError("numerical code/runtime mismatch")
            _source(sources, label + "_manifest", Path(folder) / "manifest.json")
        if cm["seeds"] != fm["seeds"]:
            raise ValueError("numerical cohorts contain different registered seeds")
        numerical = compare_timesteps(cr, fr, cm["seeds"],
                    horizon_s=p_clock["horizon_model_seconds"], coarse_dt_s=p_clock["dt_model_seconds"],
                    fine_dt_s=expected_fine["clock"]["dt_model_seconds"],
                    cohort=cm["cohort"], fine_cohort=fm["cohort"],
                    tolerance=protocol["numerical_gates"]["effect_tolerance"],
                    confidence=protocol["numerical_gates"]["confidence_level"],
                    bootstrap_resamples=protocol["analysis"]["bootstrap_resamples"],
                    bootstrap_seed=protocol["analysis"]["bootstrap_seed"])
        numeric_sources = ["coarse_manifest", "fine_manifest"]
        for label, items, folder in (("coarse", cr, coarse), ("fine", fr, fine)):
            for record in items:
                key = f"{label}:{record['seed']}:{record['arm']}:result"
                _source(sources, key, Path(folder) / "runs" / str(record["seed"]) / record["arm"] /
                        f"attempt-{record['attempt']:03d}" / "result.json")
                numeric_sources.append(key)
        # Changing dt must not change founders, external source timing/positions,
        # or the named initial RNG streams. Save their actual bytes as evidence.
        for seed in cm["seeds"]:
            for arm in ("U", "H"):
                pair_inputs = []
                for label, folder, m in (("coarse", coarse, cm), ("fine", fine, fm)):
                    path = input_path(folder, seed, arm)
                    item = read_inputs(folder, seed, arm)
                    registered = next(pair for pair in m["pairs"] if pair["seed"] == seed)
                    if digest(item) != registered["inputs"][arm]:
                        raise ValueError("numerical input no longer matches its registration")
                    key = f"{label}:{seed}:{arm}:inputs"
                    _source(sources, key, path)
                    numeric_sources.append(key)
                    pair_inputs.append(item)
                for field in ("pair_id", "geometry", "founders", "feed_schedule", "streams", "orientation_rad", "arm"):
                    if pair_inputs[0][field] != pair_inputs[1][field]:
                        raise ValueError(f"numerical refinement changed {field}")
        checks["numerical"] = {**numerical, "sources": numeric_sources}
    for name in REQUIRED_CHECKS:
        checks[name] = _bind(checks.get(name, {"status": "not_run", "sources": []}))
    report = {"schema_version": 1, "protocol_sha256": digest(protocol),
              "code_sha256": code["sha256"], "runtime": runtime,
              "test_source_sha256": _test_identity(), "checks": checks, "sources": sources}
    report["verdict"] = qualification_gate(report, protocol_sha256=digest(protocol), code_sha256=code["sha256"])
    report["report_sha256"] = canonical_sha256(report)
    atomic_json(output, report)
    return report


def verify_report(report):
    report = copy.deepcopy(report)
    expected = report.pop("report_sha256", None)
    if canonical_sha256(report) != expected:
        raise ValueError("qualification report hash mismatch")
    for entry in report["checks"].values():
        if entry != _bind(entry):
            raise ValueError("qualification check hash mismatch")
        if any(key not in report["sources"] for key in entry.get("sources", [])):
            raise ValueError("qualification references missing source")
    for source in report["sources"].values():
        if file_sha256(source["path"]) != source["sha256"]:
            raise ValueError("qualification source changed")
    report["report_sha256"] = expected
    return report


def freeze(evidence, output):
    from .runner import atomic_json, code_identity, runtime_identity
    output = Path(output)
    if output.exists():
        raise FileExistsError("freeze files are append-only")
    report = verify_report(json.loads(Path(evidence).read_text()))
    if (report["code_sha256"] != code_identity()["sha256"] or report["runtime"] != runtime_identity()
            or report["test_source_sha256"] != _test_identity()):
        raise ValueError("code/runtime/tests changed after qualification")
    verdict = qualification_gate(report, protocol_sha256=report["protocol_sha256"], code_sha256=report["code_sha256"])
    if not verdict["freeze_allowed"]:
        raise ValueError(f"qualification is {verdict['status']}; main freeze blocked")
    result = {"status": "qualified", "protocol_sha256": report["protocol_sha256"],
              "code_sha256": report["code_sha256"], "runtime": report["runtime"], "evidence": report}
    result["freeze_sha256"] = canonical_sha256(result)
    atomic_json(output, result)
    return result


def verify_freeze(path, *, protocol_sha256, code_sha256, runtime):
    """Runner entry point: recheck retained evidence before admitting main seeds."""
    frozen = copy.deepcopy(path) if isinstance(path, dict) else json.loads(Path(path).read_text())
    expected = frozen.pop("freeze_sha256", None)
    if canonical_sha256(frozen) != expected:
        raise ValueError("freeze hash mismatch")
    if frozen.get("status") != "qualified" or frozen.get("protocol_sha256") != protocol_sha256:
        raise ValueError("freeze does not qualify this protocol")
    if frozen.get("code_sha256") != code_sha256 or frozen.get("runtime") != runtime:
        raise ValueError("code/runtime differs from freeze")
    report = verify_report(frozen["evidence"])
    if report.get("test_source_sha256") != _test_identity():
        raise ValueError("qualification tests changed since freeze")
    verdict = qualification_gate(report, protocol_sha256=protocol_sha256, code_sha256=code_sha256)
    if not verdict["freeze_allowed"]:
        raise ValueError("freeze evidence does not pass all qualification gates")
    frozen["freeze_sha256"] = expected
    return frozen


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    test = commands.add_parser("tests"); test.add_argument("--output", required=True)
    qual = commands.add_parser("qualify")
    for name in ("pilot", "tests", "output"):
        qual.add_argument("--" + name, required=True)
    qual.add_argument("--coarse"); qual.add_argument("--fine")
    frozen = commands.add_parser("freeze")
    frozen.add_argument("--evidence", required=True); frozen.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    if args.command == "tests":
        result = run_tests(args.output)
        return 0 if result["returncode"] == 0 and result["unchanged_during_tests"] else 1
    if args.command == "qualify":
        report = qualify(args.pilot, args.tests, args.output, coarse=args.coarse, fine=args.fine)
        print(json.dumps(report["verdict"], indent=2))
        return 0 if report["verdict"]["freeze_allowed"] else 2
    result = freeze(args.evidence, args.output)
    print(json.dumps({"status": result["status"], "freeze_sha256": result["freeze_sha256"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
