"""Prospective v0.3 admission, complete-grid selection, and evidence audits.

The v0.3 study plan is not a new simulation schema.  This module binds that
prespecified plan to strict v0.2 execution inputs, without changing its gates.
All settings and decisions are Own-model choices; see PROTOCOL_V03.md.
"""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import json
import math
from pathlib import Path

from .protocol import digest, load_protocol, validate_protocol
from .qualification_run import _audit_pilot, _source, file_sha256

PLAN_PATH = Path(__file__).with_name("study_plan_v03.json")
PLAN_SHA256 = "c72ff2aae74ff715f4eeb4c740ad0755f976b5eff1d29216f1aebc2d983cd975"
GRID_INDEX = "grid_registration.json"


def seal(record, key="report_sha256"):
    result = copy.deepcopy(record)
    result.pop(key, None)
    result[key] = digest(result)
    return result


def read_sealed(path, key="report_sha256"):
    result = json.loads(Path(path).read_text())
    if result != seal(result, key):
        raise ValueError(f"{key} mismatch: {path}")
    return result


def write_new(path, value):
    from .runner import atomic_json
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"retained evidence is append-only: {path}")
    atomic_json(path, value)


def load_plan(path=None):
    """Reject weakened gates, unknown settings, and edits to registered inputs."""
    path = Path(path or PLAN_PATH).resolve()
    plan = json.loads(path.read_text())
    if digest(plan) != PLAN_SHA256:
        raise ValueError("v0.3 study plan differs from the prespecified canonical plan; amendment required")
    base = load_protocol(path.parent / plan["base_configuration"]["path"])
    if digest(base) != plan["base_configuration"]["canonical_sha256"]:
        raise ValueError("v0.3 base configuration changed")
    for candidate in plan["selection"]["candidates"]:
        protocol = load_protocol(path.parent / candidate["config_path"])
        if digest(protocol) != candidate["config_canonical_sha256"]:
            raise ValueError(f"prespecified candidate changed: {candidate['id']}")
    return plan


def _current_identity(record):
    from .runner import code_identity, runtime_identity
    if record["code"]["sha256"] != code_identity()["sha256"] or record["runtime"] != runtime_identity():
        raise ValueError("evidence code/runtime does not match the current candidate")


def _verify_sources(sources):
    for source in sources.values():
        if file_sha256(source["path"]) != source["sha256"]:
            raise ValueError(f"retained source changed: {source['path']}")


def _controls(path, plan):
    from .controls_v03 import verify_controls
    from .runner import code_identity, runtime_identity
    report = verify_controls(path, code_sha256=code_identity()["sha256"],
                             runtime=runtime_identity(), plan_sha256=digest(plan))
    if report["status"] != "passed":
        raise ValueError("isolated lifecycle controls did not pass; grid blocked")
    return report


def register_grid(output, controls_path, plan_path=None):
    """Register all eight complete cohorts before admitting the first arm."""
    from .runner import code_identity, register, runtime_identity
    output = Path(output).resolve()
    plan_path = Path(plan_path or PLAN_PATH).resolve()
    plan = load_plan(plan_path)
    controls = _controls(controls_path, plan)
    output.mkdir(parents=True, exist_ok=False)
    code, runtime = code_identity(), runtime_identity()
    entries = []
    for candidate in plan["selection"]["candidates"]:
        directory = output / candidate["id"]
        protocol = load_protocol(plan_path.parent / candidate["config_path"])
        manifest = register(directory, protocol, plan["selection"]["cohort_prefix"] + candidate["id"],
                            seeds=plan["selection"]["seeds"])
        entries.append({"id": candidate["id"], "rank": candidate["rank"],
                        "path": str(directory), "manifest_sha256": manifest["manifest_sha256"]})
    if code != code_identity() or runtime != runtime_identity():
        raise ValueError("code/runtime changed during whole-grid registration")
    index = seal({"study_version": "0.3", "status": "registered", "plan_path": str(plan_path),
                  "study_plan_sha256": digest(plan), "controls_path": str(Path(controls_path).resolve()),
                  "controls_report_sha256": controls["report_sha256"], "code": code, "runtime": runtime,
                  "seeds": plan["selection"]["seeds"], "candidates": entries}, "registration_sha256")
    write_new(output / GRID_INDEX, index)
    return verify_grid(output)


def verify_grid(output, *, verify_control_evidence=True):
    from .runner import read_manifest
    output = Path(output).resolve()
    index = read_sealed(output / GRID_INDEX, "registration_sha256")
    plan = load_plan(index["plan_path"])
    _current_identity(index)
    if index["study_version"] != "0.3" or index["study_plan_sha256"] != digest(plan):
        raise ValueError("grid plan identity mismatch")
    expected = plan["selection"]["candidates"]
    if index["seeds"] != plan["selection"]["seeds"] or len(index["candidates"]) != len(expected):
        raise ValueError("all eight prespecified 20-pair candidates must be registered")
    for entry, candidate in zip(index["candidates"], expected):
        if (entry["id"], entry["rank"], entry["path"]) != (
                candidate["id"], candidate["rank"], str(output / candidate["id"])):
            raise ValueError("grid candidate order/path mismatch")
        manifest = read_manifest(entry["path"])
        if (manifest["manifest_sha256"] != entry["manifest_sha256"]
                or manifest["protocol_sha256"] != candidate["config_canonical_sha256"]
                or manifest["seeds"] != index["seeds"]
                or manifest["cohort"] != plan["selection"]["cohort_prefix"] + candidate["id"]
                or manifest["code"] != index["code"] or manifest["runtime"] != index["runtime"]):
            raise ValueError("grid manifest differs from prospective whole-grid registration")
    if verify_control_evidence:
        controls = _controls(index["controls_path"], plan)
        if controls["report_sha256"] != index["controls_report_sha256"]:
            raise ValueError("controls changed after grid registration")
    return index


def verify_admission(output, manifest):
    """Runner hook: missing complete-grid or selection admission fails closed."""
    output = Path(output).resolve()
    if manifest["cohort"].startswith("diagnostic_v03_"):
        index = verify_grid(output.parent, verify_control_evidence=False)
        entry = next((c for c in index["candidates"] if c["path"] == str(output)), None)
        if entry is None or entry["manifest_sha256"] != manifest["manifest_sha256"]:
            raise ValueError("arm was not prospectively admitted in the complete v0.3 grid")
    elif any(a["change"].startswith("study_v03:") for a in manifest["protocol"]["amendments"]):
        if manifest["cohort"] == "main":
            return  # runner's v0.3 freeze verification is the main admission gate.
        admission = read_sealed(output / "study_admission.json")
        selection = read_sealed(admission["selection_path"])
        if (selection["status"] != "selected" or admission["selection_sha256"] != selection["report_sha256"]
                or admission["manifest_sha256"] != manifest["manifest_sha256"]):
            raise ValueError("selected-configuration admission mismatch")
        _current_identity(selection)
        expected = selected_protocol(selection, refined=manifest["cohort"] == "numerical_fine")
        if manifest["protocol"] != expected:
            raise ValueError("registered protocol differs from selection recorded before treatment")
        _validate_selected_seeds(manifest["cohort"], manifest["seeds"])
        if manifest["cohort"].startswith("numerical_"):
            _verify_pilot_admission(admission.get("pilot_evidence_path"), selection, full=False)
            if read_sealed(admission["pilot_evidence_path"])["report_sha256"] != admission.get("pilot_evidence_sha256"):
                raise ValueError("pilot qualification changed after numerical admission")


def audit_cohort(output, sources, *, prefix, uniform_replay=False):
    """Audit every completed arm, all raw outputs and retained failure attempts."""
    from .analysis import validate_records
    from .runner import artifact_path, collect, open_jsonl
    output = Path(output).resolve()
    manifest, records = collect(output)
    _current_identity(manifest)
    p = manifest["protocol"]
    validate_records(records, manifest["seeds"], cohort=manifest["cohort"],
                     horizon_s=p["clock"]["horizon_model_seconds"], dt_s=p["clock"]["dt_model_seconds"])
    local_sources = {}
    audited = _audit_pilot(output, manifest, records, local_sources)
    for key, value in local_sources.items():
        sources[f"{prefix}:{key}"] = value
    _source(sources, prefix + ":manifest", output / "manifest.json")
    populations = {s["role"]: s["name"] for s in p["species"]}
    by_arm, attempts, exposure_trace = [], [], []
    for record in records:
        arm_dir = output / "runs" / str(record["seed"]) / record["arm"]
        attempt = arm_dir / f"attempt-{record['attempt']:03d}"
        _source(sources, f"{prefix}:{record['seed']}:{record['arm']}:rng_final", attempt / "rng_final.json")
        births, causes, uptake = Counter(), Counter(), Counter()
        with open_jsonl(artifact_path(attempt, "events.jsonl")) as handle:
            for line in handle:
                event = json.loads(line)
                if event["kind"] == "birth":
                    births[event["population"]] += 1
                elif event["kind"] == "death":
                    causes[event["cause"]] += 1
                elif event["kind"] == "uptake":
                    uptake[(event["population"], event["region"])] += event["amount_pg"]
        with open_jsonl(artifact_path(attempt, "trajectory.jsonl")) as handle:
            last = None
            for line in handle:
                last = json.loads(line)
        trace_ok = False
        if last is not None and "exposure" in record:
            summaries = record["exposure"]["populations"]
            trace_ok = (summaries == last["tick"].get("population_exposure")
                        and record["diagnostics"].get("early_exposure_contrast") == last["tick"].get("early_exposure_contrast")
                        and record["exposure"].get("early_observed_seconds") == last["tick"].get("early_exposure_time_seconds"))
            for population in populations.values():
                for region in ("rich", "poor"):
                    observed = summaries.get(population, {}).get("uptake_pg", {}).get(region, math.nan)
                    traced = uptake[(population, region)]
                    trace_ok = trace_ok and math.isclose(observed, traced, rel_tol=p["execution"]["mass_tolerance"], abs_tol=1e-8)
        exposure_trace.append({"seed": record["seed"], "arm": record["arm"], "passed": trace_ok})
        by_arm.append({"seed": record["seed"], "arm": record["arm"], "focal_count": record["focal_count"],
                       "births": {role: births[name] for role, name in populations.items()},
                       "capacity_culling_count": causes["capacity"]})
        for directory in sorted(arm_dir.glob("attempt-*")):
            # Failed and interrupted attempts are evidence too; retries cannot erase them.
            for path in sorted(directory.iterdir()):
                if path.is_file():
                    _source(sources, f"{prefix}:attempt:{record['seed']}:{record['arm']}:{directory.name}:{path.name}", path)
            if directory != attempt:
                result_path = directory / "result.json"
                retained = json.loads(result_path.read_text()) if result_path.exists() else {}
                attempts.append({"id": f"{prefix}:attempt:{record['seed']}:{record['arm']}:{directory.name}",
                                 "status": retained.get("status", "interrupted"),
                                 "error_type": retained.get("error_type"), "error": retained.get("error")})
    passed = all(item["passed"] for name in ("food_ledger", "population_ledger", "disabled_mechanisms", "pruning")
                 for item in audited[name])
    passed = passed and all(item["passed"] for item in exposure_trace)
    replay = []
    if uniform_replay:
        indexed = {(r["seed"], r["arm"]): r for r in records}
        for seed in manifest["seeds"]:
            u, h = indexed[(seed, "U")], indexed[(seed, "H")]
            replay.append({"seed": seed, "passed": u["outputs"] == h["outputs"] and all(
                u[k] == h[k] for k in ("focal_count", "counts", "ledger", "diagnostics", "extinction_time_s", "established"))})
        passed = passed and all(r["passed"] for r in replay)
    return manifest, records, {"status": "passed" if passed else "failed", "completed_pairs": len(records) // 2,
                               "seeds": manifest["seeds"], "by_arm": by_arm, "accounting": audited,
                               "uniform_replay": replay, "exposure_trace": exposure_trace, "retained_attempts": attempts}


def candidate_decision(audit):
    """Only the 20 U runs enter the ecological eligibility denominator."""
    u = [r for r in audit["by_arm"] if r["arm"] == "U"]
    h = [r for r in audit["by_arm"] if r["arm"] == "H"]
    seeds = list(range(4000, 4020))
    complete = (audit["seeds"] == seeds and audit["completed_pairs"] == 20
                and sorted(r["seed"] for r in u) == seeds and sorted(r["seed"] for r in h) == seeds)
    birth_runs = sum(r["births"]["focal"] > 0 and r["births"]["resident"] > 0 for r in u)
    survivors = sum(r["focal_count"] > 0 for r in u)
    culls = sum(r["capacity_culling_count"] for r in audit["by_arm"])
    criteria = {"complete_pairs": complete, "matching_replay_accounting": audit["status"] == "passed",
                "no_capacity_culling": culls == 0, "U_both_population_births": birth_runs >= 10,
                "U_survival_range": 4 <= survivors <= 16}
    return {"eligible": all(criteria.values()), "criteria": criteria, "U_denominator": len(u),
            "U_runs_with_both_population_births": birth_runs, "U_survival_count": survivors,
            "capacity_culling_count": culls}


def select_candidate(grid, output):
    index = verify_grid(grid)
    sources, candidates = {}, []
    _source(sources, "grid_registration", Path(grid) / GRID_INDEX)
    _source(sources, "study_plan", index["plan_path"])
    for entry in index["candidates"]:
        manifest, records, audit = audit_cohort(entry["path"], sources, prefix=entry["id"], uniform_replay=True)
        candidates.append({"id": entry["id"], "rank": entry["rank"], "path": entry["path"],
                           "manifest_sha256": manifest["manifest_sha256"], "audit": audit,
                           **candidate_decision(audit)})
    selected = next((c["id"] for c in candidates if c["eligible"]), None)
    result = seal({"study_version": "0.3", "status": "selected" if selected else "blocked",
                   "study_plan_sha256": PLAN_SHA256, "grid_path": str(Path(grid).resolve()),
                   "grid_registration_sha256": index["registration_sha256"], "code": index["code"],
                   "runtime": index["runtime"], "candidates": candidates, "selected_candidate": selected,
                   "sources": sources, "decision": "first eligible by prespecified rank" if selected else
                   "No eligible candidate; stop. A new grid requires an amendment and fresh seeds."})
    write_new(output, result)
    return result


def verify_selection(path, *, full=True):
    result = read_sealed(path)
    _current_identity(result)
    index = verify_grid(result["grid_path"], verify_control_evidence=full)
    if result["study_plan_sha256"] != PLAN_SHA256 or result["grid_registration_sha256"] != index["registration_sha256"]:
        raise ValueError("selection does not bind the complete prospective grid")
    if full:
        _verify_sources(result["sources"])
    if [c["id"] for c in result["candidates"]] != [c["id"] for c in index["candidates"]]:
        raise ValueError("selection omits a prespecified candidate")
    for c, entry in zip(result["candidates"], index["candidates"]):
        if (c["rank"], c["path"], c["manifest_sha256"]) != (entry["rank"], entry["path"], entry["manifest_sha256"]):
            raise ValueError("selection candidate provenance mismatch")
        if any(c[key] != value for key, value in candidate_decision(c["audit"]).items()):
            raise ValueError("selection eligibility was not derived from the fixed U-only rule")
        if full:
            observed_sources = {}
            _, _, observed_audit = audit_cohort(entry["path"], observed_sources, prefix=entry["id"], uniform_replay=True)
            if observed_audit != c["audit"] or any(result["sources"].get(k) != v for k, v in observed_sources.items()):
                raise ValueError("selection audit does not reproduce the retained candidate outcomes")
    selected = next((c["id"] for c in result["candidates"] if c["eligible"]), None)
    if result["selected_candidate"] != selected or result["status"] != ("selected" if selected else "blocked"):
        raise ValueError("selection is not the first eligible candidate")
    return result


def selected_protocol(selection, *, refined=False):
    if selection["status"] != "selected":
        raise ValueError("no eligible candidate: fresh pilot and freeze blocked")
    plan = load_plan()
    candidate = next(c for c in plan["selection"]["candidates"] if c["id"] == selection["selected_candidate"])
    protocol = load_protocol(PLAN_PATH.parent / candidate["config_path"])
    protocol["feeding"]["contrast_nodes"] = 1
    protocol["replication"]["pilot_seeds"] = plan["fresh_pilot"]["seeds"]
    protocol["amendments"].append({"from_version": "0.2", "change": f"study_v03:{candidate['id']}; selected before fresh-pilot H inputs; contrast_nodes=1; pilot seeds 5000–5019.",
                                   "reason": f"Prespecified study plan {PLAN_SHA256}; selection evidence {selection['report_sha256']}"})
    if refined:
        protocol["clock"]["dt_model_seconds"] = protocol["numerical_gates"]["refined_dt"]
        protocol["clock"]["steps"] = int(round(protocol["clock"]["horizon_model_seconds"] / protocol["clock"]["dt_model_seconds"]))
    return validate_protocol(protocol)


def _validate_selected_seeds(cohort, seeds):
    if cohort == "pilot":
        if seeds != list(range(5000, 5020)):
            raise ValueError("fresh v0.3 pilot requires seeds 5000–5019")
    elif cohort in ("numerical_coarse", "numerical_fine"):
        if len(seeds) not in (50, 200, 1000) or seeds != list(range(2000, 2000 + len(seeds))):
            raise ValueError("numerical cohort must be a declared cumulative 50/200/1000-pair look")
    else:
        raise ValueError("selected registration supports fresh pilot and numerical qualification only")


def _verify_pilot_admission(path, selection, *, full=True):
    from .qualification_v03 import verify_report
    if path is None:
        raise ValueError("numerical qualification requires passing fresh-pilot evidence first")
    report = read_sealed(path)
    if full:
        verify_report(report)
    else:
        _current_identity(report)
    if (report["checks"]["selection"]["selection_sha256"] != selection["report_sha256"]
            or report["inputs"]["coarse"] is not None or report["inputs"]["fine"] is not None
            or any(check["status"] != "passed" for name, check in report["checks"].items() if name != "numerical")):
        raise ValueError("fresh-pilot gates or anomaly dispositions do not admit numerical qualification")
    return report


def register_selected(selection, output, *, cohort="pilot", seeds=None, pilot_evidence=None):
    from .runner import register
    selection_path = Path(selection).resolve()
    result = verify_selection(selection_path)
    protocol = selected_protocol(result, refined=cohort == "numerical_fine")
    seeds = list(seeds) if seeds is not None else (list(range(5000, 5020)) if cohort == "pilot" else list(range(2000, 2050)))
    _validate_selected_seeds(cohort, seeds)
    pilot_report = _verify_pilot_admission(pilot_evidence, result) if cohort.startswith("numerical_") else None
    manifest = register(output, protocol, cohort, seeds=seeds)
    admission = seal({"study_version": "0.3", "selection_path": str(selection_path),
                      "selection_sha256": result["report_sha256"], "manifest_sha256": manifest["manifest_sha256"],
                      "pilot_evidence_path": str(Path(pilot_evidence).resolve()) if pilot_report else None,
                      "pilot_evidence_sha256": pilot_report["report_sha256"] if pilot_report else None})
    write_new(Path(output) / "study_admission.json", admission)
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    reg = sub.add_parser("register-grid"); reg.add_argument("output"); reg.add_argument("--controls", required=True)
    sel = sub.add_parser("select"); sel.add_argument("grid"); sel.add_argument("output")
    fresh = sub.add_parser("register-selected"); fresh.add_argument("selection"); fresh.add_argument("output")
    fresh.add_argument("--cohort", choices=("pilot", "numerical_coarse", "numerical_fine"), default="pilot")
    fresh.add_argument("--pairs", type=int, choices=(50, 200, 1000), default=50)
    fresh.add_argument("--pilot-evidence")
    args = parser.parse_args(argv)
    if args.command == "register-grid":
        result = register_grid(args.output, args.controls)
    elif args.command == "select":
        result = select_candidate(args.grid, args.output)
    else:
        result = register_selected(args.selection, args.output, cohort=args.cohort,
                                   seeds=None if args.cohort == "pilot" else list(range(2000, 2000 + args.pairs)),
                                   pilot_evidence=args.pilot_evidence)
    print(json.dumps({k: v for k, v in result.items() if k in ("status", "selected_candidate", "report_sha256", "manifest_sha256", "registration_sha256")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
