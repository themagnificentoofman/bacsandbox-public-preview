"""Adversarial admission/qualification fixtures; never scientific evidence."""
import copy
import json
from pathlib import Path

import pytest

from . import qualification_run as legacy
from . import qualification_v03 as qv
from . import runner
from . import study_v03 as sv
from .protocol import digest, load_protocol


def audit_fixture(start=4000, *, survivors=8, birth_runs=10):
    rows = []
    for arm in ("U", "H"):
        for i, seed in enumerate(range(start, start + 20)):
            rows.append({"seed": seed, "arm": arm, "focal_count": int(i < survivors),
                         "births": {"focal": int(i < birth_runs), "resident": int(i < birth_runs)},
                         "capacity_culling_count": 0})
    return {"status": "passed", "seeds": list(range(start, start + 20)), "completed_pairs": 20,
            "by_arm": rows, "retained_attempts": []}


def test_v03_plan_rejects_rehashed_weakening(tmp_path):
    assert digest(sv.load_plan()) == sv.PLAN_SHA256
    plan = json.loads(sv.PLAN_PATH.read_text())
    plan["selection"]["eligibility"]["min_U_runs_with_both_focal_and_resident_births"] = 0
    path = tmp_path / "weakened.json"
    path.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="prespecified canonical plan"):
        sv.load_plan(path)


def test_v03_grid_registration_precedes_execution(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(sv, "_controls", lambda *_: {"status": "passed", "report_sha256": "c" * 64})
    def stub_register(path, protocol, cohort, seeds):
        path.mkdir()
        manifest = {"protocol_sha256": digest(protocol), "cohort": cohort, "seeds": seeds,
                    "code": runner.code_identity(), "runtime": runner.runtime_identity()}
        manifest["manifest_sha256"] = digest(manifest)
        runner.atomic_json(path / "manifest.json", manifest)
        with pytest.raises(FileNotFoundError):
            sv.verify_admission(path, manifest)
        calls.append(cohort)
        return manifest
    monkeypatch.setattr(runner, "register", stub_register)
    result = sv.register_grid(tmp_path / "grid", tmp_path / "controls")
    assert len(calls) == 8
    assert result["seeds"] == list(range(4000, 4020))
    for entry in result["candidates"]:
        sv.verify_admission(entry["path"], runner.read_manifest(entry["path"]))
    (tmp_path / "grid" / "C08" / "manifest.json").unlink()
    with pytest.raises(FileNotFoundError):
        sv.verify_admission(result["candidates"][0]["path"], runner.read_manifest(result["candidates"][0]["path"]))


def test_v03_eligibility_counts_only_u_both_populations():
    audit = audit_fixture(birth_runs=9)
    for row in audit["by_arm"]:
        if row["arm"] == "H":
            row["births"] = {"focal": 100, "resident": 100}
    decision = sv.candidate_decision(audit)
    assert not decision["eligible"]
    assert decision["U_denominator"] == 20
    assert decision["U_runs_with_both_population_births"] == 9
    audit = audit_fixture(birth_runs=20)
    for i, row in enumerate(audit["by_arm"][:20]):
        row["births"]["focal" if i < 10 else "resident"] = 0
    assert sv.candidate_decision(audit)["U_runs_with_both_population_births"] == 0
    audit = audit_fixture()
    audit["by_arm"][-1]["capacity_culling_count"] = 1
    assert not sv.candidate_decision(audit)["eligible"]


def selection_fixture(tmp_path, monkeypatch, eligible=()):
    plan = sv.load_plan()
    entries = [{"id": c["id"], "rank": c["rank"], "path": str(tmp_path / c["id"]),
                "manifest_sha256": str(c["rank"]) * 64} for c in plan["selection"]["candidates"]]
    index = {"plan_path": str(sv.PLAN_PATH), "candidates": entries, "registration_sha256": "r" * 64,
             "code": runner.code_identity(), "runtime": runner.runtime_identity()}
    monkeypatch.setattr(sv, "verify_grid", lambda *_, **__: index)
    monkeypatch.setattr(sv, "_source", lambda *_, **__: None)
    called = []
    def audited(path, sources, *, prefix, uniform_replay):
        called.append(prefix)
        audit = audit_fixture(survivors=8 if prefix in eligible else 0)
        return {"manifest_sha256": next(c["manifest_sha256"] for c in entries if c["id"] == prefix)}, [], audit
    monkeypatch.setattr(sv, "audit_cohort", audited)
    return called, index


def test_v03_no_eligible_candidate_stays_blocked(tmp_path, monkeypatch):
    called, _ = selection_fixture(tmp_path, monkeypatch)
    result = sv.select_candidate(tmp_path, tmp_path / "selection.json")
    assert called == [f"C{i:02d}" for i in range(1, 9)]
    assert result["status"] == "blocked" and result["selected_candidate"] is None
    with pytest.raises(ValueError, match="no eligible candidate"):
        sv.selected_protocol(result)


def test_v03_selection_uses_first_eligible_and_reaudits_rehashed_outcomes(tmp_path, monkeypatch):
    called, _ = selection_fixture(tmp_path, monkeypatch, eligible=("C02", "C04"))
    result = sv.select_candidate(tmp_path, tmp_path / "selection.json")
    assert len(called) == 8 and result["selected_candidate"] == "C02"
    tampered = copy.deepcopy(result)
    candidate = tampered["candidates"][0]
    candidate["audit"] = audit_fixture(survivors=8)
    candidate.update(sv.candidate_decision(candidate["audit"]))
    tampered["selected_candidate"] = "C01"
    path = tmp_path / "tampered.json"
    runner.atomic_json(path, sv.seal(tampered))
    with pytest.raises(ValueError, match="does not reproduce"):
        sv.verify_selection(path)


def test_v03_fresh_pilot_does_not_require_h_births():
    audit = audit_fixture(5000)
    for row in audit["by_arm"]:
        if row["arm"] == "H":
            row["focal_count"] = 0
            row["births"] = {"focal": 0, "resident": 0}
        else:
            row["focal_count"] = 1  # No holdout survival-range selection.
    assert qv.fresh_pilot_decision(audit)["status"] == "passed"
    audit["by_arm"][9]["births"]["resident"] = 0
    assert qv.fresh_pilot_decision(audit)["status"] == "failed"


def test_v03_freeze_rejects_unresolved_anomaly(tmp_path, monkeypatch):
    report = {"study_version": "0.3", "study_plan_sha256": sv.PLAN_SHA256,
              "checks": {"anomalies": legacy._bind({"status": "unresolved", "sources": []})}}
    report["verdict"] = qv.qualification_gate(report)
    report = sv.seal(report)
    assert "anomalies" in report["verdict"]["pending"]
    path = tmp_path / "qualification.json"
    runner.atomic_json(path, report)
    monkeypatch.setattr(qv, "verify_report", lambda value: value)
    with pytest.raises(ValueError, match="freeze blocked"):
        qv.freeze(path, tmp_path / "freeze.json")
    assert not (tmp_path / "freeze.json").exists()


def test_v03_numerical_requires_passing_pilot_and_prespecified_look(tmp_path):
    with pytest.raises(ValueError, match="fresh-pilot evidence"):
        sv._verify_pilot_admission(None, {})
    for seeds in (list(range(2000, 2051)), list(range(2001, 2051)), list(range(1000, 1050))):
        with pytest.raises(ValueError, match="declared cumulative"):
            sv._validate_selected_seeds("numerical_coarse", seeds)
    with pytest.raises(ValueError, match="fresh v0.3 pilot"):
        sv._validate_selected_seeds("pilot", list(range(20)))


def test_freeze_embedded_dict_tampering_rejected():
    for verifier in (legacy.verify_freeze, qv.verify_freeze):
        value = {"status": "qualified", "freeze_sha256": "0" * 64}
        original = copy.deepcopy(value)
        with pytest.raises(ValueError, match="freeze hash"):
            verifier(value, protocol_sha256="p" * 64, code_sha256="c" * 64, runtime={})
        assert value == original


def test_v03_main_rejects_legacy_freeze_and_rechecks_embedded_evidence(tmp_path, monkeypatch):
    protocol = load_protocol()
    protocol["amendments"].append({"from_version": "0.2", "change": "study_v03:C01", "reason": "synthetic unit fixture"})
    with pytest.raises(ValueError, match="v0.3 qualification freeze"):
        runner._verify_main_freeze({"status": "qualified"}, protocol, {"sha256": "c" * 64}, {})
    embedded = {"status": "qualified", "study_version": "0.3", "freeze_sha256": "0" * 64}
    manifest = {"cohort": "main", "freeze": embedded, "protocol": protocol,
                "code": runner.code_identity(), "runtime": runner.runtime_identity()}
    monkeypatch.setattr(runner, "read_manifest", lambda *_: manifest)
    with pytest.raises(ValueError, match="freeze hash"):
        runner.execute_arm(tmp_path, 1000, "U")
    assert not (tmp_path / "runs").exists()


def test_focal_exposure_zero_denominator_requires_null():
    protocol = load_protocol()
    focal = next(s for s in protocol["species"] if s["role"] == "focal")
    population = {"cell_seconds": {"rich": 0., "poor": 3., "total": 3.},
                  "uptake_pg": {"rich": 0., "poor": 6., "total": 6.},
                  "uptake_pg_per_cell_second": {"rich": None, "poor": 2., "total": 2.},
                  "rich_cell_time_fraction": 0., "living_count": 0, "initial_living_count": focal["init"],
                  "time_support": {"observed_through_seconds": 600., "first_extinction_time_seconds": 4.,
                                   "right_censored_at_observation_end": False,
                                   "first_occupied_tick_start_seconds": 0., "last_occupied_tick_end_seconds": 4.},
                  "residence_convention": "tick_start_position_and_abundance", "uptake_convention": "event_location_before_movement",
                  "cell_time_weighted_field_contrast": 0.}
    record = {"seed": 5000, "arm": "H", "focal_count": 0, "extinction_time_s": 4.,
              "diagnostics": {"early_exposure_contrast": 0.},
              "exposure": {"populations": {focal["name"]: population}, "early_window_seconds": [0, 25], "early_observed_seconds": 25.}}
    assert qv._focal_exposure_check([record], [], protocol)["status"] == "passed"
    population["uptake_pg_per_cell_second"]["rich"] = 0.
    assert qv._focal_exposure_check([record], [], protocol)["status"] == "failed"


def test_v03_numerical_gate_rejects_rehashed_rule_changes_and_reversed_intervals():
    check = {"status": "passed", "n_pairs": 1000, "registered_seeds": list(range(2000, 3000)),
             "registered_look": True, "Q_status": "passed", "Q_interval": [1., -1.],
             "tolerance": 2., "confidence_family": .50, "interval_method": "unregistered",
             "arms": {arm: {"status": "passed", "survival_change_interval": [0., 0.]} for arm in ("U", "H")}}
    report = {"study_version": "0.3", "study_plan_sha256": sv.PLAN_SHA256,
              "checks": {"numerical": legacy._bind(check)}}
    verdict = qv.qualification_gate(report)
    assert verdict["status"] == "failed"
    assert "numerical prespecified interval rule mismatch" in verdict["reasons"]
    assert "numerical interval outside tolerance" in verdict["reasons"]
