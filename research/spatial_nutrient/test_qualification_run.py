"""Evidence chain tests; synthetic JUnit fixtures are never scientific evidence."""
import json
from pathlib import Path
from unittest.mock import patch
import xml.etree.ElementTree as ET

import pytest

from . import qualification_run as qr
from . import runner
from .protocol import load_protocol
from .qualification import canonical_sha256
from .test_analysis import records_for


def xml_fixture(path, *, missing=None, failed=None):
    suite = ET.Element("testsuite")
    names = sorted({name for required in qr.TEST_REQUIREMENTS.values() for name in required})
    for name in names:
        if name == missing:
            continue
        case = ET.SubElement(suite, "testcase", name=name, classname="synthetic_unit_fixture")
        if name == failed:
            ET.SubElement(case, "failure", message="injected")
    ET.ElementTree(suite).write(path)


def test_junit_requires_exact_observed_tests_and_rejects_failure(tmp_path):
    xml = tmp_path / "fixture.xml"
    xml_fixture(xml, missing="test_duplicate_uniform_inputs_replay_exactly")
    checks = qr.junit_checks(xml)
    assert checks["replay"]["status"] == "not_run"
    assert checks["lifecycle"]["status"] == "passed"
    xml_fixture(xml, failed="test_duplicate_uniform_inputs_replay_exactly")
    assert qr.junit_checks(xml)["replay"]["status"] == "failed"


def test_qualify_without_numerical_evidence_retains_blocked_decision(tmp_path):
    pilot = tmp_path / "pilot"
    pilot.mkdir()
    (pilot / "manifest.json").write_text("synthetic unit fixture only")
    protocol = load_protocol()
    code, runtime = runner.code_identity(), runner.runtime_identity()
    manifest = {"code": code, "runtime": runtime, "cohort": "pilot", "seeds": list(range(20)),
                "protocol": protocol, "manifest_sha256": "a" * 64}
    records = records_for([(1, 1)] * 20, horizon=600.)
    for record in records:
        record["diagnostics"]["exposure_contrast"] = .5 if record["arm"] == "H" else .1
    xml = tmp_path / "fixture.xml"
    xml_fixture(xml)
    metadata = {"code_sha256": code["sha256"], "runtime": runtime,
                "test_source_sha256": qr._test_identity(), "unchanged_during_tests": True,
                "returncode": 0, "report_sha256": qr.file_sha256(xml)}
    Path(str(xml) + ".meta.json").write_text(json.dumps(metadata))
    audited = {key: [{"passed": True}] for key in ("food_ledger", "population_ledger", "disabled_mechanisms", "pruning")}
    audited["births_by_arm"] = {"U": 1, "H": 1}
    out = tmp_path / "qualification.json"
    with patch.object(runner, "collect", return_value=(manifest, records)), patch.object(qr, "_audit_pilot", return_value=audited):
        result = qr.qualify(pilot, xml, out)
    assert out.exists()
    assert result["checks"]["numerical"]["status"] == "not_run"
    assert result["checks"]["exposure"]["status"] == "passed"
    assert result["verdict"]["status"] == "inconclusive"
    assert not result["verdict"]["freeze_allowed"]
    with pytest.raises(ValueError, match="freeze blocked"):
        qr.freeze(out, tmp_path / "freeze.json")


def test_retained_source_tampering_is_detected(tmp_path):
    source = tmp_path / "source.json"
    source.write_text("original")
    check = qr._bind({"status": "passed", "sources": ["sample"]})
    report = {"checks": {"sample": check}, "sources": {"sample": {"path": str(source), "sha256": qr.file_sha256(source)}}}
    report["report_sha256"] = canonical_sha256(report)
    assert qr.verify_report(report) == report
    source.write_text("changed")
    with pytest.raises(ValueError, match="source changed"):
        qr.verify_report(report)


def test_freeze_without_qualification_provenance_rejected(tmp_path):
    path = tmp_path / "forged.json"
    path.write_text(json.dumps({"status": "qualified", "protocol_sha256": "a" * 64,
                                "code_sha256": "b" * 64, "runtime": {}}))
    with pytest.raises(ValueError, match="freeze hash"):
        qr.verify_freeze(path, protocol_sha256="a" * 64, code_sha256="b" * 64, runtime={})
