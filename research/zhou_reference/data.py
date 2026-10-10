"""Pinned-source inputs and role-restricted data for the Zhou reference.

Source zeros stay literal in the source mapping. Initialization alternatives are
explicit Own-model counterfactuals; see DOCUMENTATION.md section 10.7.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from research.data_readiness import map_benchmark as mapping

PHASES = mapping.PHASES
META = mapping.META
DEFAULT_SOURCE = mapping.HERE / "raw/supplements"


def _mapping(source_dir, species_path=None):
    return mapping.build(Path(source_dir), species_path or mapping.ROOT / "species_data.csv")


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def load_baselines(source_dir=DEFAULT_SOURCE, species_path=None):
    """Return only baseline data: follow-up observations cannot reach fitting here."""
    mapped = _mapping(source_dir, species_path)
    roster = mapped["strain_crosswalk.csv"]
    strains = tuple(row["strain_id"] for row in roster)
    target_index, = [i for i, row in enumerate(roster) if row["strain"] == "BgF5_2"]
    by_animal = {}
    for row in mapped["baseline_observations.csv"]:
        by_animal.setdefault(row["animal_id"], {})[row["strain_id"]] = row
    animals = []
    for animal_id, rows in sorted(by_animal.items()):
        ordered = [rows[strain] for strain in strains]
        if any(r["measurement_state"] not in ("finite_nonzero_reported_log",
                                              "log_zero_semantics_unresolved") for r in ordered):
            raise ValueError("unsupported baseline measurement state")
        target = ordered[target_index]
        if target["measurement_state"] != "finite_nonzero_reported_log":
            raise ValueError("target baseline must be finite, nonzero reported log")
        animals.append({"animal_id": animal_id, **{k: target[k] for k in META},
                        "baseline_log10": tuple(r["raw_log10"] for r in ordered),
                        "zero_mask": tuple(r["measurement_state"] == "log_zero_semantics_unresolved"
                                           for r in ordered),
                        "baseline_source_cells": tuple(r["source_cell"] for r in ordered),
                        "target_log10": target["raw_log10"]})
    if len(animals) != 39 or sum(sum(a["zero_mask"]) for a in animals) != 49:
        raise ValueError("baseline cohort contract changed")
    hashes = mapped["mapping_audit.json"]["source_hashes"]
    result = {"animals": tuple(animals), "strain_ids": strains,
              "target_index": target_index, "source_hashes": hashes,
              "source_fingerprint": fingerprint(hashes)}
    result["baseline_fingerprint"] = fingerprint(result)
    return result


def _is_calibration(row):
    return (row["experiment"] == "Iteration 1" and row["community"] == mapping.COMMUNITY
            and row["diet"] == mapping.DIETS[1])


def _load_role(source_dir, species_path, role):
    rows = []
    for row in _mapping(source_dir, species_path)["target_observations.csv"]:
        if row["phase_label"] == PHASES[0]:
            continue
        if row["measurement_state"] == "missing_source_cell":
            # The mapping enforces the one exact source-supported missing endpoint.
            continue
        if row["measurement_state"] != "finite_nonzero_reported_log":
            raise ValueError("unresolved target follow-up")
        allowed = _is_calibration(row)
        if allowed != (role == "calibration"):
            continue
        rows.append({"animal_id": row["animal_id"], **{k: row[k] for k in META},
                     "phase_label": row["phase_label"], "raw_log10": row["raw_log10"],
                     "source_cell": row["source_cell"], "source_sha256": row["source_sha256"],
                     "role": role})
    expected = 20 if role == "calibration" else 57
    if len(rows) != expected:
        raise ValueError(f"{role} count changed: {len(rows)} != {expected}")
    return tuple(rows)


def load_calibration(source_dir=DEFAULT_SOURCE, species_path=None):
    """Select control follow-ups by actual source metadata, never a caller role flag."""
    rows = _load_role(source_dir, species_path, "calibration")
    validate_calibration(rows, load_baselines(source_dir, species_path))
    return rows


def load_evaluation(source_dir=DEFAULT_SOURCE, species_path=None):
    """Explicit separate read for post-freeze reporting; never a fitting input."""
    return _load_role(source_dir, species_path, "evaluation")


def validate_calibration(rows, baselines):
    animals = {a["animal_id"]: a for a in baselines["animals"]}
    expected = {(a["animal_id"], phase) for a in animals.values() if _is_calibration(a)
                for phase in PHASES[1:]}
    actual = []
    for row in rows:
        animal = animals.get(row.get("animal_id"))
        if (animal is None or not _is_calibration(animal) or not _is_calibration(row)
                or row.get("role") != "calibration"
                or any(row[k] != animal[k] for k in META)):
            raise ValueError("calibration contains evaluation or mismatched identity")
        value = row["raw_log10"]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError("nonfinite calibration measurement")
        actual.append((row["animal_id"], row["phase_label"]))
    if len(actual) != 20 or len(set(actual)) != 20 or set(actual) != expected:
        raise ValueError("calibration must contain the exact 20 permitted endpoints")
    return True


def phase_index(label):
    if label not in PHASES[1:]:
        raise ValueError(f"not a follow-up phase: {label}")
    return PHASES[1:].index(label)


def get_times(experiment, chronology):
    if experiment not in ("Iteration 1", "Iteration 2"):
        raise ValueError("unknown iteration")
    if chronology in ("label", "phase_label", "phase_labels", "workbook", "workbook_label"):
        return (10.0, 28.0)
    if chronology in ("age", "age_derived"):
        return (11.0, 29.0) if experiment == "Iteration 1" else (12.0, 30.0)
    raise ValueError("unknown chronology")


def initialize(animal, scenario):
    """Apply a declared global scale without normalizing animal total biomass."""
    beta = scenario.get("biomass_scale", 1.0)
    scale = scenario.get("abundance_scale", 1e10)
    policy = scenario["zero_policy"]
    alternatives = {"absent": 0.0, "low": 1e-4, "high": 1e-2}
    if scenario.get("zero_relative", alternatives) != alternatives:
        raise ValueError("counterfactual initialization values differ from frozen contract")
    if policy not in alternatives:
        raise ValueError("unrecognized counterfactual zero policy")
    for name, value in (("biomass_scale", beta), ("abundance_scale", scale)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"invalid {name}")
    logs, zeros = animal["baseline_log10"], animal["zero_mask"]
    if len(logs) != 18 or len(zeros) != 18:
        raise ValueError("18 baseline slots required")
    result = []
    for raw, unresolved in zip(logs, zeros):
        if not isinstance(unresolved, bool) or isinstance(raw, bool) or not math.isfinite(raw):
            raise ValueError("invalid baseline value or zero mask")
        if unresolved != (raw == 0):
            raise ValueError("baseline zero semantics cannot be silently changed")
        value = alternatives[policy] if unresolved else 10.0 ** (raw - math.log10(scale))
        result.append(beta * value)
    array = np.asarray(result, dtype=float)
    if not np.all(np.isfinite(array)) or np.any(array < 0):
        raise ValueError("invalid initialized state")
    return array


def role_manifest(baselines, calibration, evaluation=None):
    validate_calibration(calibration, baselines)
    result = {"baseline_animals": len(baselines["animals"]), "baseline_slots": 702,
              "ambiguous_source_zeros": 49, "calibration_followups": len(calibration),
              "source_hashes": baselines["source_hashes"],
              "baseline_fingerprint": baselines["baseline_fingerprint"],
              "calibration_fingerprint": fingerprint(calibration)}
    if evaluation is not None:
        all_rows = tuple(calibration) + tuple(evaluation)
        keys = [(r["animal_id"], r["phase_label"]) for r in all_rows]
        if len(evaluation) != 57 or len(set(keys)) != 77:
            raise ValueError("evaluation split overlap/count mismatch")
        result["evaluation_followups"] = len(evaluation)
        result["endpoint_counts"] = dict(Counter(r["phase_label"] for r in all_rows))
        result["evaluation_fingerprint"] = fingerprint(evaluation)
    return result
