"""Qualification gates that retain uncertainty and bind evidence to code/config.

All thresholds and the staged sample sizes below are prespecified Own-model
qualification choices, not biological calibrations.  See DOCUMENTATION.md,
Spatial nutrient integrity experiment.
"""

from __future__ import annotations

import hashlib
import json
import math
import re

import numpy as np

from .analysis import (_confidence, _number, _percentile_bootstrap,
                       validate_records)


QUALIFICATION_SIZES = (50, 200, 1000)
QUALIFICATION_SEED_START = 2000
REQUIRED_CHECKS = (
    "lifecycle", "starvation", "matched_inputs", "food_ledger", "population_ledger",
    "scheduling", "replay", "null_allocator", "food_node_limit", "failure_retry",
    "endpoint", "exposure", "pilot", "numerical", "analysis", "disabled_mechanisms",
)
DIAGNOSTIC_KEYS = ("births", "deaths", "food_balance_residual_pg", "exposure_contrast")


def canonical_sha256(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode("utf-8")).hexdigest()


def bounded_discrete_mean_interval(values, support, confidence=.95):
    """Conservative finite-sample interval on a known finite support.

    Simultaneous Bonferroni Clopper--Pearson intervals enclose all category
    probabilities. Linear extrema over that box intersected with the simplex
    then enclose the expectation. Unseen categories remain possible; a zero
    sample variance never becomes proof of zero uncertainty.
    """
    from scipy.stats import beta

    confidence = _confidence(confidence)
    values = np.asarray(values, dtype=float)
    support = np.asarray(support, dtype=float)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError("values must be a nonempty finite vector")
    if support.ndim != 1 or len(support) < 2 or not np.isfinite(support).all():
        raise ValueError("support must be a finite ordered vector")
    if not np.all(np.diff(support) > 0) or not np.isin(values, support).all():
        raise ValueError("values must lie on the unique increasing support")
    n = len(values)
    counts = np.array([(values == x).sum() for x in support], dtype=int)
    tail = (1 - confidence) / (2 * len(support))
    lower = np.array([0. if k == 0 else beta.ppf(tail, k, n - k + 1) for k in counts])
    upper = np.array([1. if k == n else beta.ppf(1 - tail, k + 1, n - k) for k in counts])

    def extreme(order):
        probabilities = lower.copy()
        remaining = 1 - probabilities.sum()
        if remaining < -1e-10:
            raise ArithmeticError("confidence probability bounds have no feasible simplex")
        for index in order:
            extra = min(max(0., remaining), upper[index] - probabilities[index])
            probabilities[index] += extra
            remaining -= extra
        if remaining > 1e-10:
            raise ArithmeticError("confidence probability bounds have no feasible simplex")
        return float(probabilities @ support)

    return [extreme(range(len(support))), extreme(range(len(support) - 1, -1, -1))]


def _equivalence_status(interval, tolerance):
    if interval[0] >= -tolerance and interval[1] <= tolerance:
        return "passed"
    if interval[0] > tolerance or interval[1] < -tolerance:
        return "failed"
    return "inconclusive"


def compare_timesteps(coarse_records, fine_records, registered_seeds, *,
                      horizon_s=600., coarse_dt_s=.1, fine_dt_s=.05,
                      cohort="numerical", fine_cohort=None, tolerance=.05,
                      confidence=.95, bootstrap_resamples=10000,
                      bootstrap_seed=20261020, diagnostic_keys=DIAGNOSTIC_KEYS):
    """Compare paired H-U effects and arm survival on the same realizations.

    Three cumulative looks (50, 200, 1000) split alpha equally. Thus repeated
    attempts to qualify do not turn an uncertain first screen into an
    uncontrolled optional-stopping procedure. Arbitrary sample sizes can be
    summarized, but cannot produce a passing qualification gate.
    """
    confidence = _confidence(confidence)
    tolerance = _number(tolerance, "tolerance", 0, True)
    if tolerance > 2:
        raise ValueError("survival discrepancy tolerance exceeds its support")
    if not fine_dt_s < coarse_dt_s:
        raise ValueError("fine timestep must be smaller than coarse timestep")
    seeds = list(registered_seeds)
    coarse = validate_records(coarse_records, seeds, cohort=cohort,
                              horizon_s=horizon_s, dt_s=coarse_dt_s)
    fine = validate_records(fine_records, seeds, cohort=fine_cohort or cohort,
                            horizon_s=horizon_s, dt_s=fine_dt_s)
    n = len(coarse)
    c = np.array([[u["focal_count"] > 0, h["focal_count"] > 0] for u, h in coarse], dtype=int)
    f = np.array([[u["focal_count"] > 0, h["focal_count"] > 0] for u, h in fine], dtype=int)
    changes = f - c
    q = changes[:, 1] - changes[:, 0]
    look_confidence = 1 - (1 - confidence) / len(QUALIFICATION_SIZES)
    q_interval = bounded_discrete_mean_interval(q, [-2, -1, 0, 1, 2], look_confidence)
    q_status = _equivalence_status(q_interval, tolerance)
    arms = {}
    statuses = [q_status]
    for column, arm in enumerate(("U", "H")):
        interval = bounded_discrete_mean_interval(changes[:, column], [-1, 0, 1], look_confidence)
        status = _equivalence_status(interval, tolerance)
        statuses.append(status)
        diagnostics = {}
        for key in diagnostic_keys:
            before, after = [], []
            for cp, fp in zip(coarse, fine):
                for record, target in ((cp[column], before), (fp[column], after)):
                    try:
                        value = float(record["diagnostics"][key])
                    except (KeyError, ValueError, TypeError) as exc:
                        raise ValueError(f"missing numeric diagnostic {arm}/{key}") from exc
                    if not math.isfinite(value):
                        raise ValueError(f"nonfinite diagnostic {arm}/{key}")
                    target.append(value)
            difference = np.subtract(after, before)
            diagnostics[key] = {
                "coarse_mean": float(np.mean(before)), "fine_mean": float(np.mean(after)),
                "mean_change": float(np.mean(difference)),
                "descriptive_paired_percentile_interval": _percentile_bootstrap(
                    difference, resamples=bootstrap_resamples, seed=bootstrap_seed,
                    confidence=confidence),
                "qualification_threshold": None,
            }
        arms[arm] = {"coarse_survival": float(np.mean(c[:, column])),
                     "fine_survival": float(np.mean(f[:, column])),
                     "survival_change": float(np.mean(changes[:, column])),
                     "survival_change_interval": interval, "status": status,
                     "diagnostics": diagnostics}
    registered_look = n in QUALIFICATION_SIZES and sorted(seeds) == list(
        range(QUALIFICATION_SEED_START, QUALIFICATION_SEED_START + n))
    status = ("failed" if "failed" in statuses else
              "passed" if all(s == "passed" for s in statuses) and registered_look else
              "inconclusive")
    next_size = next((size for size in QUALIFICATION_SIZES if size > n), None)
    return {
        "status": status, "n_pairs": n, "registered_seeds": sorted(seeds),
        "registered_look": registered_look, "coarse_dt_s": coarse_dt_s,
        "fine_dt_s": fine_dt_s, "horizon_s": horizon_s, "tolerance": tolerance,
        "confidence_family": confidence, "confidence_per_look": look_confidence,
        "interval_method": "Bonferroni_Clopper_Pearson_category_simplex",
        "Q_mean": float(np.mean(q)), "Q_interval": q_interval,
        "Q_status": q_status, "Q_counts": {str(k): int(np.sum(q == k)) for k in range(-2, 3)},
        "arms": arms,
        "next_cumulative_pair_count": next_size if status == "inconclusive" else None,
        "next_registered_seeds": list(range(QUALIFICATION_SEED_START,
                                            QUALIFICATION_SEED_START + next_size))
                                 if status == "inconclusive" and next_size else [],
        "decision": ("numerical survival criteria met; inspect ecological diagnostics"
                     if status == "passed" else
                     "numerical discrepancy exceeds tolerance" if status == "failed" else
                     "uncertainty too large or sample is not a registered look; main freeze blocked"),
    }


def qualification_gate(evidence, *, protocol_sha256, code_sha256):
    """Require complete, version-bound evidence before a main-cohort freeze.

    Checks are evidence records, not free-form assurances. Their referenced
    hashes must be resolved and verified by the runner before this function is
    called. This function evaluates the record schema and decision rules.
    """
    failures, pending = [], []

    def digest(value):
        return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None

    for key, expected in (("protocol_sha256", protocol_sha256), ("code_sha256", code_sha256)):
        if not digest(expected) or evidence.get(key) != expected:
            failures.append(f"{key} does not match the candidate freeze")
    checks = evidence.get("checks", {})
    if not isinstance(checks, dict):
        return {"status": "failed", "freeze_allowed": False,
                "reasons": failures + ["checks must be a mapping"], "pending": []}
    for name in REQUIRED_CHECKS:
        record = checks.get(name)
        if not isinstance(record, dict):
            pending.append(name)
            continue
        status = record.get("status")
        if status == "failed":
            failures.append(f"{name} failed")
        elif status != "passed":
            pending.append(name)
        if not digest(record.get("evidence_sha256")):
            failures.append(f"{name} lacks a retained evidence digest")
    pilot = checks.get("pilot", {})
    if pilot.get("status") == "passed":
        if pilot.get("seeds") != list(range(20)) or pilot.get("completed_pairs") != 20:
            failures.append("pilot must include all 20 separate registered pairs")
    numerical = checks.get("numerical", {})
    if numerical.get("status") == "passed":
        n = numerical.get("n_pairs")
        if type(n) is not int or n not in QUALIFICATION_SIZES or numerical.get("registered_seeds") != list(
                range(QUALIFICATION_SEED_START, QUALIFICATION_SEED_START + n)):
            failures.append("numerical evidence is not a prespecified cumulative look")
        if numerical.get("registered_look") is not True or numerical.get("Q_status") != "passed":
            failures.append("numerical uncertainty criterion was not met")
        for arm in ("U", "H"):
            if numerical.get("arms", {}).get(arm, {}).get("status") != "passed":
                failures.append(f"{arm} numerical survival criterion was not met")
    exposure = checks.get("exposure", {})
    if exposure.get("status") == "passed":
        criterion = exposure.get("criterion", {})
        if not isinstance(criterion, dict):
            criterion = {}
        metric = criterion.get("metric")
        # This gate concerns resource geography. It must never select settings
        # because they produced the desired persistence effect.
        if not isinstance(metric, str) or not metric or any(
                word in metric.lower() for word in ("surviv", "persistence", "delta_p")):
            failures.append("exposure criterion must concern measured resource geography")
        try:
            threshold = _number(criterion["threshold"], "exposure threshold")
            measured = _number(exposure["measured"], "exposure measurement")
            if criterion.get("operator") != ">=" or measured < threshold:
                failures.append("exposure criterion was not met")
        except (KeyError, TypeError, ValueError):
            failures.append("exposure needs a prespecified threshold and finite measurement")
        if not digest(exposure.get("criterion_sha256")) or exposure.get("criterion_sha256") != canonical_sha256(criterion):
            failures.append("exposure criterion digest missing or mismatched")
    status = "failed" if failures else "inconclusive" if pending else "passed"
    return {"status": status, "freeze_allowed": status == "passed",
            "reasons": failures, "pending": pending,
            "protocol_sha256": protocol_sha256, "code_sha256": code_sha256}
