"""Prespecified paired finite-time persistence analysis.

The unit of resampling is a complete registered pair.  Numerical failures are
never outcomes.  These intervals describe Monte Carlo uncertainty conditional
on the declared synthetic model, not biological parameter uncertainty.
See DOCUMENTATION.md, Spatial nutrient integrity experiment.
"""

from __future__ import annotations

from collections import Counter
import math
from numbers import Integral, Real
from statistics import NormalDist
from typing import Iterable, Mapping

import numpy as np


class IncompleteExperimentError(ValueError):
    """The registered denominator cannot be reconstructed without omissions."""


def _integer(value, name, minimum=0):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _number(value, name, minimum=0, strict=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite number")
    value = float(value)
    if not math.isfinite(value) or value < minimum or (strict and value == minimum):
        raise ValueError(f"{name} is outside its finite domain")
    return value


def _confidence(value):
    value = _number(value, "confidence", 0, True)
    if value >= 1:
        raise ValueError("confidence must be between zero and one")
    return value


def validate_records(records: Iterable[Mapping], registered_seeds: Iterable[int], *,
                     cohort: str, horizon_s: float, dt_s: float):
    """Return seed-sorted (U, H) records after strict denominator validation.

The caller passes one final attempt per arm, selected by the runner's retained
completion markers. Passing attempt history directly is deliberately an error.
Optional ``initial_focal_count=0`` is rejected as not inoculated.
"""
    horizon_s = _number(horizon_s, "horizon_s", 0, True)
    dt_s = _number(dt_s, "dt_s", 0, True)
    steps = horizon_s / dt_s
    if not math.isclose(steps, round(steps), rel_tol=0, abs_tol=1e-8):
        raise ValueError("horizon_s must be an integral number of ticks")
    if not isinstance(cohort, str) or not cohort:
        raise ValueError("cohort must be explicit")
    seeds = [_integer(s, "seed") for s in registered_seeds]
    if not seeds or len(seeds) != len(set(seeds)):
        raise ValueError("registered seeds must be nonempty and unique")
    expected = {(s, a) for s in seeds for a in ("U", "H")}
    indexed = {}
    for record in records:
        if not isinstance(record, Mapping):
            raise IncompleteExperimentError("arm record must be a mapping")
        try:
            seed = _integer(record["seed"], "seed")
            arm = record["arm"]
            key = (seed, arm)
            if key not in expected:
                raise IncompleteExperimentError(f"unregistered arm {key!r}")
            if key in indexed:
                raise IncompleteExperimentError(f"duplicate arm {key!r}")
            if record["cohort"] != cohort:
                raise IncompleteExperimentError(f"mixed cohorts at {key!r}")
            if record["status"] != "completed" or record["validated"] is not True:
                raise IncompleteExperimentError(f"arm {key!r} is not validated and completed")
            reached = _number(record["reached_time_s"], "reached_time_s")
            if not math.isclose(reached, horizon_s, rel_tol=0, abs_tol=1e-9):
                raise IncompleteExperimentError(f"arm {key!r} did not reach the registered horizon")
            count = _integer(record["focal_count"], "focal_count")
            if "initial_focal_count" in record:
                _integer(record["initial_focal_count"], "initial_focal_count", 1)
            extinct = record["extinction_time_s"]
            if count > 0:
                if extinct is not None:
                    raise IncompleteExperimentError(f"survivor {key!r} has an extinction event")
            else:
                extinct = _number(extinct, "extinction_time_s", 0, True)
                if extinct > horizon_s + 1e-9:
                    raise IncompleteExperimentError(f"extinction after horizon at {key!r}")
                tick = extinct / dt_s
                if not math.isclose(tick, round(tick), rel_tol=0, abs_tol=1e-7):
                    raise IncompleteExperimentError(f"extinction not on a completed tick at {key!r}")
            indexed[key] = record
        except (KeyError, TypeError, ValueError) as exc:
            if isinstance(exc, IncompleteExperimentError):
                raise
            raise IncompleteExperimentError(f"invalid arm record: {exc}") from exc
    missing = sorted(expected - indexed.keys())
    if missing:
        raise IncompleteExperimentError(f"missing registered arms: {missing!r}")
    for field in ("protocol_sha256", "code_sha256", "manifest_sha256"):
        if any(field in r for r in indexed.values()):
            if not all(field in r for r in indexed.values()) or len({r[field] for r in indexed.values()}) != 1:
                raise IncompleteExperimentError(f"mixed or missing {field} provenance")
    if any("tail_min_focal_count" in r for r in indexed.values()):
        for record in indexed.values():
            try:
                minimum = _integer(record["tail_min_focal_count"], "tail_min_focal_count")
            except (KeyError, ValueError) as exc:
                raise IncompleteExperimentError("partial or invalid establishment window record") from exc
            if minimum > record["focal_count"]:
                raise IncompleteExperimentError("tail minimum cannot exceed the endpoint count")
    return [(indexed[(s, "U")], indexed[(s, "H")]) for s in sorted(seeds)]


def wilson_interval(successes: int, n: int, confidence: float = .95):
    n = _integer(n, "n", 1)
    successes = _integer(successes, "successes")
    if successes > n:
        raise ValueError("successes exceeds n")
    confidence = _confidence(confidence)
    z = NormalDist().inv_cdf((1 + confidence) / 2)
    p = successes / n
    denominator = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denominator
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return [max(0., center - half), min(1., center + half)]


def _percentile_bootstrap(values, *, resamples, seed, confidence):
    """Resample complete-pair contrasts in bounded-memory batches."""
    resamples = _integer(resamples, "bootstrap_resamples", 1)
    seed = _integer(seed, "bootstrap_seed")
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    means = np.empty(resamples)
    # The categorical multinomial is exactly the empirical-pair bootstrap for
    # a mean and avoids allocating resamples by n when n is large.
    unique, counts = np.unique(values, return_counts=True)
    probabilities = counts / len(values)
    for start in range(0, resamples, 1024):
        stop = min(start + 1024, resamples)
        frequencies = rng.multinomial(len(values), probabilities, size=stop - start)
        means[start:stop] = frequencies @ unique / len(values)
    alpha = 1 - confidence
    bounds = np.quantile(means, [alpha / 2, 1 - alpha / 2], method="linear")
    return [float(x) for x in bounds]


def analyze_pairs(records, registered_seeds, *, cohort, horizon_s, dt_s,
                  confidence=.95, bootstrap_resamples=10000,
                  bootstrap_seed=20261020):
    """Analyze all registered pairs, with no drop/replacement/arm matching."""
    confidence = _confidence(confidence)
    pairs = validate_records(records, registered_seeds, cohort=cohort,
                             horizon_s=horizon_s, dt_s=dt_s)
    counts = Counter({"both_survive": 0, "H_only": 0, "U_only": 0, "neither": 0})
    contrasts = []
    absolute_counts = []
    extinction_records = []
    for uniform, heterogeneous in pairs:
        u, h = int(uniform["focal_count"] > 0), int(heterogeneous["focal_count"] > 0)
        counts[{(1, 1): "both_survive", (0, 1): "H_only",
                (1, 0): "U_only", (0, 0): "neither"}[(u, h)]] += 1
        contrasts.append(h - u)
        absolute_counts.append({"seed": uniform["seed"], "U": uniform["focal_count"],
                                "H": heterogeneous["focal_count"]})
        for record in (uniform, heterogeneous):
            event = record["focal_count"] == 0
            extinction_records.append({"seed": record["seed"], "arm": record["arm"],
                                       "time_s": record["extinction_time_s"] if event else horizon_s,
                                       "event_observed": event, "right_censored": not event})
    n = len(pairs)
    estimate = sum(contrasts) / n
    half = math.sqrt(2 * math.log(2 / (1 - confidence)) / n)
    conservative = [max(-1., estimate - half), min(1., estimate + half)]
    percentile = _percentile_bootstrap(contrasts, resamples=bootstrap_resamples,
                                      seed=bootstrap_seed, confidence=confidence)
    degenerate = len(set(contrasts)) == 1
    arms = {}
    for arm in ("U", "H"):
        successes = counts["both_survive"] + counts[f"{arm}_only"]
        arms[arm] = {"survivors": successes, "n": n, "probability": successes / n,
                     "wilson_interval": wilson_interval(successes, n, confidence)}
    result = {
        "cohort": cohort, "n_pairs": n, "registered_seeds": [u["seed"] for u, _ in pairs],
        "horizon_s": horizon_s, "dt_s": dt_s, "confidence": confidence,
        "estimand": "P(focal survives at T | H) - P(focal survives at T | U)",
        "delta_p": estimate, "contingency": dict(counts), "arms": arms,
        "bootstrap": {"method": "paired_percentile", "unit": "whole_pair",
                      "resamples": bootstrap_resamples, "seed": bootstrap_seed,
                      "interval": percentile, "degenerate": degenerate},
        "hoeffding_interval": conservative,
        "reported_interval": conservative if degenerate else percentile,
        "reported_interval_method": "hoeffding_degenerate_fallback" if degenerate else "paired_percentile",
        "absolute_focal_counts": absolute_counts, "extinction_records": extinction_records,
        "interpretation": "Monte Carlo uncertainty conditional on the model and parameters; no biological validation",
    }
    if "tail_min_focal_count" in pairs[0][0]:
        result["secondary_establishment"] = {
            "window_fraction": [.8, 1.], "requires_all_completed_ticks": True,
            "thresholds": {
                str(threshold): {arm: sum(record["tail_min_focal_count"] >= threshold
                                         for pair in pairs for record in pair if record["arm"] == arm)
                                 for arm in ("U", "H")}
                for threshold in (1, 5, 10)
            },
        }
    return result
