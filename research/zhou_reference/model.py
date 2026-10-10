"""Deterministic, scaled Zhou reference; see DOCUMENTATION.md section 10.7.

This research model is independent of the production pg-based simulator. Two
log growth factors represent the target and the exactly tied background group;
all 18 original labels are reconstructed without abundance rounding.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
import math
from collections.abc import Mapping, Sequence
from numbers import Real

import numpy as np
from scipy.integrate import solve_ivp


N_LABELS = 18
ACCOUNTING_TOLERANCE = 1.0e-7
LOG_CONVERGENCE_TOLERANCE = 1.0e-6
STAGE_DOMAIN_RETRY_MAX_STEPS = (1.0, 0.5, 0.25)
TIGHT_STAGE_DOMAIN_RETRY_MAX_STEPS = (0.5, 0.25, 0.125)


class ModelNumericalError(RuntimeError):
    """A retained, reportable numerical failure, never a clipped trajectory."""

    def __init__(self, message: str, integration_attempts=None):
        super().__init__(message)
        self.integration_attempts = list(integration_attempts or [])


class StageDomainError(ModelNumericalError):
    """A finite internal trial crossed a Monod pole; no accepted state is implied."""

    def __init__(self, time: float, state: np.ndarray, half: np.ndarray):
        super().__init__("Finite internal solver stage has a nonpositive Monod denominator")
        self.stage = {"time_day": float(time), "state": np.asarray(state).tolist(),
                      "resources": np.asarray(state[2:4]).tolist(),
                      "half_saturation": np.asarray(half).tolist(),
                      "denominators": (half + state[2:4]).tolist()}


@dataclass(frozen=True)
class ModelParameters:
    """Own-model parameters in declared scaled units and physical days.

    Resource and biomass equivalents are distinct units; yield converts them.
    These defaults are an exploratory center, not measured mouse kinetics.
    Supplementation adds only ``input_G`` after ``glucomannan_start_day``.
    ``initial_G`` applies to both arms, independently of the diet label.
    """

    qmax_R_target: float = 1.0
    qmax_R_background: float = 1.0
    qmax_G_target: float = 1.0
    qmax_G_background: float = 0.5
    K_R: float = 1.0
    K_G: float = 1.0
    Y_R_target: float = 1.0
    Y_R_background: float = 1.0
    Y_G_target: float = 1.0
    Y_G_background: float = 1.0
    loss: float = 0.5
    input_R: float = 1.0
    input_G: float = 0.1
    washout_R: float = 1.0
    washout_G: float = 1.0
    initial_R: float = 1.0
    initial_G: float = 0.0
    glucomannan_start_day: float = 0.0

    def __post_init__(self) -> None:
        for field in fields(self):
            value = getattr(self, field.name)
            if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
                    or not math.isfinite(value)):
                raise ValueError(f"{field.name} must be a finite real number")
            if value < 0 or (field.name.startswith("K_") and value == 0):
                raise ValueError(f"{field.name} has an invalid nonpositive value")


@dataclass(frozen=True)
class SimulationResult:
    times: np.ndarray
    biomass: np.ndarray
    resources: np.ndarray
    target_log_change: np.ndarray
    log_factors: np.ndarray
    accounting: dict[str, np.ndarray]
    diagnostics: dict[str, object]


def _parameters(value: ModelParameters | Mapping[str, float]) -> ModelParameters:
    if isinstance(value, ModelParameters):
        return value
    if isinstance(value, Mapping):
        return ModelParameters(**value)
    raise TypeError("parameters must be ModelParameters or a parameter mapping")


def _positive_values(log_initial: np.ndarray, log_factors: np.ndarray) -> np.ndarray:
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result = np.exp(log_initial + log_factors)
    if not np.all(np.isfinite(result)) or np.any(result <= 0):
        raise ModelNumericalError("Positive biomass overflowed or underflowed; no extinction was imposed")
    return result


def simulate(
    initial_biomass: Sequence[float],
    target_index: int,
    times: Sequence[float],
    parameters: ModelParameters | Mapping[str, float],
    supplemented: bool,
    tighter: bool = False,
) -> SimulationResult:
    """Integrate from day zero, returning predictions at increasing ``times``.

    Accepted solver states AND requested dense-output states must be finite and
    nonnegative. All retained positive biomass must remain representable. State
    failures are raised, with no clipping, implicit floor or alternate solver.
    Only finite internal Monod-pole trials may restart the same segment, with
    the bounded step-limit sequence declared in DOCUMENTATION.md section 10.7.
    A discontinuous supplement start is an exact integration boundary.
    """
    p = _parameters(parameters)
    initial = np.asarray(initial_biomass, dtype=float)
    requested = np.asarray(times, dtype=float)
    if initial.shape != (N_LABELS,) or not np.all(np.isfinite(initial)) or np.any(initial < 0):
        raise ValueError("initial_biomass must contain 18 finite nonnegative values")
    if (isinstance(target_index, (bool, np.bool_))
            or not isinstance(target_index, (int, np.integer))
            or not 0 <= target_index < N_LABELS):
        raise ValueError("target_index must be an integer from 0 through 17")
    if initial[target_index] <= 0:
        raise ValueError("target initialization must be positive for the log observation")
    if (requested.ndim != 1 or not len(requested) or not np.all(np.isfinite(requested))
            or np.any(requested < 0) or np.any(np.diff(requested) <= 0)):
        raise ValueError("times must be a nonempty strictly increasing finite nonnegative sequence")
    if not isinstance(supplemented, (bool, np.bool_)) or not isinstance(tighter, (bool, np.bool_)):
        raise ValueError("supplemented and tighter must be boolean")

    background = np.arange(N_LABELS) != target_index
    group_initial = np.array([initial[target_index], np.sum(initial[background])])
    if not np.all(np.isfinite(group_initial)):
        raise ValueError("initial group biomass totals must be finite")
    positive_groups = group_initial > 0
    group_log_initial = np.zeros(2)
    group_log_initial[positive_groups] = np.log(group_initial[positive_groups])
    label_groups = np.ones(N_LABELS, dtype=int)
    label_groups[target_index] = 0
    positive_labels = initial > 0
    label_log_initial = np.log(initial[positive_labels])

    # State: target/background log factors, R/G, then cumulative inputs,
    # resource outflows, resource uptakes, group growth and group losses.
    state = np.zeros(14)
    state[2:4] = (p.initial_R, p.initial_G)
    initial_resources = state[2:4].copy()
    qmax = np.array([[p.qmax_R_target, p.qmax_G_target],
                     [p.qmax_R_background, p.qmax_G_background]])
    yields = np.array([[p.Y_R_target, p.Y_G_target],
                      [p.Y_R_background, p.Y_G_background]])
    half = np.array([p.K_R, p.K_G])
    washout = np.array([p.washout_R, p.washout_G])
    rtol, atol, max_step = (1e-10, 1e-12, 0.5) if tighter else (1e-8, 1e-10, 1.0)
    outputs = np.zeros((len(requested), len(state)))
    accepted = [state.reshape(1, -1).copy()]
    evaluations = 0
    messages: list[str] = []
    integration_attempts: list[dict[str, object]] = []
    step_limits = TIGHT_STAGE_DOMAIN_RETRY_MAX_STEPS if tighter else STAGE_DOMAIN_RETRY_MAX_STEPS

    def group_biomass(factors: np.ndarray) -> np.ndarray:
        shape = factors.shape
        values = np.zeros(shape)
        values[..., positive_groups] = _positive_values(
            group_log_initial[positive_groups], factors[..., positive_groups])
        return values

    def rhs(time: float, current: np.ndarray, inputs: np.ndarray) -> np.ndarray:
        if not np.all(np.isfinite(current)):
            raise ModelNumericalError("Nonfinite solver stage")
        x = group_biomass(current[:2])
        resource = current[2:4]
        denominator = half + resource
        if np.any(denominator <= 0):
            raise StageDomainError(time, current, half)
        q = qmax * (resource / denominator)
        uptake = (q * x[:, None]).sum(axis=0)
        growth_rates = (yields * q).sum(axis=1)
        result = np.empty(14)
        result[:2] = growth_rates - p.loss
        result[2:4] = inputs - washout * resource - uptake
        result[4:6] = inputs
        result[6:8] = washout * resource
        result[8:10] = uptake
        result[10:12] = growth_rates * x
        result[12:14] = p.loss * x
        return result

    last_time = float(requested[-1])
    boundaries = [0.0]
    if 0 < p.glucomannan_start_day < last_time and supplemented and p.input_G > 0:
        boundaries.append(p.glucomannan_start_day)
    boundaries.append(last_time)
    if last_time == 0:
        outputs[0] = state
    for start, stop in zip(boundaries[:-1], boundaries[1:]):
        if stop == start:
            continue
        inputs = np.array([p.input_R, p.input_G if supplemented and start >= p.glucomannan_start_day else 0.0])
        for attempt_number, step_limit in enumerate(step_limits, start=1):
            attempt = {"segment_start_day": float(start), "segment_stop_day": float(stop),
                       "attempt": attempt_number, "max_step_days": step_limit,
                       "rhs_calls": 0, "success": False}

            def counted_rhs(t, y):
                attempt["rhs_calls"] += 1
                return rhs(t, y, inputs)

            try:
                # A retry starts from the identical segment boundary, never from
                # the failed trial. Equations and error tolerances are unchanged.
                solution = solve_ivp(
                    counted_rhs, (start, stop), state.copy(),
                    method="DOP853", rtol=rtol, atol=atol, max_step=step_limit,
                    dense_output=True,
                )
            except StageDomainError as exc:
                attempt.update(failure_type=type(exc).__name__, failure=str(exc), stage=exc.stage)
                integration_attempts.append(attempt)
                evaluations += attempt["rhs_calls"]
                if attempt_number == len(step_limits):
                    raise ModelNumericalError(
                        "Bounded internal-stage domain retries exhausted", integration_attempts) from exc
                continue
            except Exception as exc:
                attempt.update(failure_type=type(exc).__name__, failure=str(exc))
                integration_attempts.append(attempt)
                if isinstance(exc, ModelNumericalError):
                    exc.integration_attempts = list(integration_attempts)
                raise
            attempt.update(success=bool(solution.success), message=solution.message)
            integration_attempts.append(attempt)
            evaluations += attempt["rhs_calls"]
            break
        messages.append(solution.message)
        if not solution.success:
            raise ModelNumericalError(f"DOP853 failed: {solution.message}", integration_attempts)
        accepted.append(solution.y.T)
        mask = (requested >= start) & (requested <= stop)
        if np.any(mask):
            outputs[mask] = solution.sol(requested[mask]).T
        state = solution.y[:, -1]

    checked = np.concatenate([*accepted, outputs], axis=0)
    if not np.all(np.isfinite(checked)):
        raise ModelNumericalError("Nonfinite accepted or requested state", integration_attempts)
    if np.any(checked[:, 2:] < 0):
        raise ModelNumericalError("Negative accepted or requested resource/accounting state; not clipped", integration_attempts)
    try:
        checked_groups = group_biomass(checked[:, :2])
        # Check rare label representability at every accepted or requested state.
        _positive_values(label_log_initial, checked[:, label_groups[positive_labels]])
    except ModelNumericalError as exc:
        exc.integration_attempts = list(integration_attempts)
        raise
    resource_residual = checked[:, 2:4] - initial_resources - checked[:, 4:6] + checked[:, 6:8] + checked[:, 8:10]
    biomass_residual = checked_groups - group_initial - checked[:, 10:12] + checked[:, 12:14]
    resource_scale = np.maximum(1.0, initial_resources + checked[:, 4:6])
    biomass_scale = np.maximum(1.0, group_initial + checked[:, 10:12])
    max_resource_error = float(np.max(np.abs(resource_residual) / resource_scale))
    max_biomass_error = float(np.max(np.abs(biomass_residual) / biomass_scale))
    if max(max_resource_error, max_biomass_error) > ACCOUNTING_TOLERANCE:
        raise ModelNumericalError(
            f"Accounting residual exceeds {ACCOUNTING_TOLERANCE}: "
            f"resource={max_resource_error}, biomass={max_biomass_error}", integration_attempts)

    biomass = np.zeros((len(requested), N_LABELS))
    biomass[:, positive_labels] = _positive_values(label_log_initial, outputs[:, label_groups[positive_labels]])
    # Preserve the caller's exact initial representation at the exact baseline.
    biomass[requested == 0] = initial
    diagnostics: dict[str, object] = {
        "method": "DOP853", "rtol": rtol, "atol": atol, "max_step_days": max_step,
        "nfev": evaluations, "accepted_steps": sum(len(chunk) - 1 for chunk in accepted[1:]),
        "checked_states": len(checked), "minimum_resource": float(checked[:, 2:4].min()),
        "minimum_positive_biomass": float(biomass[:, positive_labels].min()),
        "max_scaled_resource_accounting_residual": max_resource_error,
        "max_scaled_biomass_accounting_residual": max_biomass_error,
        "accounting_tolerance": ACCOUNTING_TOLERANCE,
        "exact_zero_labels": int(np.sum(~positive_labels)), "messages": messages,
        "integration_attempts": integration_attempts,
        "stage_domain_retries": sum("stage" in attempt for attempt in integration_attempts),
    }
    return SimulationResult(
        times=requested.copy(), biomass=biomass, resources=outputs[:, 2:4].copy(),
        target_log_change=outputs[:, 0] / math.log(10), log_factors=outputs[:, :2].copy(),
        accounting={"input": outputs[:, 4:6].copy(), "outflow": outputs[:, 6:8].copy(),
                    "uptake": outputs[:, 8:10].copy(), "growth": outputs[:, 10:12].copy(),
                    "loss": outputs[:, 12:14].copy()},
        diagnostics=diagnostics,
    )


def check_convergence(standard: SimulationResult, tighter: SimulationResult) -> dict[str, object]:
    """Compare scored target log changes; contrast checks belong to the runner."""
    if not np.array_equal(standard.times, tighter.times):
        raise ValueError("Convergence comparison requires identical requested times")
    difference = float(np.max(np.abs(standard.target_log_change - tighter.target_log_change)))
    return {"max_target_log10_difference": difference,
            "log10_tolerance": LOG_CONVERGENCE_TOLERANCE,
            "passed": bool(difference <= LOG_CONVERGENCE_TOLERANCE)}
