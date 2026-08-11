"""Transactional integration between live spatial state and reduced FBA batches."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import logging
import math
from numbers import Integral, Real
from typing import Any, Callable, Dict, Mapping, Optional, Protocol, Sequence, Tuple

import numpy as np

from .diet_translator import DietTranslationEngine
from .ecm import ECM, SubstrateLayerSnapshot, SubstrateLayerSpec
from .metabolism import MetabolicOrchestrator, SUBSTRATE_ORDER


GridCoordinate = Tuple[int, int]
BIOMASS_CONSERVATION_ATOL_GDW = 1.0e-18


class FbaConservationError(RuntimeError):
    """Raised when a live batch violates represented-substrate conservation."""


class FbaLifecyclePolicy(Protocol):
    """Own-model policy for initializing and later partitioning live FBA biomass."""

    def initialize_agent(self, bacterium: Any) -> "FbaAgentState":
        ...

    def evaluate(
        self,
        bacterium: Any,
        state: "FbaAgentState",
        *,
        dt_seconds: float,
        environmental_fitness: float,
        carrying_capacity_factor: float,
    ) -> "FbaPhysiology":
        ...

    def partition_for_fission(
        self,
        parent: Any,
        child: Any,
        state: "FbaAgentState",
    ) -> Tuple["FbaAgentState", "FbaAgentState"]:
        ...


@dataclass(frozen=True)
class FbaSpatialConfig:
    """Validated live-FBA transport and adaptive-integration settings."""

    environment_volume_l: float
    substrate_layers: Tuple[SubstrateLayerSpec, ...]
    low_substrate_threshold_m: float = 1.0e-4
    depletion_safety: float = 0.9
    max_substeps: int = 256
    cache_enabled: bool = True
    conservation_atol_m: float = 1.0e-15

    def __post_init__(self) -> None:
        numeric_values = (
            ("environment_volume_l", self.environment_volume_l, True),
            ("low_substrate_threshold_m", self.low_substrate_threshold_m, False),
            ("depletion_safety", self.depletion_safety, True),
            ("conservation_atol_m", self.conservation_atol_m, False),
        )
        for name, value, strictly_positive in numeric_values:
            if isinstance(value, bool) or not isinstance(value, Real):
                raise TypeError(f"{name} must be a real number")
            numeric = float(value)
            if not math.isfinite(numeric) or numeric < 0.0 or (
                strictly_positive and numeric == 0.0
            ):
                relation = "positive" if strictly_positive else "nonnegative"
                raise ValueError(f"{name} must be finite and {relation}")
        if self.depletion_safety > 1.0:
            raise ValueError("depletion_safety must be in (0, 1]")
        if isinstance(self.max_substeps, bool) or not isinstance(self.max_substeps, Integral):
            raise TypeError("max_substeps must be an integer")
        if int(self.max_substeps) <= 0:
            raise ValueError("max_substeps must be positive")
        if not isinstance(self.cache_enabled, bool):
            raise TypeError("cache_enabled must be a boolean")
        if tuple(spec.name for spec in self.substrate_layers) != SUBSTRATE_ORDER:
            raise ValueError(f"substrate layer order must exactly match {SUBSTRATE_ORDER}")


@dataclass(frozen=True)
class FbaAgentState:
    """Committed reduced-metabolism state kept separate from legacy pg reserves."""

    biomass_gdw: float
    growth_rate_s: float = 0.0
    biomass_delta_gdw: float = 0.0
    substrate_consumed_mol: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    acetate_secreted_mol: float = 0.0


@dataclass(frozen=True)
class FbaPhysiology:
    """Lifecycle and rendering decisions derived from committed FBA state."""

    render_scale: float
    motility_scale: float
    chemotactic_drive: float
    division_probability: float
    starved: bool


@dataclass(frozen=True)
class GrowthOnlyFbaLifecyclePolicy:
    """Own-model biomass-doubling policy with no maintenance or starvation claim.

    The reduced LP can add biomass but has no maintenance, damage, or death reaction.
    Consequently this policy permits binary fission after solver-owned biomass doubles
    and intentionally never classifies an agent as starved.
    """

    initial_biomass_gdw: float
    evidence_class: str = "Own (BacSandbox model)"

    def __post_init__(self) -> None:
        if (
            isinstance(self.initial_biomass_gdw, bool)
            or not isinstance(self.initial_biomass_gdw, Real)
        ):
            raise TypeError("initial_biomass_gdw must be a real number")
        if not math.isfinite(float(self.initial_biomass_gdw)) or self.initial_biomass_gdw <= 0.0:
            raise ValueError("initial_biomass_gdw must be finite and positive")
        if self.evidence_class != "Own (BacSandbox model)":
            raise ValueError("growth-only FBA lifecycle must retain its Own-model evidence class")

    @property
    def division_biomass_gdw(self) -> float:
        return 2.0 * float(self.initial_biomass_gdw)

    def initialize_agent(self, bacterium: Any) -> FbaAgentState:
        return FbaAgentState(biomass_gdw=float(self.initial_biomass_gdw))

    def evaluate(
        self,
        bacterium: Any,
        state: FbaAgentState,
        *,
        dt_seconds: float,
        environmental_fitness: float,
        carrying_capacity_factor: float,
    ) -> FbaPhysiology:
        numeric_inputs = (
            ("dt_seconds", dt_seconds),
            ("environmental_fitness", environmental_fitness),
            ("carrying_capacity_factor", carrying_capacity_factor),
        )
        for name, value in numeric_inputs:
            if isinstance(value, bool) or not isinstance(value, Real):
                raise TypeError(f"{name} must be a real number")
            if not math.isfinite(float(value)) or float(value) < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")
        division_ready = (
            state.biomass_gdw >= self.division_biomass_gdw
            and float(carrying_capacity_factor) > 0.0
        )
        return FbaPhysiology(
            render_scale=1.0,
            motility_scale=1.0,
            chemotactic_drive=1.0,
            division_probability=1.0 if division_ready else 0.0,
            starved=False,
        )

    def partition_for_fission(
        self,
        parent: Any,
        child: Any,
        state: FbaAgentState,
    ) -> Tuple[FbaAgentState, FbaAgentState]:
        partitioned = FbaAgentState(
            biomass_gdw=state.biomass_gdw * 0.5,
            growth_rate_s=state.growth_rate_s,
            biomass_delta_gdw=state.biomass_delta_gdw * 0.5,
            substrate_consumed_mol=tuple(
                value * 0.5 for value in state.substrate_consumed_mol
            ),
            acetate_secreted_mol=state.acetate_secreted_mol * 0.5,
        )
        return partitioned, partitioned


class FbaAgentStateStore:
    """Explicit sidecar ownership compatible with pooled ``Bacterium`` objects."""

    def __init__(self, initializer: Callable[[Any], FbaAgentState]) -> None:
        if not callable(initializer):
            raise TypeError("initializer must be callable")
        self._initializer = initializer
        self._states: Dict[Any, FbaAgentState] = {}

    def reconcile(self, agents: Sequence[Any]) -> None:
        agents = tuple(agents)
        if len(set(agents)) != len(agents):
            raise ValueError("live FBA agent order contains duplicate objects")
        live = set(agents)
        candidate = {agent: state for agent, state in self._states.items() if agent in live}
        for agent in agents:
            if agent not in candidate:
                candidate[agent] = self._validate_state(self._initializer(agent))
        self._states = candidate

    def state_for(self, agent: Any) -> FbaAgentState:
        try:
            return self._states[agent]
        except KeyError as error:
            raise KeyError("agent is not registered in the FBA state store") from error

    def remove(self, agent: Any) -> None:
        self._states.pop(agent, None)

    def clear(self) -> None:
        """Remove all committed sidecar states for a fresh simulation run."""
        self._states.clear()

    def partition_for_fission(
        self,
        parent: Any,
        child: Any,
        lifecycle_policy: FbaLifecyclePolicy,
    ) -> None:
        if parent not in self._states:
            raise KeyError("parent is not registered in the FBA state store")
        parent_state, child_state = lifecycle_policy.partition_for_fission(
            parent,
            child,
            self._states[parent],
        )
        parent_state = self._validate_state(parent_state)
        child_state = self._validate_state(child_state)
        candidate = dict(self._states)
        candidate[parent] = parent_state
        candidate[child] = child_state
        self._states = candidate

    def snapshot(self) -> Dict[Any, FbaAgentState]:
        return dict(self._states)

    def restore(self, snapshot: Mapping[Any, FbaAgentState]) -> None:
        restored = {
            agent: self._validate_state(state)
            for agent, state in snapshot.items()
        }
        self._states = restored

    def commit_batch(self, batch: "SpatialMetabolicBatch") -> None:
        candidate = dict(self._states)
        for agent_index, agent in enumerate(batch.agents):
            state = FbaAgentState(
                biomass_gdw=float(batch.biomass_gdw[agent_index]),
                growth_rate_s=float(batch.growth_rate_s[agent_index]),
                biomass_delta_gdw=float(batch.biomass_delta_gdw[agent_index]),
                substrate_consumed_mol=tuple(
                    float(value) for value in batch.substrate_consumed_mol[agent_index]
                ),
                acetate_secreted_mol=float(batch.acetate_secreted_mol[agent_index]),
            )
            candidate[agent] = self._validate_state(state)
        self._states = candidate

    @staticmethod
    def _validate_state(state: FbaAgentState) -> FbaAgentState:
        if not isinstance(state, FbaAgentState):
            raise TypeError("FBA initializer must return FbaAgentState")
        values = (
            state.biomass_gdw,
            state.growth_rate_s,
            state.biomass_delta_gdw,
            *state.substrate_consumed_mol,
            state.acetate_secreted_mol,
        )
        if len(state.substrate_consumed_mol) != len(SUBSTRATE_ORDER):
            raise ValueError("substrate_consumed_mol must contain three values")
        if not all(math.isfinite(float(value)) for value in values):
            raise ValueError("FBA state values must be finite")
        nonnegative = (
            state.biomass_gdw,
            state.biomass_delta_gdw,
            *state.substrate_consumed_mol,
            state.acetate_secreted_mol,
        )
        if any(float(value) < 0.0 for value in nonnegative):
            raise ValueError("FBA amount and biomass state must be nonnegative")
        return state


@dataclass
class SpatialMetabolicBatch:
    """Caller-owned arrays matching ``MetabolicOrchestrator.advance_batch``."""

    agents: Tuple[Any, ...]
    coordinates: Tuple[GridCoordinate, ...]
    phenotype_id: np.ndarray
    environment_id: np.ndarray
    concentrations_m: np.ndarray
    environment_volume_l: np.ndarray
    biomass_gdw: np.ndarray
    growth_rate_s: np.ndarray
    biomass_delta_gdw: np.ndarray
    substrate_consumed_mol: np.ndarray
    acetate_secreted_mol: np.ndarray


@dataclass(frozen=True)
class DietBolus:
    """One dry-feed bolus distributed over normalized grid-cell weights."""

    event_id: str
    diet_name: str
    feed_mass_g: float
    footprint: Tuple[Tuple[GridCoordinate, float], ...]


@dataclass(frozen=True)
class DietBolusDeposited:
    event_id: str
    total_substrate_mol: Tuple[float, float, float]
    affected_coordinates: Tuple[GridCoordinate, ...]


class FbaTickDisposition(str, Enum):
    COMMITTED = "committed"
    NO_AGENTS = "no_agents"
    ROLLED_BACK = "rolled_back"


@dataclass(frozen=True)
class FbaTickResult:
    disposition: FbaTickDisposition
    tick_index: int
    dt_seconds: float
    backend_used: Optional[str]
    agent_count: int
    environment_count: int
    kernel_summary: Optional[Tuple[int, int, int, int, int, float]]
    substrate_consumed_mol: Tuple[float, float, float]
    biomass_delta_gdw: float
    acetate_secreted_mol: float
    max_conservation_error_mol: float


@dataclass(frozen=True)
class FbaTransactionEvent:
    disposition: FbaTickDisposition
    tick_index: int
    phase: str
    backend_used: Optional[str]
    error_type: Optional[str]
    error_message: Optional[str]


def deposit_diet_bolus(
    engine: DietTranslationEngine,
    ecm: ECM,
    bolus: DietBolus,
    *,
    environment_volume_l: float,
    dt_seconds: float,
) -> DietBolusDeposited:
    """Deposit a bolus while conserving translated amount across its footprint."""
    if not isinstance(bolus, DietBolus):
        raise TypeError("bolus must be a DietBolus")
    if not isinstance(bolus.event_id, str) or not bolus.event_id.strip():
        raise ValueError("bolus event_id must be a non-empty string")
    if isinstance(environment_volume_l, bool) or not isinstance(environment_volume_l, Real):
        raise TypeError("environment_volume_l must be a real number")
    environment_volume_l = float(environment_volume_l)
    if not math.isfinite(environment_volume_l) or environment_volume_l <= 0.0:
        raise ValueError("environment_volume_l must be finite and positive")
    if not bolus.footprint:
        raise ValueError("bolus footprint must contain at least one grid cell")

    coordinates = tuple(coordinate for coordinate, _ in bolus.footprint)
    ecm.sample_substrate_batch_m(coordinates)
    weights = []
    for _, weight in bolus.footprint:
        if isinstance(weight, bool) or not isinstance(weight, Real):
            raise TypeError("bolus footprint weights must be real numbers")
        weight = float(weight)
        if not math.isfinite(weight) or weight <= 0.0:
            raise ValueError("bolus footprint weights must be finite and positive")
        weights.append(weight)
    if not math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("bolus footprint weights must sum to one")

    translation_volume_l = environment_volume_l * len(coordinates)
    concentration_increment = engine.get_molar_substrate_bounds(
        bolus.diet_name,
        bolus.feed_mass_g,
        translation_volume_l,
        dt_seconds,
    )
    total_substrate_mol = concentration_increment * translation_volume_l
    existing = ecm.sample_substrate_batch_m(coordinates)
    increments = np.asarray(weights, dtype=np.float64)[:, None] * (
        total_substrate_mol[None, :] / environment_volume_l
    )
    updated = np.array(existing + increments, dtype=np.float64, order="C", copy=True)
    ecm.replace_substrate_batch_m(coordinates, updated)
    return DietBolusDeposited(
        event_id=bolus.event_id,
        total_substrate_mol=tuple(float(value) for value in total_substrate_mol),
        affected_coordinates=coordinates,
    )


def build_coordinate_environment_batch(
    species: Sequence[Any],
    ecm: ECM,
    states: FbaAgentStateStore,
    phenotype_index: Mapping[Any, int],
    environment_volume_l: float,
) -> SpatialMetabolicBatch:
    """Build deterministic agent and environment rows from live coordinates."""
    agents = tuple(agent for profile in species for agent in profile.members)
    states.reconcile(agents)
    agent_coordinates = tuple(ecm.world_to_cell(agent.x, agent.y) for agent in agents)
    coordinates = tuple(sorted(set(agent_coordinates)))
    environment_lookup = {coordinate: index for index, coordinate in enumerate(coordinates)}
    environment_id = np.asarray(
        tuple(environment_lookup[coordinate] for coordinate in agent_coordinates),
        dtype=np.int32,
    )
    phenotype_values = []
    for profile in species:
        if profile not in phenotype_index:
            raise KeyError("every live species must have an FBA phenotype index")
        value = phenotype_index[profile]
        if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
            raise ValueError("phenotype indices must be nonnegative integers")
        phenotype_values.extend([int(value)] * len(profile.members))
    phenotype_id = np.asarray(phenotype_values, dtype=np.int32)
    concentrations_m = ecm.sample_substrate_batch_m(coordinates)
    environment_volume = np.full(len(coordinates), environment_volume_l, dtype=np.float64)
    biomass = np.asarray(
        tuple(states.state_for(agent).biomass_gdw for agent in agents),
        dtype=np.float64,
    )
    agent_count = len(agents)
    return SpatialMetabolicBatch(
        agents=agents,
        coordinates=coordinates,
        phenotype_id=phenotype_id,
        environment_id=environment_id,
        concentrations_m=concentrations_m,
        environment_volume_l=environment_volume,
        biomass_gdw=biomass,
        growth_rate_s=np.zeros(agent_count, dtype=np.float64),
        biomass_delta_gdw=np.zeros(agent_count, dtype=np.float64),
        substrate_consumed_mol=np.zeros((agent_count, len(SUBSTRATE_ORDER)), dtype=np.float64),
        acetate_secreted_mol=np.zeros(agent_count, dtype=np.float64),
    )


class FbaIntegrationAdapter:
    """Own one atomic live-FBA transaction across grid and sidecar state."""

    def __init__(
        self,
        orchestrator: MetabolicOrchestrator,
        diet_engine: DietTranslationEngine,
        config: FbaSpatialConfig,
        phenotype_index: Mapping[Any, int],
        lifecycle_policy: FbaLifecyclePolicy,
        *,
        logger: Optional[logging.Logger] = None,
        event_sink: Optional[Callable[[FbaTransactionEvent], None]] = None,
    ) -> None:
        if getattr(orchestrator, "mode", None) != "fba":
            raise ValueError("live integration requires a MetabolicOrchestrator in fba mode")
        if not isinstance(diet_engine, DietTranslationEngine):
            raise TypeError("diet_engine must be a DietTranslationEngine")
        if not isinstance(config, FbaSpatialConfig):
            raise TypeError("config must be FbaSpatialConfig")
        if not hasattr(lifecycle_policy, "initialize_agent"):
            raise TypeError("lifecycle_policy must define initialize_agent")
        if event_sink is not None and not callable(event_sink):
            raise TypeError("event_sink must be callable")
        self.orchestrator = orchestrator
        self.diet_engine = diet_engine
        self.config = config
        self.phenotype_index = dict(phenotype_index)
        self.lifecycle_policy = lifecycle_policy
        self.state_store = FbaAgentStateStore(lifecycle_policy.initialize_agent)
        self.logger = logger or logging.getLogger(__name__)
        self.event_sink = event_sink

    def advance_tick(
        self,
        *,
        species: Sequence[Any],
        ecm: ECM,
        tick_index: int,
        dt_seconds: float,
        boluses: Sequence[DietBolus] = (),
        defer_events: bool = False,
    ) -> FbaTickResult:
        if isinstance(tick_index, bool) or not isinstance(tick_index, Integral) or tick_index < 0:
            raise ValueError("tick_index must be a nonnegative integer")
        if isinstance(dt_seconds, bool) or not isinstance(dt_seconds, Real):
            raise TypeError("dt_seconds must be a real number")
        dt_seconds = float(dt_seconds)
        if not math.isfinite(dt_seconds) or dt_seconds < 0.0:
            raise ValueError("dt_seconds must be finite and nonnegative")

        substrate_snapshot: SubstrateLayerSnapshot = ecm.snapshot_substrate_layers()
        state_snapshot = self.state_store.snapshot()
        cache_snapshot = self.orchestrator.snapshot_cache()
        phase = "registration"
        try:
            ecm.register_substrate_layers(self.config.substrate_layers)
            phase = "reconcile"
            agents = tuple(agent for profile in species for agent in profile.members)
            self.state_store.reconcile(agents)
            phase = "bolus"
            for bolus in boluses:
                deposit_diet_bolus(
                    self.diet_engine,
                    ecm,
                    bolus,
                    environment_volume_l=self.config.environment_volume_l,
                    dt_seconds=dt_seconds,
                )
            phase = "diffusion"
            ecm.tick_substrates(dt_seconds)
            phase = "batch"
            batch = build_coordinate_environment_batch(
                species,
                ecm,
                self.state_store,
                self.phenotype_index,
                self.config.environment_volume_l,
            )
            if not batch.agents:
                result = self._result(
                    FbaTickDisposition.NO_AGENTS,
                    int(tick_index),
                    dt_seconds,
                    batch,
                    None,
                    0.0,
                )
                if not defer_events:
                    self.emit_success(result)
                return result

            initial_concentrations = batch.concentrations_m.copy()
            initial_biomass = batch.biomass_gdw.copy()
            phase = "solve"
            summary = self.orchestrator.advance_batch(
                phenotype_id=batch.phenotype_id,
                environment_id=batch.environment_id,
                concentrations_m=batch.concentrations_m,
                environment_volume_l=batch.environment_volume_l,
                biomass_gdw=batch.biomass_gdw,
                growth_rate_s=batch.growth_rate_s,
                biomass_delta_gdw=batch.biomass_delta_gdw,
                substrate_consumed_mol=batch.substrate_consumed_mol,
                acetate_secreted_mol=batch.acetate_secreted_mol,
                dt_seconds=dt_seconds,
                low_substrate_threshold_m=self.config.low_substrate_threshold_m,
                depletion_safety=self.config.depletion_safety,
                max_substeps=int(self.config.max_substeps),
                cache_enabled=self.config.cache_enabled,
            )
            phase = "validate"
            max_error = self._validate_batch(
                batch,
                initial_concentrations,
                initial_biomass,
            )
            phase = "commit"
            ecm.replace_substrate_batch_m(batch.coordinates, batch.concentrations_m)
            self.state_store.commit_batch(batch)
            result = self._result(
                FbaTickDisposition.COMMITTED,
                int(tick_index),
                dt_seconds,
                batch,
                tuple(summary),
                max_error,
            )
            if not defer_events:
                self.emit_success(result)
            return result
        except Exception as error:
            ecm.restore_substrate_layers(substrate_snapshot)
            self.state_store.restore(state_snapshot)
            self.orchestrator.restore_cache(cache_snapshot)
            self.logger.exception("live FBA tick rolled back at phase %s", phase)
            if not defer_events:
                self.emit_failure(int(tick_index), phase, error)
            raise

    def emit_success(self, result: FbaTickResult) -> None:
        """Emit one success event after the owning transaction has fully committed."""
        self._emit(FbaTransactionEvent(
            disposition=result.disposition,
            tick_index=result.tick_index,
            phase="commit",
            backend_used=result.backend_used,
            error_type=None,
            error_message=None,
        ))

    def emit_failure(self, tick_index: int, phase: str, error: BaseException) -> None:
        """Emit one rollback event without retaining the exception or traceback graph."""
        self._emit(FbaTransactionEvent(
            disposition=FbaTickDisposition.ROLLED_BACK,
            tick_index=int(tick_index),
            phase=phase,
            backend_used=getattr(self.orchestrator, "backend_used", None),
            error_type=type(error).__name__,
            error_message=str(error),
        ))

    def _validate_batch(
        self,
        batch: SpatialMetabolicBatch,
        initial_concentrations: np.ndarray,
        initial_biomass: np.ndarray,
    ) -> float:
        arrays = (
            batch.concentrations_m,
            batch.biomass_gdw,
            batch.growth_rate_s,
            batch.biomass_delta_gdw,
            batch.substrate_consumed_mol,
            batch.acetate_secreted_mol,
        )
        if not all(np.all(np.isfinite(array)) for array in arrays):
            raise ValueError("live FBA result arrays must contain only finite values")
        for array in (
            batch.concentrations_m,
            batch.biomass_gdw,
            batch.biomass_delta_gdw,
            batch.substrate_consumed_mol,
            batch.acetate_secreted_mol,
        ):
            if np.any(array < 0.0):
                raise ValueError("live FBA amount and biomass arrays must be nonnegative")
        if np.any(batch.concentrations_m > initial_concentrations):
            raise FbaConservationError("FBA solve increased a substrate concentration")
        biomass_errors = np.abs(
            (batch.biomass_gdw - initial_biomass) - batch.biomass_delta_gdw
        )
        if np.any(biomass_errors > BIOMASS_CONSERVATION_ATOL_GDW):
            maximum = float(np.max(biomass_errors)) if biomass_errors.size else 0.0
            raise FbaConservationError(
                f"FBA solve biomass delta mismatch; maximum error={maximum!r} gDW"
            )

        consumed_by_environment = np.zeros_like(batch.concentrations_m)
        np.add.at(
            consumed_by_environment,
            batch.environment_id,
            batch.substrate_consumed_mol,
        )
        depleted_mol = (
            initial_concentrations - batch.concentrations_m
        ) * batch.environment_volume_l[:, None]
        errors = np.abs(depleted_mol - consumed_by_environment)
        allowed = self.config.conservation_atol_m * batch.environment_volume_l[:, None]
        if np.any(errors > allowed):
            maximum = float(np.max(errors)) if errors.size else 0.0
            raise FbaConservationError(
                f"represented substrate mass mismatch; maximum error={maximum!r} mol"
            )
        return float(np.max(errors)) if errors.size else 0.0

    def _result(
        self,
        disposition: FbaTickDisposition,
        tick_index: int,
        dt_seconds: float,
        batch: SpatialMetabolicBatch,
        summary: Optional[Tuple[int, int, int, int, int, float]],
        max_error: float,
    ) -> FbaTickResult:
        consumed = tuple(
            float(value)
            for value in np.sum(batch.substrate_consumed_mol, axis=0)
        )
        return FbaTickResult(
            disposition=disposition,
            tick_index=tick_index,
            dt_seconds=dt_seconds,
            backend_used=getattr(self.orchestrator, "backend_used", None),
            agent_count=len(batch.agents),
            environment_count=len(batch.coordinates),
            kernel_summary=summary,
            substrate_consumed_mol=consumed,
            biomass_delta_gdw=float(np.sum(batch.biomass_delta_gdw)),
            acetate_secreted_mol=float(np.sum(batch.acetate_secreted_mol)),
            max_conservation_error_mol=max_error,
        )

    def _emit(self, event: FbaTransactionEvent) -> None:
        if self.event_sink is None:
            return
        try:
            self.event_sink(event)
        except Exception:
            self.logger.exception("live FBA event sink failed")