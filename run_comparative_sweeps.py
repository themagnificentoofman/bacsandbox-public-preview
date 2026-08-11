#!/usr/bin/env python3
"""Headless paired-seed sweeps for the implemented legacy pg production model."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator, Mapping, Sequence, Tuple

os.environ.setdefault("PYGAME_HIDE_SUPPORT_PROMPT", "1")

import numpy as np

import bacteria.core_sim as core_sim
from bacteria.constants import (
    ATMOSPHERIC_O2_FRACTION,
    DEFAULT_AFF_CALORIE,
    DEFAULT_AFF_GLYCAN,
    DEFAULT_AFF_MDCF,
    DEFAULT_BARANYI_Q0,
    DEFAULT_ENV_TEMPERATURE_C,
    DEFAULT_HGT_COMPETENCE,
    DEFAULT_PERSISTER_ENTER_THRESH,
    DEFAULT_PERSISTER_EXIT_THRESH,
    DEFAULT_PERSISTER_SWITCH,
    DEFAULT_PHAGE_RESISTANCE,
    DIET_PROFILES,
    ECM_CELL_SIZE_UM,
    INITIAL_SPAWN_CENTRAL_RADIUS_FRACTION,
    INITIAL_SPAWN_CLUSTER_COUNT,
    INITIAL_SPAWN_CLUSTER_SPREAD_FRACTION,
    INITIAL_SPAWN_LUMEN_BIAS,
    INITIAL_SPAWN_LUMEN_RADIUS_FRACTION,
    INITIAL_SPAWN_MODE,
    INITIAL_SPECIES_ARRANGEMENT,
    LUMEN_COMPETITION_CULL_WEIGHT,
    LUMEN_FOOD_WASHOUT_RATE_PER_SEC,
    MUCUS_DROP_AMOUNT_FACTOR,
    MUCUS_FOOD_WASHOUT_RATE_PER_SEC,
    MUCUS_LAYER_RADIUS_FRACTION,
    MUCUS_MOTILITY_FACTOR,
    NEUTRAL_PH,
    PETRI_RADIUS_PADDING_UM,
)
from bacteria.core_sim import BacteriaSpecies, Simulation
from bacteria.data import SpeciesDataManager
from bacteria.metabolism import MODEL_ID, MODEL_REVISION, PRIOR_ID, PRIOR_REVISION


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_SPECIES_DATA_PATH = PROJECT_ROOT / "species_data.csv"
DEFAULT_OUTPUT_DIRECTORY = Path("/workspace/scratch")
CONDITION_HOMOGENEOUS = "homogeneous"
CONDITION_SPATIAL_NICHE = "spatial_niche"
CONDITION_ORDER = (CONDITION_HOMOGENEOUS, CONDITION_SPATIAL_NICHE)
SPATIAL_EVIDENCE_CLASS = "Own (BacSandbox model)"
LOW_ABUNDANCE_EVIDENCE_CLASS = "Own (BacSandbox analysis definition)"
FBA_LIMITATION = (
    "The reduced live FBA LP has no maintenance, damage, or starvation reaction. "
    "These extinction records come only from the implemented legacy pg production model "
    "and do not validate reduced-FBA extinction dynamics."
)


RELATIVE_ABUNDANCE_FIELDS = (
    "sweep_id",
    "pair_id",
    "seed",
    "diet",
    "condition",
    "step",
    "elapsed_seconds",
    "species",
    "abundance_count",
    "total_count",
    "relative_abundance",
)
JACCARD_FIELDS = (
    "sweep_id",
    "pair_id",
    "seed",
    "diet",
    "step",
    "elapsed_seconds",
    "left_condition",
    "right_condition",
    "left_present_species_count",
    "right_present_species_count",
    "intersection_size",
    "union_size",
    "jaccard_community_distance",
)
EXTINCTION_FIELDS = (
    "sweep_id",
    "pair_id",
    "seed",
    "diet",
    "condition",
    "species",
    "initial_count",
    "final_count",
    "low_abundance_max_initial_count",
    "low_abundance_definition",
    "low_abundance_evidence_class",
    "event_observed",
    "status",
    "time_to_event_seconds",
    "censor_time_seconds",
)
METADATA_FIELDS = (
    "sweep_id",
    "pair_id",
    "seed",
    "diet",
    "condition",
    "dt_seconds",
    "step_count",
    "timepoint_count",
    "simulation_engine",
    "backend",
    "physics_backend_used",
    "ecm_backend_used",
    "metabolic_model_id",
    "metabolic_model_revision",
    "metabolic_prior_id",
    "metabolic_prior_revision",
    "reduced_model_active",
    "available_fba_lifecycle_policy",
    "fba_extinction_dynamics_validated",
    "fba_lifecycle_limitation",
    "species_data_path",
    "species_data_sha256",
    "founder_state_sha256",
    "post_founder_rng_state_sha256",
    "world_width_um",
    "world_height_um",
    "carrying_capacity_cells",
    "drop_frequency_s_inv",
    "drop_range_min_pg",
    "drop_range_max_pg",
    "low_abundance_max_initial_count",
    "low_abundance_definition",
    "low_abundance_evidence_class",
    "spatial_parameter_evidence_class",
    "spatial_partition_enabled",
    "mucus_layer_radius_fraction",
    "mucus_motility_anchoring_factor",
    "lumen_food_washout_rate_s_inv",
    "mucus_food_washout_rate_s_inv",
    "mucus_drop_amount_factor",
    "lumen_competition_cull_weight",
    "initial_spawn_mode",
    "initial_species_arrangement",
    "initial_spawn_lumen_bias",
    "initial_spawn_lumen_radius_fraction",
    "initial_spawn_cluster_count",
    "initial_spawn_cluster_spread_fraction",
    "initial_spawn_central_radius_fraction",
    "petri_radius_padding_um",
    "ecm_cell_size_um",
)
CSV_SCHEMAS = {
    "relative_abundance.csv": RELATIVE_ABUNDANCE_FIELDS,
    "jaccard_distance.csv": JACCARD_FIELDS,
    "extinction_censoring.csv": EXTINCTION_FIELDS,
    "run_metadata.csv": METADATA_FIELDS,
}


@dataclass(frozen=True)
class SpatialControl:
    condition: str
    partition_enabled: bool
    mucus_layer_radius_fraction: float
    mucus_motility_factor: float
    lumen_food_washout_rate_s_inv: float
    mucus_food_washout_rate_s_inv: float
    mucus_drop_amount_factor: float
    lumen_competition_cull_weight: float

    def __post_init__(self) -> None:
        if self.condition not in CONDITION_ORDER:
            raise ValueError(f"unsupported condition: {self.condition}")
        values = (
            self.mucus_layer_radius_fraction,
            self.mucus_motility_factor,
            self.lumen_food_washout_rate_s_inv,
            self.mucus_food_washout_rate_s_inv,
            self.mucus_drop_amount_factor,
            self.lumen_competition_cull_weight,
        )
        if not all(math.isfinite(value) and value >= 0.0 for value in values):
            raise ValueError("spatial control values must be finite and nonnegative")


HOMOGENEOUS_CONTROL = SpatialControl(
    condition=CONDITION_HOMOGENEOUS,
    partition_enabled=False,
    mucus_layer_radius_fraction=1.0,
    mucus_motility_factor=1.0,
    lumen_food_washout_rate_s_inv=LUMEN_FOOD_WASHOUT_RATE_PER_SEC,
    mucus_food_washout_rate_s_inv=LUMEN_FOOD_WASHOUT_RATE_PER_SEC,
    mucus_drop_amount_factor=1.0,
    lumen_competition_cull_weight=1.0,
)
SPATIAL_NICHE_CONTROL = SpatialControl(
    condition=CONDITION_SPATIAL_NICHE,
    partition_enabled=True,
    mucus_layer_radius_fraction=MUCUS_LAYER_RADIUS_FRACTION,
    mucus_motility_factor=MUCUS_MOTILITY_FACTOR,
    lumen_food_washout_rate_s_inv=LUMEN_FOOD_WASHOUT_RATE_PER_SEC,
    mucus_food_washout_rate_s_inv=MUCUS_FOOD_WASHOUT_RATE_PER_SEC,
    mucus_drop_amount_factor=MUCUS_DROP_AMOUNT_FACTOR,
    lumen_competition_cull_weight=LUMEN_COMPETITION_CULL_WEIGHT,
)
SPATIAL_CONTROLS = {
    control.condition: control
    for control in (HOMOGENEOUS_CONTROL, SPATIAL_NICHE_CONTROL)
}


@dataclass(frozen=True)
class SweepConfig:
    species_data_path: Path = DEFAULT_SPECIES_DATA_PATH
    output_directory: Path = DEFAULT_OUTPUT_DIRECTORY
    seeds: Tuple[int, ...] = (0,)
    diets: Tuple[str, ...] = tuple(DIET_PROFILES)
    steps: int = 1_000
    dt_seconds: float = 0.1
    world_width_um: int = 800
    world_height_um: int = 600
    carrying_capacity_cells: int = 500
    drop_frequency_s_inv: float = 1.0
    drop_range_pg: Tuple[float, float] = (1.0, 50.0)
    backend: str = "python"
    low_abundance_max_initial_count: int = 10
    sweep_id: str = "comparative_sweep"

    def __post_init__(self) -> None:
        if not self.seeds or any(isinstance(seed, bool) or seed < 0 for seed in self.seeds):
            raise ValueError("seeds must contain nonnegative integers")
        if not self.diets or any(diet not in DIET_PROFILES for diet in self.diets):
            raise ValueError(f"diets must be selected from {tuple(DIET_PROFILES)}")
        if isinstance(self.steps, bool) or self.steps <= 0:
            raise ValueError("steps must be a positive integer")
        if not math.isfinite(self.dt_seconds) or self.dt_seconds <= 0.0:
            raise ValueError("dt_seconds must be finite and positive")
        if self.world_width_um <= 0 or self.world_height_um <= 0:
            raise ValueError("world dimensions must be positive")
        if self.carrying_capacity_cells <= 0:
            raise ValueError("carrying_capacity_cells must be positive")
        if not math.isfinite(self.drop_frequency_s_inv) or self.drop_frequency_s_inv < 0.0:
            raise ValueError("drop_frequency_s_inv must be finite and nonnegative")
        if (
            len(self.drop_range_pg) != 2
            or not all(math.isfinite(value) and value >= 0.0 for value in self.drop_range_pg)
            or self.drop_range_pg[0] > self.drop_range_pg[1]
        ):
            raise ValueError("drop_range_pg must be a finite nonnegative ordered pair")
        if self.backend not in {"auto", "python", "cpp"}:
            raise ValueError("backend must be auto, python, or cpp")
        if self.low_abundance_max_initial_count <= 0:
            raise ValueError("low_abundance_max_initial_count must be positive")
        if not self.sweep_id.strip():
            raise ValueError("sweep_id must be nonempty")


@dataclass(frozen=True)
class ConditionRun:
    condition: str
    seed: int
    diet: str
    founder_state_sha256: str
    post_founder_rng_state_sha256: str
    elapsed_seconds: Tuple[float, ...]
    population_counts: Tuple[Mapping[str, int], ...]
    initial_counts: Mapping[str, int]
    extinction_times: Mapping[str, float | None]
    physics_backend_used: str
    ecm_backend_used: str


@dataclass(frozen=True)
class SweepResult:
    sweep_id: str
    relative_abundance_rows: Tuple[Mapping[str, object], ...]
    jaccard_rows: Tuple[Mapping[str, object], ...]
    extinction_rows: Tuple[Mapping[str, object], ...]
    metadata_rows: Tuple[Mapping[str, object], ...]


def jaccard_community_distance(
    left_counts: Mapping[str, int],
    right_counts: Mapping[str, int],
) -> float:
    """Return presence/absence Jaccard distance; two empty communities have distance zero."""
    left_present = {name for name, count in left_counts.items() if count > 0}
    right_present = {name for name, count in right_counts.items() if count > 0}
    union = left_present | right_present
    if not union:
        return 0.0
    return 1.0 - len(left_present & right_present) / len(union)


@contextmanager
def isolated_rng(seed: int) -> Iterator[None]:
    """Run against seeded global RNGs, then restore both caller-owned states exactly."""
    python_state = random.getstate()
    numpy_state = np.random.get_state()
    random.seed(seed)
    np.random.seed(seed % (2**32))
    try:
        yield
    finally:
        random.setstate(python_state)
        np.random.set_state(numpy_state)


@contextmanager
def applied_spatial_control(control: SpatialControl) -> Iterator[None]:
    """Temporarily apply one sequential headless control to production-model globals."""
    replacements = {
        "MUCUS_LAYER_RADIUS_FRACTION": control.mucus_layer_radius_fraction,
        "MUCUS_MOTILITY_FACTOR": control.mucus_motility_factor,
        "LUMEN_FOOD_WASHOUT_RATE_PER_SEC": control.lumen_food_washout_rate_s_inv,
        "MUCUS_FOOD_WASHOUT_RATE_PER_SEC": control.mucus_food_washout_rate_s_inv,
        "MUCUS_DROP_AMOUNT_FACTOR": control.mucus_drop_amount_factor,
        "LUMEN_COMPETITION_CULL_WEIGHT": control.lumen_competition_cull_weight,
    }
    original = {name: getattr(core_sim, name) for name in replacements}
    try:
        for name, value in replacements.items():
            setattr(core_sim, name, value)
        yield
    finally:
        for name, value in original.items():
            setattr(core_sim, name, value)


def _row_value(row: Mapping[str, str], field: str, default: float) -> float:
    value = row.get(field)
    return float(default if value is None or value == "" else value)


def load_species(species_data_path: Path) -> list[BacteriaSpecies]:
    """Load the same parameter surface as the UI without importing presentation code."""
    manager = SpeciesDataManager(str(species_data_path))
    colors = ((200, 50, 50), (50, 200, 50), (50, 50, 200), (200, 200, 50))
    profiles = []
    for index, row in enumerate(manager.species):
        name = row.get("name", "").strip()
        if not name:
            raise ValueError(f"species row {index + 2} has an empty name")
        profiles.append(BacteriaSpecies(
            name=name,
            init=int(_row_value(row, "init", 0.0)),
            rep=_row_value(row, "rep", 0.0),
            aging=_row_value(row, "aging", 0.0),
            so=_row_value(row, "so", 0.0),
            sf=_row_value(row, "sf", 0.0),
            speed=_row_value(row, "speed", 40.0),
            ph=_row_value(row, "ph", 7.0),
            color=colors[index % len(colors)],
            food_capacity=_row_value(row, "food_capacity", 10.0),
            metabolic_rate=_row_value(row, "metabolic_rate", 0.1),
            consumption_speed=_row_value(row, "consumption_speed", 1.0),
            greediness=_row_value(row, "greediness", 0.5),
            t6ss_activity=_row_value(row, "t6ss_activity", 0.5),
            bacteriocin_production=_row_value(row, "bacteriocin_production", 0.5),
            immunity_factor=_row_value(row, "immunity_factor", 0.5),
            rep_thresh=_row_value(row, "rep_thresh", 0.8),
            share_thresh=_row_value(row, "share_thresh", 0.6),
            predation_thresh=_row_value(row, "predation_thresh", 0.9),
            base_size=_row_value(row, "base_size", 1.0),
            sensor_radius=_row_value(row, "sensor_radius", 200.0),
            temp_min=_row_value(row, "temp_min", 4.0),
            temp_max=_row_value(row, "temp_max", 47.0),
            acetate_secretion=_row_value(row, "acetate_secretion", 0.0),
            acetate_uptake=_row_value(row, "acetate_uptake", 0.0),
            baranyi_q0=_row_value(row, "baranyi_q0", DEFAULT_BARANYI_Q0),
            persister_switch=_row_value(row, "persister_switch", DEFAULT_PERSISTER_SWITCH),
            persister_enter_thresh=_row_value(
                row,
                "persister_enter_thresh",
                DEFAULT_PERSISTER_ENTER_THRESH,
            ),
            persister_exit_thresh=_row_value(
                row,
                "persister_exit_thresh",
                DEFAULT_PERSISTER_EXIT_THRESH,
            ),
            phage_resistance=_row_value(row, "phage_resistance", DEFAULT_PHAGE_RESISTANCE),
            hgt_competence=_row_value(row, "hgt_competence", DEFAULT_HGT_COMPETENCE),
            aff_calorie=_row_value(row, "aff_calorie", DEFAULT_AFF_CALORIE),
            aff_glycan=_row_value(row, "aff_glycan", DEFAULT_AFF_GLYCAN),
            aff_mdcf=_row_value(row, "aff_mdcf", DEFAULT_AFF_MDCF),
        ))
    names = [profile.name for profile in profiles]
    if len(names) != len(set(names)):
        raise ValueError("species names must be unique for community-distance output")
    if not profiles or not any(profile.init_count > 0 for profile in profiles):
        raise ValueError("species data must include at least one positive founder count")
    return profiles


def _founder_state_hash(simulation: Simulation) -> str:
    state = []
    fields = (
        "x",
        "y",
        "base_size",
        "food_capacity",
        "food_level",
        "metabolic_rate",
        "consumption_speed",
        "substrate_affinity",
        "speed",
        "size",
        "age",
        "lifespan_variance",
        "heading",
        "base_run_duration",
        "state_timer",
        "is_tumbling",
        "q",
        "is_persister",
        "phage_resistant",
        "infected",
        "lysis_timer",
        "metabolic_burden_factor",
    )
    for profile in simulation.species:
        agents = []
        for bacterium in profile.members:
            agents.append({field: getattr(bacterium, field) for field in fields})
        state.append({"species": profile.name, "agents": agents})
    payload = json.dumps(state, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _rng_state_hash() -> str:
    payload = repr((random.getstate(), np.random.get_state())).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _population_counts(simulation: Simulation) -> Dict[str, int]:
    return {profile.name: len(profile.members) for profile in simulation.species}


def _run_condition(config: SweepConfig, seed: int, diet: str, control: SpatialControl) -> ConditionRun:
    with isolated_rng(seed):
        species = load_species(config.species_data_path)
        environment = {
            "ph": NEUTRAL_PH,
            "o2": ATMOSPHERIC_O2_FRACTION,
            "co2": 0.0,
            "temperature": DEFAULT_ENV_TEMPERATURE_C,
            "drop_frequency": config.drop_frequency_s_inv,
            "drop_range": config.drop_range_pg,
            "phage_inflow": 0.0,
            "diet_profile": diet,
        }
        simulation = Simulation(
            species,
            environment,
            cap=config.carrying_capacity_cells,
            steps=config.steps,
            seed=seed,
            physics_backend=config.backend,
            ecm_backend=config.backend,
        )
        simulation._ensure_petri_geometry(config.world_width_um, config.world_height_um)
        simulation._seed_initial_population()
        founder_state_sha256 = _founder_state_hash(simulation)
        post_founder_rng_state_sha256 = _rng_state_hash()
        initial_counts = _population_counts(simulation)
        population_counts: list[Mapping[str, int]] = [initial_counts]
        elapsed_seconds = [0.0]
        simulation.running = True
        with applied_spatial_control(control):
            for _ in range(config.steps):
                simulation.update(
                    config.world_width_um,
                    config.world_height_um,
                    config.dt_seconds,
                )
                population_counts.append(_population_counts(simulation))
                elapsed_seconds.append(float(simulation.stats["elapsed_seconds"][-1]))
        extinction_times = {
            profile.name: profile.extinction_time
            for profile in simulation.species
        }
        return ConditionRun(
            condition=control.condition,
            seed=seed,
            diet=diet,
            founder_state_sha256=founder_state_sha256,
            post_founder_rng_state_sha256=post_founder_rng_state_sha256,
            elapsed_seconds=tuple(elapsed_seconds),
            population_counts=tuple(population_counts),
            initial_counts=dict(initial_counts),
            extinction_times=extinction_times,
            physics_backend_used=simulation.physics.backend_used or "not_exercised",
            ecm_backend_used=simulation.ecm.backend_used or "not_exercised",
        )


def _pair_id(seed: int, diet: str) -> str:
    return f"{diet}:seed={seed}"


def _relative_abundance_rows(
    sweep_id: str,
    pair_id: str,
    run: ConditionRun,
) -> Iterable[Mapping[str, object]]:
    for step, (elapsed, counts) in enumerate(zip(run.elapsed_seconds, run.population_counts)):
        total = sum(counts.values())
        for species, count in counts.items():
            yield {
                "sweep_id": sweep_id,
                "pair_id": pair_id,
                "seed": run.seed,
                "diet": run.diet,
                "condition": run.condition,
                "step": step,
                "elapsed_seconds": elapsed,
                "species": species,
                "abundance_count": count,
                "total_count": total,
                "relative_abundance": count / total if total else 0.0,
            }


def _jaccard_rows(
    sweep_id: str,
    pair_id: str,
    homogeneous: ConditionRun,
    spatial: ConditionRun,
) -> Iterable[Mapping[str, object]]:
    if homogeneous.elapsed_seconds != spatial.elapsed_seconds:
        raise RuntimeError("paired conditions produced different time grids")
    for step, (elapsed, left, right) in enumerate(zip(
        homogeneous.elapsed_seconds,
        homogeneous.population_counts,
        spatial.population_counts,
    )):
        left_present = {name for name, count in left.items() if count > 0}
        right_present = {name for name, count in right.items() if count > 0}
        intersection = left_present & right_present
        union = left_present | right_present
        yield {
            "sweep_id": sweep_id,
            "pair_id": pair_id,
            "seed": homogeneous.seed,
            "diet": homogeneous.diet,
            "step": step,
            "elapsed_seconds": elapsed,
            "left_condition": homogeneous.condition,
            "right_condition": spatial.condition,
            "left_present_species_count": len(left_present),
            "right_present_species_count": len(right_present),
            "intersection_size": len(intersection),
            "union_size": len(union),
            "jaccard_community_distance": jaccard_community_distance(left, right),
        }


def _extinction_rows(
    config: SweepConfig,
    pair_id: str,
    run: ConditionRun,
) -> Iterable[Mapping[str, object]]:
    censor_time = run.elapsed_seconds[-1]
    final_counts = run.population_counts[-1]
    for species, initial_count in run.initial_counts.items():
        if not 0 < initial_count <= config.low_abundance_max_initial_count:
            continue
        extinction_time = run.extinction_times[species]
        observed = extinction_time is not None
        yield {
            "sweep_id": config.sweep_id,
            "pair_id": pair_id,
            "seed": run.seed,
            "diet": run.diet,
            "condition": run.condition,
            "species": species,
            "initial_count": initial_count,
            "final_count": final_counts[species],
            "low_abundance_max_initial_count": config.low_abundance_max_initial_count,
            "low_abundance_definition": "0 < initial_count <= low_abundance_max_initial_count",
            "low_abundance_evidence_class": LOW_ABUNDANCE_EVIDENCE_CLASS,
            "event_observed": observed,
            "status": "extinct" if observed else "right_censored",
            "time_to_event_seconds": extinction_time if observed else censor_time,
            "censor_time_seconds": censor_time,
        }


def _metadata_row(
    config: SweepConfig,
    pair_id: str,
    run: ConditionRun,
    control: SpatialControl,
    species_data_sha256: str,
) -> Mapping[str, object]:
    return {
        "sweep_id": config.sweep_id,
        "pair_id": pair_id,
        "seed": run.seed,
        "diet": run.diet,
        "condition": run.condition,
        "dt_seconds": config.dt_seconds,
        "step_count": config.steps,
        "timepoint_count": config.steps + 1,
        "simulation_engine": "legacy_pg_production",
        "backend": config.backend,
        "physics_backend_used": run.physics_backend_used,
        "ecm_backend_used": run.ecm_backend_used,
        "metabolic_model_id": MODEL_ID,
        "metabolic_model_revision": MODEL_REVISION,
        "metabolic_prior_id": PRIOR_ID,
        "metabolic_prior_revision": PRIOR_REVISION,
        "reduced_model_active": False,
        "available_fba_lifecycle_policy": "GrowthOnlyFbaLifecyclePolicy",
        "fba_extinction_dynamics_validated": False,
        "fba_lifecycle_limitation": FBA_LIMITATION,
        "species_data_path": str(config.species_data_path.resolve()),
        "species_data_sha256": species_data_sha256,
        "founder_state_sha256": run.founder_state_sha256,
        "post_founder_rng_state_sha256": run.post_founder_rng_state_sha256,
        "world_width_um": config.world_width_um,
        "world_height_um": config.world_height_um,
        "carrying_capacity_cells": config.carrying_capacity_cells,
        "drop_frequency_s_inv": config.drop_frequency_s_inv,
        "drop_range_min_pg": config.drop_range_pg[0],
        "drop_range_max_pg": config.drop_range_pg[1],
        "low_abundance_max_initial_count": config.low_abundance_max_initial_count,
        "low_abundance_definition": "0 < initial_count <= low_abundance_max_initial_count",
        "low_abundance_evidence_class": LOW_ABUNDANCE_EVIDENCE_CLASS,
        "spatial_parameter_evidence_class": SPATIAL_EVIDENCE_CLASS,
        "spatial_partition_enabled": control.partition_enabled,
        "mucus_layer_radius_fraction": control.mucus_layer_radius_fraction,
        "mucus_motility_anchoring_factor": control.mucus_motility_factor,
        "lumen_food_washout_rate_s_inv": control.lumen_food_washout_rate_s_inv,
        "mucus_food_washout_rate_s_inv": control.mucus_food_washout_rate_s_inv,
        "mucus_drop_amount_factor": control.mucus_drop_amount_factor,
        "lumen_competition_cull_weight": control.lumen_competition_cull_weight,
        "initial_spawn_mode": INITIAL_SPAWN_MODE,
        "initial_species_arrangement": INITIAL_SPECIES_ARRANGEMENT,
        "initial_spawn_lumen_bias": INITIAL_SPAWN_LUMEN_BIAS,
        "initial_spawn_lumen_radius_fraction": INITIAL_SPAWN_LUMEN_RADIUS_FRACTION,
        "initial_spawn_cluster_count": INITIAL_SPAWN_CLUSTER_COUNT,
        "initial_spawn_cluster_spread_fraction": INITIAL_SPAWN_CLUSTER_SPREAD_FRACTION,
        "initial_spawn_central_radius_fraction": INITIAL_SPAWN_CENTRAL_RADIUS_FRACTION,
        "petri_radius_padding_um": PETRI_RADIUS_PADDING_UM,
        "ecm_cell_size_um": ECM_CELL_SIZE_UM,
    }


def run_comparative_sweep(config: SweepConfig) -> SweepResult:
    """Run each diet/seed pair for exactly ``config.steps`` updates per control."""
    species_data_sha256 = hashlib.sha256(config.species_data_path.read_bytes()).hexdigest()
    relative_rows = []
    jaccard_rows = []
    extinction_rows = []
    metadata_rows = []
    for diet in config.diets:
        for seed in config.seeds:
            pair_id = _pair_id(seed, diet)
            runs = {
                condition: _run_condition(config, seed, diet, SPATIAL_CONTROLS[condition])
                for condition in CONDITION_ORDER
            }
            homogeneous = runs[CONDITION_HOMOGENEOUS]
            spatial = runs[CONDITION_SPATIAL_NICHE]
            if homogeneous.founder_state_sha256 != spatial.founder_state_sha256:
                raise RuntimeError(f"paired founder states differ for {pair_id}")
            if (
                homogeneous.post_founder_rng_state_sha256
                != spatial.post_founder_rng_state_sha256
            ):
                raise RuntimeError(f"paired post-founder RNG states differ for {pair_id}")
            for condition in CONDITION_ORDER:
                run = runs[condition]
                relative_rows.extend(_relative_abundance_rows(config.sweep_id, pair_id, run))
                extinction_rows.extend(_extinction_rows(config, pair_id, run))
                metadata_rows.append(_metadata_row(
                    config,
                    pair_id,
                    run,
                    SPATIAL_CONTROLS[condition],
                    species_data_sha256,
                ))
            jaccard_rows.extend(_jaccard_rows(
                config.sweep_id,
                pair_id,
                homogeneous,
                spatial,
            ))
    return SweepResult(
        sweep_id=config.sweep_id,
        relative_abundance_rows=tuple(relative_rows),
        jaccard_rows=tuple(jaccard_rows),
        extinction_rows=tuple(extinction_rows),
        metadata_rows=tuple(metadata_rows),
    )


def allocate_output_directory(output_directory: Path, directory_name: str) -> Path:
    """Atomically allocate a fresh run directory without replacing prior outputs."""
    output_directory.mkdir(parents=True, exist_ok=True)
    safe_name = "".join(
        character if character.isalnum() or character in {"-", "_"} else "_"
        for character in directory_name
    ).strip("_")
    if not safe_name:
        raise ValueError("output directory name must contain a letter or number")
    for index in range(10_000):
        suffix = "" if index == 0 else f"-{index:03d}"
        candidate = output_directory / f"{safe_name}{suffix}"
        try:
            candidate.mkdir()
        except FileExistsError:
            continue
        return candidate
    raise FileExistsError("could not allocate a non-overwriting output directory")


def _write_csv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_sweep_outputs(
    result: SweepResult,
    output_directory: Path,
    directory_name: str | None = None,
) -> Mapping[str, Path]:
    """Write four schema-stable CSVs in a newly allocated run directory."""
    run_directory = allocate_output_directory(
        output_directory,
        directory_name or result.sweep_id,
    )
    rows_by_file = {
        "relative_abundance.csv": result.relative_abundance_rows,
        "jaccard_distance.csv": result.jaccard_rows,
        "extinction_censoring.csv": result.extinction_rows,
        "run_metadata.csv": result.metadata_rows,
    }
    paths = {}
    for filename, fields in CSV_SCHEMAS.items():
        path = run_directory / filename
        _write_csv(path, fields, rows_by_file[filename])
        paths[filename] = path
    return paths


def _default_sweep_id() -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return f"comparative_sweep_{timestamp}"


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-data", type=Path, default=DEFAULT_SPECIES_DATA_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIRECTORY)
    parser.add_argument("--seeds", nargs="+", type=int, default=[0])
    parser.add_argument("--diets", nargs="+", choices=tuple(DIET_PROFILES), default=list(DIET_PROFILES))
    parser.add_argument("--steps", type=int, default=1_000)
    parser.add_argument("--dt", type=float, default=0.1, dest="dt_seconds")
    parser.add_argument("--width", type=int, default=800, dest="world_width_um")
    parser.add_argument("--height", type=int, default=600, dest="world_height_um")
    parser.add_argument("--capacity", type=int, default=500, dest="carrying_capacity_cells")
    parser.add_argument("--drop-frequency", type=float, default=1.0, dest="drop_frequency_s_inv")
    parser.add_argument("--drop-min", type=float, default=1.0)
    parser.add_argument("--drop-max", type=float, default=50.0)
    parser.add_argument("--backend", choices=("auto", "python", "cpp"), default="python")
    parser.add_argument(
        "--low-abundance-max-initial-count",
        type=int,
        default=10,
    )
    parser.add_argument("--sweep-id", default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    config = SweepConfig(
        species_data_path=args.species_data,
        output_directory=args.output_dir,
        seeds=tuple(args.seeds),
        diets=tuple(args.diets),
        steps=args.steps,
        dt_seconds=args.dt_seconds,
        world_width_um=args.world_width_um,
        world_height_um=args.world_height_um,
        carrying_capacity_cells=args.carrying_capacity_cells,
        drop_frequency_s_inv=args.drop_frequency_s_inv,
        drop_range_pg=(args.drop_min, args.drop_max),
        backend=args.backend,
        low_abundance_max_initial_count=args.low_abundance_max_initial_count,
        sweep_id=args.sweep_id or _default_sweep_id(),
    )
    result = run_comparative_sweep(config)
    paths = write_sweep_outputs(result, config.output_directory)
    for filename, path in paths.items():
        print(f"{filename}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())