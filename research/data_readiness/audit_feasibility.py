"""Read-only ZHOU_BM001 feasibility diagnostics; no ecological simulation or fit.

Configuration probes use actual interfaces. Complete-baseline quota calculations
are resolution diagnostics, never runnable founders. See DOCUMENTATION.md 10.6.
"""
import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
import csv
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform

os.environ.setdefault("PYGAME_HIDE_SUPPORT_PROMPT", "1")

import pygame

from bacteria import constants
from bacteria.diet_translator import DietTranslationEngine
from bacteria.metabolic_integration import FbaAgentState, GrowthOnlyFbaLifecyclePolicy
from bacteria.metabolism import build_species_phenotype_config, concentration_dependent_uptake_bound
from research.data_readiness import map_benchmark
from run_comparative_sweeps import SweepConfig, HOMOGENEOUS_CONTROL, _environment

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BUDGETS = (500, 5000, 50000)  # Diagnostic examples, not selected model counts.
CODE_INPUTS = (
    "species_data.csv", "run_comparative_sweeps.py", "bacteria/core_sim.py",
    "bacteria/constants.py", "bacteria/data.py", "bacteria/metabolism.py",
    "bacteria/metabolic_integration.py", "bacteria/diet_translator.py",
    "bacteria/models/diet_profiles.json", "bacteria/models/reduced_metabolism_v1.json",
    "bacteria/models/metabolic_priors_v1.json", "research/data_readiness/map_benchmark.py",
    "research/data_readiness/extract_zhou.py", "research/data_readiness/provenance/supplements.json",
    "research/data_readiness/audit_feasibility.py",
)


def runtime_dependencies():
    """Name the imported Pygame provider, not another installed distribution.

    Both distributions expose ``pygame``. IS_CE identifies the loaded
    implementation; metadata availability alone cannot resolve co-installation.
    See DOCUMENTATION.md 10.6 for the audit provenance boundary.
    """
    provider = "pygame-ce" if getattr(pygame, "IS_CE", False) else "pygame"
    dependencies = {name: importlib.metadata.version(name)
                    for name in ("numpy", provider, "openpyxl")}
    installed = {provider: dependencies[provider]}
    other = "pygame" if provider == "pygame-ce" else "pygame-ce"
    try:
        installed[other] = importlib.metadata.version(other)
    except importlib.metadata.PackageNotFoundError:
        pass
    return dependencies, {
        "import_name": "pygame", "distribution": provider,
        "distribution_version": dependencies[provider],
        "module_version": pygame.version.ver,
        "module_matches_distribution_version": pygame.version.ver == dependencies[provider],
        "installed_distributions": installed,
        "coinstalled_distributions": len(installed) > 1,
    }


def founder_resolution(baselines, budgets):
    if not budgets or any(type(n) is not int or n <= 0 for n in budgets):
        raise ValueError("agent budgets must be positive integers")
    grouped = defaultdict(list)
    for row in baselines:
        grouped[row["animal_id"]].append(row)
    results, excluded = [], []
    for animal_id, rows in grouped.items():
        if len(rows) != 18 or len({r["strain_id"] for r in rows}) != 18:
            raise ValueError("resolution diagnostic requires 18 unique strain slots")
        unknown = [r for r in rows if r["measurement_state"] != "finite_nonzero_reported_log"]
        if unknown:
            excluded.append({"animal_id": animal_id, "unresolved_cells": [r["source_cell"] for r in unknown]})
            continue
        abundances = [r["derived_genome_units_per_g"] for r in rows]
        if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x <= 0
               for x in abundances):
            raise ValueError("invalid complete baseline abundance")
        total = math.fsum(abundances)
        if not math.isfinite(total):
            raise ValueError("nonfinite baseline total")
        fractions = [x / total for x in abundances]
        rare = min(range(len(rows)), key=fractions.__getitem__)
        target, = [i for i, r in enumerate(rows) if r["taxon"] == "S. copri"]
        for n in budgets:
            results.append({"animal_id": animal_id, **{k: rows[0][k] for k in map_benchmark.META},
                            "diagnostic_budget": n, "rarest_taxon": rows[rare]["taxon"],
                            "rarest_source_cell": rows[rare]["source_cell"],
                            "minimum_fraction": fractions[rare], "minimum_expected_quota": n * fractions[rare],
                            "taxa_with_expected_quota_below_one": sum(n * p < 1 for p in fractions),
                            "target_fraction": fractions[target], "target_expected_quota": n * fractions[target],
                            "budget_for_minimum_expected_quota_one": math.ceil(1 / fractions[rare]),
                            "forcing_one_rarest_agent_fraction_inflation": 1 / (n * fractions[rare]),
                            "evidence_class": "arithmetic_complete_baseline_subset_only",
                            "runnable_founders": False})
    return results, excluded


def diet_probe():
    translator = DietTranslationEngine()
    results = []
    for diet in map_benchmark.DIETS:
        for path in ("legacy_sweep", "molar_translator"):
            try:
                if path == "legacy_sweep":
                    SweepConfig(diets=(diet,))
                else:
                    # Positive placeholders test name acceptance, not study exposure.
                    translator.get_molar_substrate_bounds(diet, feed_mass_g=1.0, volume_L=1.0, dt_s=1.0)
            except ValueError as error:
                results.append({"diet": diet, "path": path, "accepted": False, "reason": str(error)})
            else:
                results.append({"diet": diet, "path": path, "accepted": True, "reason": "accepted by current interface"})
    return {"legacy_diets": list(constants.DIET_PROFILES), "molar_diets": list(translator.diet_names),
            "requested_diet_probes": results}


def profiles_audit(species_path):
    with species_path.open(newline="") as stream:
        profiles = list(csv.DictReader(stream))
    indexed = map_benchmark.index_unique(profiles, ["name"])
    reference = indexed[("a",)]
    fields = [f for f in reference if f not in ("name", "init")]
    config = build_species_phenotype_config(profiles)
    results = []
    for i, profile in enumerate(profiles):
        differences = {f: {"reference_a": reference[f], "profile": profile[f]}
                       for f in fields if profile[f] != reference[f]}
        results.append({"name": profile["name"], "configured_init": int(profile["init"]),
                        "non_name_non_init_fields": len(fields),
                        "fields_equal_to_synthetic_a": len(fields) - len(differences),
                        "differing_fields": json.dumps(differences, sort_keys=True),
                        "uptake_vmax_mol_gdw_s": json.dumps(config["uptake_vmax"][i].tolist()),
                        "uptake_km_molar": json.dumps(config["uptake_km"][i].tolist()),
                        "evidence_class": "current_configuration_not_empirical_measurement"})
    return results


def structural_probes():
    diet_asset = json.loads((ROOT / "bacteria/models/diet_profiles.json").read_text())
    routes = defaultdict(list)
    for component, properties in diet_asset["components"].items():
        routes[properties["substrate"]].append(component)
    # Own-model counterexample: simultaneous C/K scaling preserves the uptake bound.
    first = {"vmax": 1.8e-6, "concentration": 1e-4, "km": 1.5e-4}
    second = {"vmax": first["vmax"], "concentration": first["concentration"] * 10,
              "km": first["km"] * 10}
    counterexample = [{**inputs, "uptake_bound": concentration_dependent_uptake_bound(**inputs)}
                      for inputs in (first, second)]
    policy = GrowthOnlyFbaLifecyclePolicy(initial_biomass_gdw=1e-12)
    lifecycle = []
    for biomass in (0.5e-12, 1e-12, 2e-12):
        for capacity in (0.0, 1.0):
            result = policy.evaluate(None, FbaAgentState(biomass_gdw=biomass), dt_seconds=1.0,
                                     environmental_fitness=1.0, carrying_capacity_factor=capacity)
            lifecycle.append({"biomass_gdw": biomass, "capacity_factor": capacity, **asdict(result)})
    config = SweepConfig()
    return {"substrate_component_groups": dict(routes),
            "glucomannan_component_present": "glucomannan" in diet_asset["components"],
            "uptake_scale_counterexample": counterexample,
            "counterexample_scope": "equal local uptake bounds only; not whole-model invariance or empirical doses",
            "growth_only_policy_probes": lifecycle,
            "growth_only_scope": "policy has no starvation death; live simulator separately culls over-capacity agents",
            "default_harness": {"steps": config.steps, "dt_seconds": config.dt_seconds,
                                "horizon_seconds": config.steps * config.dt_seconds,
                                "capacity": config.carrying_capacity_cells,
                                "environment": _environment(config, config.diets[0]),
                                "homogeneous_control": asdict(HOMOGENEOUS_CONTROL)},
            "probe_only": True, "ecological_simulation_executed": False}


def time_workload(time_rows, dt_seconds):
    if not math.isfinite(dt_seconds) or dt_seconds <= 0:
        raise ValueError("dt_seconds must be positive and finite")
    output = []
    for row in time_rows:
        if row["phase_label"] == map_benchmark.PHASES[0]:
            continue
        for basis, field in (("literal_phase_label", "nominal_days_from_phase_label"),
                             ("reported_age_subtraction", "days_from_reported_age_subtraction")):
            days = row[field]
            output.append({"experiment": row["experiment"], "phase_label": row["phase_label"],
                           "interpretation": basis, "days": days, "dt_seconds": dt_seconds,
                           "literal_seconds": days * 86400,
                           "ticks_if_literal_seconds": math.ceil(days * 86400 / dt_seconds),
                           "clock_selected": False, "wall_time_measured": False})
    return output


def build(source_dir, budgets=DEFAULT_BUDGETS):
    before = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in CODE_INPUTS}
    mapping = map_benchmark.build(source_dir, ROOT / "species_data.csv")
    resolution, excluded = founder_resolution(mapping["baseline_observations.csv"], budgets)
    profiles = profiles_audit(ROOT / "species_data.csv")
    structural = structural_probes()
    targets = mapping["target_observations.csv"]
    baselines = {r["animal_id"]: r["raw_log10"] for r in targets if r["phase_label"] == map_benchmark.PHASES[0]}
    declines = []
    for phase in map_benchmark.PHASES[1:]:
        for diet in map_benchmark.DIETS:
            selected = [r for r in targets if r["phase_label"] == phase and r["diet"] == diet
                        and r["measurement_state"] == "finite_nonzero_reported_log"]
            declines.append({"phase": phase, "diet": diet, "n": len(selected),
                             "observed_below_own_baseline": sum(r["raw_log10"] < baselines[r["animal_id"]]
                                                               for r in selected)})
    after = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in CODE_INPUTS}
    if before != after:
        raise ValueError("code or asset changed during feasibility audit")
    dependencies, pygame_provenance = runtime_dependencies()
    audit = {"schema_version": 1, "benchmark_id": "ZHOU_BM001", "evidence_class": "configuration_and_arithmetic_feasibility",
             "code_sha256": before, "source_hashes": mapping["mapping_audit.json"]["source_hashes"],
             "python": platform.python_version(), "dependencies": dependencies,
             "pygame_provenance": pygame_provenance,
             "complete_baseline_animals": len({r["animal_id"] for r in resolution}),
             "complete_baseline_iterations": dict(Counter(r["experiment"] for r in mapping["animal_mapping.csv"]
                                                    if r["composition_values_complete"])),
             "excluded_incomplete_baselines": excluded,
             "observed_declines": declines, "runtime_performance_measured": False,
             "ecological_simulation_executed": False, "parameters_fitted": False,
             "benchmark_execution_ready": False,
             "note": "Excluded baselines are excluded only from quota diagnostics, not the benchmark cohort."}
    return {"feasibility_audit.json": audit, "diet_interface_probes.json": diet_probe(),
            "structural_probes.json": structural, "profile_provenance.csv": profiles,
            "founder_resolution.csv": resolution,
            "time_workload.csv": time_workload(mapping["time_crosswalk.csv"], SweepConfig().dt_seconds)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=map_benchmark.HERE / "raw/supplements")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        map_benchmark.write_outputs(build(args.source_dir), args.output)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
