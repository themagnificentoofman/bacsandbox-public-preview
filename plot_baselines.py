#!/usr/bin/env python3
"""Phase I reproducible-baseline plotting (BacSandbox).

Reads an exported ``simulation_stats.csv`` (the multi-section CSV produced by
``Simulation.export``) and renders two baseline figures used for the malnutrition
extinction/abundance analysis deliverable:

  1. ``abundance_trajectories.png`` -- per-species population vs elapsed sim-time.
  2. ``extinction_survival.png``    -- a Kaplan-Meier-style species survival curve
     plus a bar chart of each species' time-to-extinction.
  3. ``zone_occupancy.png``         -- lumen vs mucus population over time (Phase I
     spatial heterogeneity), skipped when the CSV predates zone tracking.
  4. ``substrate_consumption.png``  -- cumulative per-substrate mass ingested by the
     community (Phase I diet analytics), skipped when those columns are absent.

The script is standalone (pandas + matplotlib only) and does not import the
simulation, so baselines can be regenerated from any archived CSV.

Usage:
    python plot_baselines.py --input simulation_stats.csv --outdir .
"""
from __future__ import annotations

import argparse
import os

import matplotlib

matplotlib.use("Agg")  # headless backend: write PNGs without a display server
import matplotlib.pyplot as plt
import pandas as pd


# Section markers written by Simulation.export(); parsing keys off these literals.
# GLOBAL_HEADER lists only the core (always-present) columns; it is used as a prefix to locate
# the Global Stats block. Newer exports append optional zone/diet-analytics columns that are
# discovered dynamically from the actual header row so legacy CSVs still parse.
GLOBAL_HEADER = ["time", "total", "births", "deaths", "elapsed_seconds",
                 "dt_seconds", "tickrate_hz", "overall_births", "overall_deaths"]
# Optional columns added in Phase I refinements (zone occupancy + per-substrate consumption).
OPTIONAL_GLOBAL_COLUMNS = ["lumen", "mucus",
                          "consumed_calorie", "consumed_glycan", "consumed_mdcf"]
PERSPECIES_HEADER = ["Index", "Species", "Step", "Count", "Births", "Deaths",
                     "Overall Births", "Overall Deaths"]
EXTINCTION_MARKER = "--- Extinction Events ---"
EXTINCTION_HEADER = ["Species", "Time_To_Extinction_s", "Extinct", "Final_Count"]
ENV_CONFIG_MARKER = "--- Environment Configuration ---"


def _read_rows(path: str) -> list[list[str]]:
    import csv
    with open(path, "r", newline="") as f:
        return [row for row in csv.reader(f)]


def _find_block(rows: list[list[str]], header: list[str]) -> list[list[str]]:
    """Return the contiguous data rows immediately following the row equal to ``header``."""
    out: list[list[str]] = []
    capture = False
    for row in rows:
        if not capture:
            if row[: len(header)] == header:
                capture = True
            continue
        if not row or not row[0] or row[0].startswith("---"):
            break
        out.append(row)
    return out


def _find_block_with_header(rows: list[list[str]], header_prefix: list[str]):
    """Locate a block whose header row starts with ``header_prefix``.

    Returns ``(matched_header_row, data_rows)``. The matched header row captures any optional
    trailing columns (e.g. lumen/mucus/consumed_*) so newer exports parse with their full schema
    while legacy CSVs (lacking those columns) still match on the core prefix.
    """
    out: list[list[str]] = []
    capture = False
    matched_header: list[str] | None = None
    for row in rows:
        if not capture:
            if row[: len(header_prefix)] == header_prefix:
                capture = True
                matched_header = row
            continue
        if not row or not row[0] or row[0].startswith("---"):
            break
        out.append(row)
    return matched_header, out


def _read_env_config(rows: list[list[str]]) -> dict[str, str]:
    """Parse the ``--- Environment Configuration ---`` key/value block into a dict (best-effort)."""
    cfg: dict[str, str] = {}
    capture = False
    for row in rows:
        if not capture:
            if row and row[0] == ENV_CONFIG_MARKER:
                capture = True
            continue
        if not row or not row[0] or row[0].startswith("---"):
            break
        if len(row) >= 2:
            cfg[row[0]] = row[1]
    return cfg


def parse_csv(path: str):
    """Parse the export into (global_df, species_df, extinction_df, env_config)."""
    rows = _read_rows(path)

    # Discover the actual Global Stats columns (core prefix + any optional trailing columns) so
    # newer exports with zone/consumption data parse fully while legacy CSVs still load.
    global_header, global_rows = _find_block_with_header(rows, GLOBAL_HEADER)
    columns = global_header if global_header else list(GLOBAL_HEADER)
    width = len(columns)
    # Normalise row widths (pad short rows, trim overlong ones) to match the header.
    norm_rows = [r[:width] + [""] * (width - len(r)) for r in global_rows]
    global_df = pd.DataFrame(norm_rows, columns=columns)
    for col in columns:
        global_df[col] = pd.to_numeric(global_df[col], errors="coerce")

    species_rows = _find_block(rows, PERSPECIES_HEADER)
    species_df = pd.DataFrame(species_rows, columns=PERSPECIES_HEADER)
    for col in ("Index", "Step", "Count", "Births", "Deaths",
                "Overall Births", "Overall Deaths"):
        species_df[col] = pd.to_numeric(species_df[col], errors="coerce")

    # Map per-species Step -> elapsed seconds via the global stats (same row order).
    step_to_elapsed = dict(zip(range(len(global_df)), global_df["elapsed_seconds"]))
    species_df["elapsed_seconds"] = species_df["Step"].map(step_to_elapsed)

    ext_rows = _find_block(rows, EXTINCTION_HEADER)
    extinction_df = pd.DataFrame(ext_rows, columns=EXTINCTION_HEADER)
    if not extinction_df.empty:
        extinction_df["Time_To_Extinction_s"] = pd.to_numeric(
            extinction_df["Time_To_Extinction_s"], errors="coerce")
        extinction_df["Extinct"] = pd.to_numeric(extinction_df["Extinct"], errors="coerce").fillna(0).astype(int)
        extinction_df["Final_Count"] = pd.to_numeric(extinction_df["Final_Count"], errors="coerce")

    env_config = _read_env_config(rows)
    return global_df, species_df, extinction_df, env_config


def plot_abundance(species_df: pd.DataFrame, outdir: str) -> str:
    fig, ax = plt.subplots(figsize=(10, 6))
    for name, grp in species_df.groupby("Species"):
        grp = grp.sort_values("elapsed_seconds")
        ax.plot(grp["elapsed_seconds"], grp["Count"], label=str(name), linewidth=1.5)
    ax.set_xlabel("Elapsed simulation time (s)")
    ax.set_ylabel("Population (cells)")
    ax.set_title("Species Abundance Trajectories")
    ax.legend(title="Species", loc="best")
    ax.grid(True, alpha=0.3)
    path = os.path.join(outdir, "abundance_trajectories.png")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_survival(extinction_df: pd.DataFrame, global_df: pd.DataFrame, outdir: str) -> str:
    fig, (ax_surv, ax_bar) = plt.subplots(1, 2, figsize=(13, 5))

    total_species = max(len(extinction_df), 1)
    # Kaplan-Meier-style step survival: fraction of species still extant over time.
    # Extinct species contribute a downward step at their time-to-extinction; survivors never drop.
    t_max = float(global_df["elapsed_seconds"].max()) if not global_df.empty else 0.0
    ext_times = sorted(
        t for t in extinction_df.loc[extinction_df["Extinct"] == 1, "Time_To_Extinction_s"].tolist()
        if pd.notna(t)
    )
    xs = [0.0]
    ys = [1.0]
    surviving = total_species
    for t in ext_times:
        xs.append(t); ys.append(surviving / total_species)   # horizontal up to the event
        surviving -= 1
        xs.append(t); ys.append(surviving / total_species)   # vertical drop at the event
    xs.append(max(t_max, xs[-1])); ys.append(surviving / total_species)
    ax_surv.step(xs, ys, where="post", color="#b03030", linewidth=2)
    ax_surv.set_ylim(-0.02, 1.05)
    ax_surv.set_xlabel("Elapsed simulation time (s)")
    ax_surv.set_ylabel("Fraction of species surviving")
    ax_surv.set_title("Species Survival Curve")
    ax_surv.grid(True, alpha=0.3)

    # Per-species time-to-extinction bar (survivors shown at the run horizon, hatched).
    names = extinction_df["Species"].tolist()
    times = []
    colors = []
    for _, r in extinction_df.iterrows():
        if int(r["Extinct"]) == 1 and pd.notna(r["Time_To_Extinction_s"]):
            times.append(float(r["Time_To_Extinction_s"]))
            colors.append("#b03030")
        else:
            times.append(t_max)            # survived to the end of the recorded run
            colors.append("#3070b0")
    ax_bar.bar(range(len(names)), times, color=colors)
    ax_bar.set_xticks(range(len(names)))
    ax_bar.set_xticklabels(names, rotation=45, ha="right")
    ax_bar.set_ylabel("Time-to-extinction (s)  /  run horizon")
    ax_bar.set_title("Time-to-Extinction by Species")
    ax_bar.grid(True, axis="y", alpha=0.3)

    path = os.path.join(outdir, "extinction_survival.png")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_zone_occupancy(global_df: pd.DataFrame, outdir: str) -> str | None:
    """Plot lumen vs mucus population over time. Returns the PNG path, or None if the CSV lacks
    zone columns (legacy export predating Phase I spatial heterogeneity tracking)."""
    if global_df.empty or "lumen" not in global_df.columns or "mucus" not in global_df.columns:
        return None
    if global_df["lumen"].notna().sum() == 0 and global_df["mucus"].notna().sum() == 0:
        return None
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(global_df["elapsed_seconds"], global_df["lumen"], label="Lumen", color="#3070b0", linewidth=1.5)
    ax.plot(global_df["elapsed_seconds"], global_df["mucus"], label="Mucus", color="#b07030", linewidth=1.5)
    ax.set_xlabel("Elapsed simulation time (s)")
    ax.set_ylabel("Population (cells)")
    ax.set_title("Zone Occupancy: Lumen vs Mucus")
    ax.legend(title="Spatial zone", loc="best")
    ax.grid(True, alpha=0.3)
    path = os.path.join(outdir, "zone_occupancy.png")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_substrate_consumption(global_df: pd.DataFrame, outdir: str,
                               diet_profile: str | None = None) -> str | None:
    """Plot cumulative per-substrate mass ingested by the community over time. Returns the PNG path,
    or None if the CSV lacks the consumption columns (legacy export predating diet analytics)."""
    cols = ["consumed_calorie", "consumed_glycan", "consumed_mdcf"]
    if global_df.empty or not all(c in global_df.columns for c in cols):
        return None
    if all(global_df[c].notna().sum() == 0 for c in cols):
        return None
    fig, ax = plt.subplots(figsize=(10, 6))
    # Substrate render colours mirror the in-sim SUBSTRATE_COLORS (calorie / glycan / MDCF).
    series = [
        ("consumed_calorie", "Generic calorie", "#c8c846"),
        ("consumed_glycan", "Complex glycan", "#d28c32"),
        ("consumed_mdcf", "MDCF prebiotic", "#8c5ac8"),
    ]
    for col, label, color in series:
        ax.plot(global_df["elapsed_seconds"], global_df[col], label=label, color=color, linewidth=1.5)
    ax.set_xlabel("Elapsed simulation time (s)")
    ax.set_ylabel("Cumulative mass ingested (pg)")
    title = "Cumulative Substrate Consumption"
    if diet_profile:
        title += f"  (diet: {diet_profile})"
    ax.set_title(title)
    ax.legend(title="Substrate class", loc="best")
    ax.grid(True, alpha=0.3)
    path = os.path.join(outdir, "substrate_consumption.png")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot BacSandbox Phase I baselines from an exported CSV.")
    ap.add_argument("--input", default="simulation_stats.csv", help="Path to exported simulation_stats.csv")
    ap.add_argument("--outdir", default=".", help="Directory to write the PNG figures into")
    args = ap.parse_args()

    if not os.path.exists(args.input):
        raise SystemExit(f"Input CSV not found: {args.input}")
    os.makedirs(args.outdir, exist_ok=True)

    global_df, species_df, extinction_df, env_config = parse_csv(args.input)

    if species_df.empty:
        raise SystemExit("No per-species history found in the CSV; nothing to plot.")

    p1 = plot_abundance(species_df, args.outdir)
    print(f"Wrote {p1}")
    if not extinction_df.empty:
        p2 = plot_survival(extinction_df, global_df, args.outdir)
        print(f"Wrote {p2}")
    else:
        print("No extinction-events section found; skipped survival curve.")

    p3 = plot_zone_occupancy(global_df, args.outdir)
    if p3:
        print(f"Wrote {p3}")
    else:
        print("No lumen/mucus columns found; skipped zone-occupancy plot.")

    p4 = plot_substrate_consumption(global_df, args.outdir, env_config.get("diet_profile"))
    if p4:
        print(f"Wrote {p4}")
    else:
        print("No per-substrate consumption columns found; skipped substrate plot.")


if __name__ == "__main__":
    main()
