"""Plot published observations only; no predictions, fitting, or success labels.

Optional plotting dependency: matplotlib==3.10.8. Run preparation first.
Horizontal bars are available-observation means, not uncertainty intervals.
"""
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from research.data_readiness.prepare_observations import OUT, PHASES, read_csv


def main():
    rows = read_csv(OUT / "observations.csv")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    colors = ["#355C8C", "#BD663D"]
    counts = []
    for axis, section, item in zip(axes.flat,
            ["s1c_growth", "s1c_growth", "s9a_copri_fecal", "s9a_copri_fecal"],
            ["DSM 18205", "BgF5_2", "Iteration 1", "Iteration 2"]):
        culture = section == "s1c_growth"
        subset = [r for r in rows if r["source_section"] == section and r["taxon_or_strain" if culture else "experiment"] == item]
        conditions = ["1% glucose", "1.5% glucomannan"] if culture else ["Mirpur-18", "Mirpur-18 + Glucomannan"]
        phases = ["12", "24", "36", "48"] if culture else PHASES
        timefield = "time_hours" if culture else "phase_label"
        for i, condition in enumerate(conditions):
            for x, phase in enumerate(phases):
                records = [r for r in subset if r["condition"] == condition and r[timefield] == phase]
                values = [float(r["numeric_value"]) for r in records if r["measurement_state"] == "observed_finite"]
                center = x + (i - 0.5) * 0.30
                positions = [center + (j - (len(values) - 1) / 2) * 0.018 for j in range(len(values))]
                axis.scatter(positions, values, s=20, alpha=0.75, color=colors[i], label=condition if x == 0 else None)
                axis.plot([center - 0.09, center + 0.09], [statistics.mean(values)] * 2, color=colors[i], lw=2)
                counts.append({"section": section, "panel": item, "condition": condition, "phase": phase,
                               "n_finite": len(values), "n_slots": len(records)})
        axis.set_title(("S1C · " if culture else "S9A · Copri · ") + item, loc="left", fontweight="bold")
        axis.set_xticks(range(len(phases)), phases if culture else ["Weaning", "+10 days", "+4 weeks"])
        axis.set_xlabel("Hours (source labels)" if culture else "Published phase labels; actual ages unresolved")
        axis.set_ylabel("Published normalized OD600" if culture else "Published log10(genome units / g feces)")
        # One shared range per assay, derived from every finite value; no clipping.
        assay_values = [float(r["numeric_value"]) for r in rows if r["source_section"] == section and r["measurement_state"] == "observed_finite"]
        low, high = min(assay_values), max(assay_values)
        padding = (high - low) * 0.12 or 0.1
        axis.set_ylim(low - padding, high + padding)
        axis.grid(axis="y", alpha=0.2)
        axis.legend(fontsize=8, loc="upper left" if culture else "lower right", frameon=False)
        if item == "Iteration 2":
            axis.text(0.02, 0.97, "+4 weeks GM: 9/10 available; one blank retained in table", transform=axis.transAxes, fontsize=8, va="top")
    fig.suptitle("Zhou published measurements · descriptive preparation", x=0.075, ha="left", fontsize=17, fontweight="bold")
    fig.text(0.075, 0.018, "Points: published observations (small deterministic horizontal offsets). Bars: available-observation means.\n"
             "S1C: 3 labels per group; normalization and replicate independence unresolved; substrate concentrations differ.\n"
             "S9A: iterations kept separate. No model predictions, uncertainty estimates, endpoint choice, or success criterion.", fontsize=9)
    fig.tight_layout(rect=(0.04, 0.095, 1, 0.95), h_pad=2)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / ("published_observations." + suffix), dpi=170,
                    metadata={"Creator": "BacSandbox source-readiness diagnostics", "CreationDate": None, "ModDate": None} if suffix == "pdf" else None)
    plt.close(fig)
    (OUT / "plot_counts.json").write_text(json.dumps({"matplotlib": matplotlib.__version__, "counts": counts}, indent=2) + "\n")


if __name__ == "__main__":
    main()
