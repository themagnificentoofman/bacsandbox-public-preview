"""Extract pinned published tables without biological transformations or imputation.

Requires openpyxl==3.1.5. Run from any directory. Original workbooks are read-only.
Only metadata is propagated, within explicitly inspected table ranges. Every
propagated label retains its anchor cell; measurement blanks/errors stay literal.
"""
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import openpyxl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from audit_benchmark_table import audit

RAW = HERE / "raw/supplements"
MANIFEST = HERE / "provenance/supplements.json"
OUT = HERE / "results/extracted"
META = ["experiment", "community", "diet", "sex", "mouse_id"]


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_source(name, entries):
    data = (RAW / name).read_bytes()
    entry = entries[name]
    if len(data) != entry["size"] or hashlib.sha256(data).hexdigest() != entry["sha256"] or hashlib.md5(data).hexdigest() != entry["md5"]:
        raise ValueError(f"source integrity failure: {name}")
    return openpyxl.load_workbook(RAW / name, data_only=False)


def metadata(sheet, start, end, columns):
    """Section-bounded display-label propagation; no measurement forward fill."""
    state = {}
    for row in range(start, end + 1):
        result, anchors = {}, {}
        for column, name in columns.items():
            cell = sheet.cell(row, column)
            if cell.value is not None:
                state[name] = (cell.value, cell.coordinate)
            # Identifiers are never inherited, unlike hierarchical display labels.
            if name in ("mouse_id", "replicate") and cell.value is None:
                raise ValueError(f"missing identifier at {cell.coordinate}")
            if name not in state:
                raise ValueError(f"metadata has no anchor at {cell.coordinate}")
            result[name], anchors[name] = state[name]
        yield row, result, anchors


def numeric_value(cell):
    if cell.data_type == "f":
        raise ValueError(f"formula requires reviewed handling: {cell.coordinate}")
    # Excel errors remain strings; no dropping, numeric coercion or imputation.
    return cell.value


def summarize(rows, field):
    counts = Counter()
    issues = []
    for row in rows:
        value = row[field]
        if value is None:
            counts["blank"] += 1
            issues.append({"cell": row["value_cell"], "kind": "blank"})
        elif isinstance(value, (int, float)) and math.isfinite(value):
            counts["finite"] += 1
            counts["zero"] += value == 0
            counts["negative"] += value < 0
        else:
            counts["error_or_nonnumeric"] += 1
            issues.append({"cell": row["value_cell"], "kind": "error_or_nonnumeric", "raw": value})
    return {"records": len(rows), **{k: counts[k] for k in ("finite", "blank", "zero", "negative", "error_or_nonnumeric")}, "issues": issues}


def key(row, columns):
    return tuple(str(row[c]) for c in columns)


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    entries = {Path(e["key"]).name: e for e in json.loads(MANIFEST.read_text())["files"]}
    books = {n: load_source(n, entries) for n in entries if n.endswith(".xlsx")}
    # Also verify the SI PDF, without rewriting it.
    for n, entry in entries.items():
        if hashlib.sha256((RAW / n).read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"source integrity failure: {n}")
    s1 = books["pnas.2521522122.sd01.xlsx"]["Dataset S1"]
    s7 = books["pnas.2521522122.sd07.xlsx"]["Dataset S7"]
    s8 = books["pnas.2521522122.sd08.xlsx"]["Dataset S8"]
    s9 = books["pnas.2521522122.sd09.xlsx"]["Dataset S9"]
    if s1["E53"].value != "Optical Density (600 nm, normalized to starting OD600)":
        raise ValueError("S1C measurement header changed")
    if [s1.cell(54, c).value for c in range(5, 9)] != ["12 hr", "24 hrs", "36 hrs", "48 hrs"]:
        raise ValueError("S1C time headers changed")
    datasets = {}
    growth = []
    for r, meta, anchors in metadata(s1, 55, 66, {2: "strain", 3: "treatment", 4: "replicate"}):
        for c, hours in zip(range(5, 9), (12, 24, 36, 48)):
            cell = s1.cell(r, c)
            growth.append({**meta, "hours": hours, "value": numeric_value(cell), "units": s1["E53"].value,
                           "value_cell": cell.coordinate, "source_sheet": s1.title, "source_file": "pnas.2521522122.sd01.xlsx",
                           "metadata_anchors": json.dumps(anchors, sort_keys=True)})
    datasets["s1c_growth"] = growth
    roster = []
    for r, meta, anchors in metadata(s7, 4, 21, {2: "gavage_group", 3: "taxon", 4: "strain", 5: "reference"}):
        roster.append({**meta, "source_range": f"B{r}:E{r}", "metadata_anchors": json.dumps(anchors, sort_keys=True)})
    write_csv(OUT / "s7_roster.csv", roster)
    mice = []
    columns = dict(zip(range(2, 7), META))
    for r, meta, anchors in metadata(s8, 4, 80, columns):
        mice.append({**meta, "weight_g": numeric_value(s8.cell(r, 7)), "normalized_weight": numeric_value(s8.cell(r, 8)),
                     "source_range": f"B{r}:H{r}", "metadata_anchors": json.dumps(anchors, sort_keys=True)})
    write_csv(OUT / "s8_mouse_weights.csv", mice)
    a = []
    for r, meta, anchors in metadata(s9, 7, 45, columns):
        for c, logc in zip(range(7, 10), range(11, 14)):
            cell = s9.cell(r, c)
            a.append({**meta, "phase": s9.cell(6, c).value, "specimen": "feces", "organism": "S. copri",
                      "value": numeric_value(cell), "reported_log10": numeric_value(s9.cell(r, logc)),
                      "units": "genome units / g feces", "value_cell": cell.coordinate,
                      "log10_cell": s9.cell(r, logc).coordinate, "source_sheet": s9.title,
                      "source_file": "pnas.2521522122.sd09.xlsx", "metadata_anchors": json.dumps(anchors, sort_keys=True)})
    datasets["s9a_copri_fecal"] = a
    for label, start, end, header, last_col, specimen, unit_cell in [
        ("b", 83, 199, 82, 24, "feces", "G81"),
        ("c", 226, 264, 225, 25, "cecal contents", "G224"),
        ("d", 290, 403, 289, 24, "feces", "G288"),
        ("e", 429, 466, 428, 24, "cecal contents", "G427")]:
        if s9.cell(header, 7).value != "Timepoint":
            raise ValueError(f"S9{label} schema changed")
        rows = []
        for r, meta, anchors in metadata(s9, start, end, columns | {7: "phase"}):
            for c in range(8, last_col + 1):
                cell = s9.cell(r, c)
                rows.append({**meta, "specimen": specimen, "organism": s9.cell(header, c).value,
                             "value": numeric_value(cell), "units": s9[unit_cell].value,
                             "value_cell": cell.coordinate, "source_sheet": s9.title,
                             "source_file": "pnas.2521522122.sd09.xlsx",
                             "metadata_anchors": json.dumps(anchors, sort_keys=True)})
        datasets[f"s9{label}_community"] = rows
    report = {"source_hashes_verified": {n: e["sha256"] for n, e in entries.items()},
              "workbooks": {n: [{"sheet": s.title, "rows": s.max_row, "columns": s.max_column,
                                  "merged_ranges": len(s.merged_cells.ranges), "formula_cells": sum(c.data_type == "f" for row in s for c in row)} for s in b] for n, b in books.items()},
              "roster_count": len(roster), "gavage_groups": dict(Counter(r["gavage_group"] for r in roster)),
              "mouse_rows": len(mice), "datasets": {}, "scientific_readiness": "mentor_and_observation_decisions_pending"}
    for name, rows in datasets.items():
        write_csv(OUT / (name + ".csv"), rows)
        keys = ["strain", "treatment", "replicate", "hours"] if name == "s1c_growth" else META + ["phase", "specimen", "organism"]
        contract = {"required_columns": keys + ["value"], "key_columns": keys, "value_columns": ["value"],
                    "missing_tokens": [], "units": {"value": rows[0]["units"]},
                    "time_mapping": "literal hours" if name == "s1c_growth" else "phase labels preserved; postnatal-day discrepancies unresolved",
                    "source_locator": {"workbook": rows[0]["source_file"], "sheet": rows[0]["source_sheet"], "cells": "value_cell in every record"}}
        if name != "s1c_growth":
            contract["join_columns"] = META
        (OUT / (name + ".contract.json")).write_text(json.dumps(contract, indent=2) + "\n")
        result = audit(OUT / (name + ".csv"), contract, None if name == "s1c_growth" else OUT / "s8_mouse_weights.csv")
        (OUT / (name + ".audit.json")).write_text(json.dumps(result, indent=2) + "\n")
        report["datasets"][name] = summarize(rows, "value") | {"structural_ok": result["structural_ok"], "join_data_keys_without_metadata": [] if name == "s1c_growth" else result["metadata_join"]["data_keys_without_metadata"]}
    # Identity requires experimental arm as well as mouse label; diagnose collisions.
    report["mouse_id_collisions"] = {}
    for fields in (["sex", "mouse_id"], ["experiment", "sex", "mouse_id"], META):
        counts = Counter(key(r, fields) for r in mice)
        report["mouse_id_collisions"]["+".join(fields)] = {"distinct_keys": len(counts), "duplicate_key_groups": sum(n > 1 for n in counts.values())}
    report["animal_groups"] = [dict(zip(META[:3], k), count=n) for k, n in sorted(Counter(key(r, META[:3]) for r in mice).items())]
    report["roster_to_observed_taxon_alias"] = {"S. copri": "Segatella copri"}
    observed = {r["organism"] for name, rows in datasets.items() if name.startswith("s9") for r in rows}
    normalized = {report["roster_to_observed_taxon_alias"].get(n, n) for n in observed}
    report["unmatched_observed_taxa"] = sorted(normalized - {r["taxon"] for r in roster})
    report["unmatched_roster_taxa"] = sorted({r["taxon"] for r in roster} - normalized)
    # Compare provided scales without altering either; tolerance is not an outcome threshold.
    report["s9a_log10_max_abs_difference"] = max(abs(math.log10(r["value"]) - r["reported_log10"]) for r in a if r["value"] is not None)
    report["s9a_paired_missing_scales"] = sum(r["value"] is None and r["reported_log10"] is None for r in a)
    samples = {}
    for name, rows in datasets.items():
        if not name.startswith("s9"):
            continue
        for row in rows:
            k = key(row, META + ["phase", "specimen"])
            samples.setdefault(k, []).append(row)
    report["sample_slots"] = len(samples)
    report["sample_slots_with_any_value"] = sum(any(r["value"] is not None for r in rs) for rs in samples.values())
    report["entirely_blank_samples"] = [dict(zip(META + ["phase", "specimen"], k), cells=[r["value_cell"] for r in rows]) for k, rows in samples.items() if all(r["value"] is None for r in rows)]
    report["extracted_csv_hashes"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.glob("*.csv"))}
    (OUT / "inventory.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    report = build()
    print(json.dumps({"roster_count": report["roster_count"], "mouse_rows": report["mouse_rows"],
                      "datasets": {n: {k: v for k, v in s.items() if k != "issues"} for n, s in report["datasets"].items()},
                      "sample_slots": report["sample_slots"], "sample_slots_with_any_value": report["sample_slots_with_any_value"],
                      "entirely_blank_samples": report["entirely_blank_samples"], "mouse_id_collisions": report["mouse_id_collisions"]}, indent=2))
