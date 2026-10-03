"""Read-only CSV readiness checks; never fit, impute, threshold, or score biology.

Contracts describe an explicitly exported table, not a guessed workbook layout.
Run with --help. No simulator imports or third-party dependencies are required.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path


def read_table(path):
    """Keep literal cells; reject ambiguous headers and ragged CSV records."""
    path = Path(path)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, strict=True)
        header = next(reader, [])
        if not header or any(not h.strip() for h in header):
            raise ValueError("missing or blank header")
        if len(set(header)) != len(header):
            raise ValueError("duplicate header")
        rows = []
        for values in reader:
            if len(values) != len(header):
                raise ValueError(f"ragged CSV record ending at physical line {reader.line_num}")
            rows.append((reader.line_num, dict(zip(header, values))))
    return header, rows


def require_columns(header, columns):
    absent = sorted(set(columns) - set(header))
    if absent:
        raise ValueError(f"missing required columns: {absent}")


def key_index(rows, columns, missing_tokens):
    index = {}
    duplicates, missing = [], []
    for line, row in rows:
        key = tuple(row[c] for c in columns)
        if any(v.strip() in missing_tokens for v in key):
            missing.append(line)
        elif key in index:
            duplicates.append({"key": list(key), "first_line": index[key], "line": line})
        else:
            index[key] = line
    return index, duplicates, missing


def audit(path, contract, metadata=None):
    """Return structural evidence, including unresolved scientific metadata."""
    if not isinstance(contract, dict):
        raise ValueError("contract must be a JSON object")
    required = ("required_columns", "key_columns", "value_columns", "missing_tokens",
                "units", "time_mapping", "source_locator")
    if any(name not in contract for name in required):
        raise ValueError(f"contract must contain {required}")
    for field in ("required_columns", "key_columns", "value_columns", "missing_tokens"):
        value = contract[field]
        if not isinstance(value, list) or any(not isinstance(x, str) for x in value):
            raise ValueError(f"{field} must be a list of strings")
        if len(set(value)) != len(value):
            raise ValueError(f"duplicate entries in {field}")
    if not contract["key_columns"] or not contract["value_columns"]:
        raise ValueError("explicit key_columns and value_columns required")
    if set(contract["key_columns"]) & set(contract["value_columns"]):
        raise ValueError("key and value columns must be distinct")
    if not isinstance(contract["units"], dict):
        raise ValueError("units must map value columns to declared units or null")
    header, rows = read_table(path)
    require_columns(header, contract["required_columns"] + contract["key_columns"] + contract["value_columns"])
    missing_tokens = set(contract["missing_tokens"]) | {""}
    index, duplicates, missing_keys = key_index(rows, contract["key_columns"], missing_tokens)
    columns = {}
    invalid = []
    for col in contract["value_columns"]:
        counts = Counter()
        for line, row in rows:
            raw = row[col]
            if raw.strip() in missing_tokens:
                counts["missing"] += 1
                continue
            try:
                value = float(raw)
            except ValueError:
                counts["nonnumeric"] += 1
                invalid.append({"line": line, "column": col, "raw": raw, "reason": "nonnumeric"})
                continue
            if not math.isfinite(value):
                counts["nonfinite"] += 1
                invalid.append({"line": line, "column": col, "raw": raw, "reason": "nonfinite"})
                continue
            counts["finite"] += 1
            if value == 0:
                counts["zero"] += 1
            if value < 0:
                counts["negative"] += 1
        columns[col] = {k: counts[k] for k in ("missing", "finite", "zero", "negative", "nonnumeric", "nonfinite")}
    join = None
    if metadata is not None:
        join_columns = contract.get("join_columns")
        if not isinstance(join_columns, list) or not join_columns or any(not isinstance(x, str) for x in join_columns) or len(set(join_columns)) != len(join_columns):
            raise ValueError("metadata join requires explicit unique join_columns")
        mh, mr = read_table(metadata)
        require_columns(header, join_columns)
        require_columns(mh, join_columns)
        mi, md, mm = key_index(mr, join_columns, missing_tokens)
        di, _, dm = key_index(rows, join_columns, missing_tokens)
        join = {"columns": join_columns, "metadata_sha256": hashlib.sha256(Path(metadata).read_bytes()).hexdigest(),
                "duplicate_metadata_keys": md, "missing_metadata_keys": mm, "missing_data_keys": dm,
                "data_keys_without_metadata": [list(k) for k in sorted(di.keys() - mi.keys())],
                "metadata_keys_without_data": [list(k) for k in sorted(mi.keys() - di.keys())]}
    structural_ok = bool(rows) and not (duplicates or missing_keys or invalid)
    if join:
        structural_ok = structural_ok and not any(join[k] for k in ("duplicate_metadata_keys", "missing_metadata_keys", "missing_data_keys", "data_keys_without_metadata"))
    return {"schema_version": 1, "purpose": "data readiness only; no empirical evaluation",
            "source_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            "contract": contract, "row_count": len(rows), "header": header,
            "duplicate_keys": duplicates, "missing_key_lines": missing_keys,
            "columns": columns, "invalid_values": invalid, "metadata_join": join,
            "structural_ok": bool(structural_ok),
            "scientific_readiness": "not_assessed",
            "unknown_units": [c for c in contract["value_columns"] if not contract["units"].get(c)],
            "interpretation": "Zero, negative and missing values are counted, not classified as extinction. No blank correction, log transform, unit conversion, interpolation or row exclusion is performed."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output exists; choose a new path")
    try:
        result = audit(args.table, json.loads(args.contract.read_text()), args.metadata)
        # Validate the whole report before reserving a non-overwritable output.
        payload = json.dumps(result, indent=2, allow_nan=False) + "\n"
    except (OSError, ValueError, csv.Error) as exc:
        parser.error(str(exc))
    with args.output.open("x", encoding="utf-8") as stream:
        stream.write(payload)
    return 0 if result["structural_ok"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
