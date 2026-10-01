"""Inventory CSV configuration, loader defaults and fixed-temperature behavior.

This is a read-only research audit. A value present in CSV is not evidence that it
was measured. The tool neither imputes empirical parameters nor changes the loader.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import sys

os.environ.setdefault("PYGAME_HIDE_SUPPORT_PROMPT", "1")

from bacteria.constants import DEFAULT_ENV_TEMPERATURE_C
from bacteria.core_sim import BacteriaSpecies
from bacteria.data import SpeciesDataManager


_ALIASES = {"t6ss_activity": "aggression", "immunity_factor": "defense"}
_ATTRIBUTES = {
    "init": "init_count", "rep": "base_rep", "aging": "base_aging",
    "so": "o2_sens", "sf": "food_sens", "ph": "preferred_ph",
}


def audit_species_parameters(
    species_data_path: str | Path,
    temperature_c: float = DEFAULT_ENV_TEMPERATURE_C,
) -> dict:
    """Return a JSON-safe inventory without modifying inputs or running a simulation.

    Numeric validity concerns parsing/finite values only, not biological validity.
    Constructor clamping is reported separately rather than silently accepted as
    measured data. Duplicate CSV headers are rejected because DictReader loses them.
    """
    if not math.isfinite(temperature_c):
        raise ValueError("temperature_c must be finite")
    path = Path(species_data_path)
    source = path.read_bytes()
    reader = csv.DictReader(io.StringIO(source.decode("utf-8"), newline=""))
    headers = reader.fieldnames
    if not headers:
        raise ValueError("species CSV must have a header")
    if len(headers) != len(set(headers)):
        raise ValueError("duplicate CSV headers cannot be audited unambiguously")
    raw_rows = [(reader.line_num, row) for row in reader]
    manager = SpeciesDataManager(str(path))
    if source != path.read_bytes():
        raise ValueError("species CSV changed during the audit; retry on a stable input")
    if len(raw_rows) != len(manager.species):
        raise ValueError("raw and normalized row counts differ")
    fields = manager.fieldnames
    profiles = []
    issues = []
    origins = {"csv_explicit": 0, "legacy_alias": 0, "compatibility_default": 0}
    for index, ((line_end, raw), normalized) in enumerate(zip(raw_rows, manager.species), start=2):
        row_issues = []
        if None in raw:
            row_issues.append("row has more values than CSV headers")
        if not (raw.get("name") or "").strip():
            row_issues.append("species name is missing or blank")
        parameters = []
        numeric = {}
        for field in fields:
            alias = _ALIASES.get(field)
            if raw.get(field):
                origin, source_field = "csv_explicit", field
            elif alias and raw.get(alias):
                origin, source_field = "legacy_alias", alias
            else:
                origin, source_field = "compatibility_default", None
            origins[origin] += 1
            entry = {
                "field": field,
                "raw_value": raw.get(field),
                "source_field": source_field,
                "value_origin": origin,
                "loader_value": normalized[field],
                "empirical_provenance": "not_recorded_per_value",
            }
            if field != "name":
                try:
                    value = float(normalized[field])
                    if not math.isfinite(value):
                        raise ValueError("nonfinite value")
                    numeric[field] = value
                    entry["numeric_value"] = value
                except (ValueError, TypeError, OverflowError):
                    row_issues.append(f"{field}: loader value is not a finite number")
                    entry["numeric_value"] = None
            parameters.append(entry)
        profile = {
            "csv_row": index,
            "csv_line_end": line_end,
            "name": normalized["name"],
            "parameters": parameters,
            "issues": row_issues,
        }
        if not row_issues:
            # The production constructor converts initial counts to int and applies
            # bounds to several traits; observe those transformations directly.
            species = BacteriaSpecies(
                name=normalized["name"].strip(), color=(0, 0, 0), **numeric,
            )
            for entry in parameters:
                field = entry["field"]
                value = getattr(species, _ATTRIBUTES.get(field, field))
                entry["constructor_value"] = value
                before = entry.get("numeric_value", entry["loader_value"])
                entry["constructor_changed_value"] = value != before
            profile.update({
                "founder_count": species.init_count,
                "temperature_c": temperature_c,
                "temperature_factor": species.temperature_factor(temperature_c),
            })
        # The interactive app parses the loader's string directly with int(),
        # whereas the sweep/constructor path accepts float-convertible text.
        # Report this entrypoint difference without changing either loader.
        try:
            profile["interactive_app_init_count"] = max(0, int(normalized["init"]))
            profile["interactive_app_init_compatible"] = True
        except (ValueError, TypeError, OverflowError):
            profile["interactive_app_init_count"] = None
            profile["interactive_app_init_compatible"] = False
            row_issues.append("init: interactive app requires integer text and would skip this profile")
        issues.extend(f"record {index} (ending at CSV line {line_end}): {message}" for message in row_issues)
        profiles.append(profile)
    names = [profile["name"].strip() for profile in profiles]
    if len(names) != len(set(names)):
        issues.append("duplicate species names are ambiguous for community analysis")
    if not profiles:
        issues.append("CSV contains no species profiles")
    valid_profiles = [profile for profile in profiles if "founder_count" in profile]
    return {
        "schema_version": 1,
        "source_path": str(path.resolve()),
        "source_sha256": hashlib.sha256(source).hexdigest(),
        "interpretation": {
            "purpose": "configuration and implementation audit; not empirical validation",
            "provenance": "raw/default/alias origin is not scientific provenance",
            "empirical_provenance": "CSV has no per-value measured/inferred/assumed source ledger",
            "temperature": "constant temperature can retain a non-unit, parameter-dependent multiplier",
            "missing_values": "compatibility defaults are existing loader behavior, not an empirical imputation policy",
            "founder_counts": "constructor diagnostics, not executed populations; interactive app init parsing is checked separately",
            "csv_locations": "csv_row counts logical records including the header; csv_line_end is the physical ending line",
        },
        "csv_header_count": len(headers),
        "loader_field_count": len(fields),
        "numeric_field_count": len(fields) - 1,
        "trait_field_count_excluding_name_and_init": len(fields) - 2,
        "profile_count": len(profiles),
        "unknown_csv_headers": [field for field in headers if field not in fields and field not in _ALIASES.values()],
        "value_origin_counts": origins,
        "numeric_and_structural_issues": issues,
        "all_profiles_evaluated": len(valid_profiles) == len(profiles),
        "total_founders": sum(profile["founder_count"] for profile in valid_profiles) if len(valid_profiles) == len(profiles) else None,
        "inoculated_profile_count": sum(profile["founder_count"] > 0 for profile in valid_profiles) if len(valid_profiles) == len(profiles) else None,
        "profiles": profiles,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-data", type=Path, default=Path("species_data.csv"))
    parser.add_argument("--temperature", type=float, default=DEFAULT_ENV_TEMPERATURE_C)
    parser.add_argument("--output", type=Path, help="new JSON path; an existing file is never overwritten")
    parser.add_argument("--require-explicit", action="store_true", help="fail audit if loader compatibility defaults were used")
    args = parser.parse_args(argv)
    try:
        report = audit_species_parameters(args.species_data, args.temperature)
        output = json.dumps(report, indent=2, allow_nan=False) + "\n"
        if args.output:
            with args.output.open("x", encoding="utf-8") as destination:
                destination.write(output)
        else:
            sys.stdout.write(output)
    except (OSError, UnicodeError, ValueError) as error:
        print(f"parameter audit failed: {error}", file=sys.stderr)
        return 2
    if report["numeric_and_structural_issues"]:
        return 1
    if args.require_explicit and report["value_origin_counts"]["compatibility_default"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
