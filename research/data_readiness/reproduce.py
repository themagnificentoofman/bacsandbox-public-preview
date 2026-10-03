"""Rebuild the shipped-configuration ledger without inferring empirical evidence.

Run from the repository root. No simulation, parameter change, or biological score.
"""
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def main():
    (HERE / "results").mkdir(exist_ok=True)
    species_path = ROOT / "species_data.csv"
    digest = hashlib.sha256(species_path.read_bytes()).hexdigest()
    with species_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        fields = reader.fieldnames
    with (HERE / "results/parameter_ledger.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["profile", "field", "configured_value", "software_origin", "evidence_class", "exact_strain_match", "empirical_source", "source_units", "condition", "source_csv_sha256", "source_location", "disposition"])
        for line, row in enumerate(rows, 2):
            for field in fields:
                writer.writerow([row["name"], field, row[field], "shipped_csv_explicit",
                                 "unknown", "unestablished", "", "", "", digest,
                                 f"species_data.csv:{line}:{field}",
                                 "configuration only; not an empirical measurement or adopted benchmark parameter"])
    print(f"Wrote {len(rows) * len(fields)} configuration-ledger rows; no empirical mapping inferred.")



if __name__ == "__main__":
    main()
