"""Focused mapping checks: ambiguous states, arm identity and pinned source cells."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from research.data_readiness import map_benchmark as mapping


class MappingContractTests(unittest.TestCase):
    def test_distinct_measurement_states_and_inverse_log(self):
        cases = [(None, "missing_source_cell"), ("#NUM!", "source_excel_error"),
                 (0, "log_zero_semantics_unresolved"), (False, "nonnumeric"),
                 ("0", "nonnumeric"), (float("nan"), "nonfinite"),
                 (float("inf"), "nonfinite"), (400, "inverse_log_out_of_range"),
                 (-400, "inverse_log_out_of_range")]
        for value, state in cases:
            with self.subTest(value=value):
                self.assertEqual(mapping.classify_log(value), (state, None))
        self.assertEqual(mapping.classify_log(-1), ("finite_nonzero_reported_log", 0.1))
        self.assertEqual(mapping.classify_log(8), ("finite_nonzero_reported_log", 1e8))

    def test_compound_identity_preserves_arm_and_literal_identifier(self):
        rows = [{"diet": d, "mouse_id": m} for d in mapping.DIETS for m in ("2", "2L", "02")]
        self.assertEqual(len(mapping.index_unique(rows, ["diet", "mouse_id"])), 6)
        with self.assertRaises(ValueError):
            mapping.index_unique(rows, ["mouse_id"])
        for value in (None, "", " "):
            with self.assertRaises(ValueError):
                mapping.index_unique([{"id": value}], ["id"])

    def test_output_refuses_overwrite_and_invalid_json_writes_nothing(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "mapping"
            with self.assertRaises(ValueError):
                mapping.write_outputs({"bad.json": {"value": float("nan")}}, output)
            self.assertFalse(output.exists())
            mapping.write_outputs({"audit.json": {"ok": True}}, output)
            before = (output / "audit.json").read_bytes()
            hashes = json.loads((output / "checksums.json").read_text())
            self.assertEqual(hashes["audit.json"], hashlib.sha256(before).hexdigest())
            with self.assertRaises(ValueError):
                mapping.write_outputs({"audit.json": {"ok": False}}, output)
            self.assertEqual((output / "audit.json").read_bytes(), before)


class PublishedMappingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = mapping.HERE / "raw/supplements"
        cls.before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in cls.raw.iterdir() if p.is_file()}
        cls.outputs = mapping.build(cls.raw, mapping.ROOT / "species_data.csv")

    def test_all_702_baseline_slots_retain_source_values_and_zeros(self):
        sheets, _ = mapping.load_verified(self.raw)
        rows = self.outputs["baseline_observations.csv"]
        self.assertEqual(len(rows), 702)
        for row in rows:
            self.assertEqual(row["raw_log10"], sheets[9][row["source_cell"]].value)
            if row["raw_log10"] == 0:
                self.assertIsNone(row["derived_genome_units_per_g"])
        self.assertEqual(Counter(r["measurement_state"] for r in rows),
                         {"finite_nonzero_reported_log": 653, "log_zero_semantics_unresolved": 49})
        audit = self.outputs["mapping_audit.json"]
        self.assertEqual(audit["animals_with_unresolved_baseline_values"], 31)
        self.assertFalse(audit["simulation_ready"])
        self.assertTrue(all(not r["agent_counts_available"] for r in self.outputs["animal_mapping.csv"]))

    def test_missing_endpoint_does_not_delete_earlier_or_control_arm(self):
        rows = self.outputs["endpoint_eligibility.csv"]
        counts = {p: sum(r["endpoint_eligible"] for r in rows if r["phase_label"] == p)
                  for p in mapping.PHASES}
        self.assertEqual(counts, dict(zip(mapping.PHASES, (0, 39, 38))))
        treated = [r for r in rows if mapping.key(r) == mapping.MISSING_KEY]
        self.assertEqual([r["measurement_available"] for r in treated], [True, True, False])
        control_key = (*mapping.MISSING_KEY[:2], mapping.DIETS[1], *mapping.MISSING_KEY[3:])
        control = [r for r in rows if mapping.key(r) == control_key]
        self.assertEqual(len(control), 3)
        self.assertTrue(all(r["measurement_available"] for r in control))
        self.assertTrue(control[-1]["paired_change_eligible"])

    def test_missing_endpoint_corroborates_other_tables_without_substitution(self):
        evidence, = self.outputs["missing_endpoint.json"]
        self.assertEqual(evidence["fecal_target_log_cell"], "M35")
        self.assertEqual(evidence["fecal_target_linear_cell"], "I35")
        self.assertEqual(evidence["companion_cells"], [f"{chr(c)}189" for c in range(72, 89)])
        self.assertTrue(evidence["all_companions_missing"])
        self.assertEqual(evidence["cecal_evidence"]["source_cell"], "H254")
        self.assertAlmostEqual(evidence["cecal_evidence"]["raw_log10"], 10.2051111301871)
        self.assertEqual(evidence["s8_weight_evidence"]["weight_g"], 19.3)
        target, = [r for r in self.outputs["target_observations.csv"] if r["source_cell"] == "M35"]
        self.assertIsNone(target["raw_log10"])
        self.assertIsNone(target["derived_genome_units_per_g"])

    def test_roster_and_time_mapping_do_not_claim_parameters_or_seconds(self):
        roster = self.outputs["strain_crosswalk.csv"]
        self.assertEqual(len({r["strain_id"] for r in roster}), 18)
        candidates = [r for r in roster if r["model_profile_candidate"]]
        self.assertEqual([r["strain"] for r in candidates], ["BgF5_2"])
        self.assertEqual(candidates[0]["mapping_status"], "label_candidate_parameters_unvalidated")
        clock = self.outputs["time_crosswalk.csv"]
        self.assertEqual(len(clock), 6)
        self.assertTrue(all(r["model_elapsed_seconds"] is None for r in clock))
        self.assertEqual([r["days_from_reported_age_subtraction"] for r in clock], [0, 11, 29, 0, 12, 30])

    def test_sources_unchanged_and_mismatched_bytes_rejected(self):
        after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in self.raw.iterdir() if p.is_file()}
        self.assertEqual(after, self.before)
        with tempfile.TemporaryDirectory() as temp:
            name = "pnas.2521522122.sd07.xlsx"
            (Path(temp) / name).write_bytes((self.raw / name).read_bytes() + b"changed")
            with self.assertRaisesRegex(ValueError, "integrity failure"):
                mapping.load_verified(Path(temp))


if __name__ == "__main__":
    unittest.main()
