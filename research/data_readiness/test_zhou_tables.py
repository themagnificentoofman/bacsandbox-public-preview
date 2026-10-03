"""Readiness checks against the actual checksum-pinned published workbooks."""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from research.data_readiness import extract_zhou


class PublishedTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        original = extract_zhou.OUT
        cls.addClassCleanup(setattr, extract_zhou, "OUT", original)
        extract_zhou.OUT = Path(cls.temp.name)
        cls.before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in extract_zhou.RAW.iterdir() if p.is_file()}
        cls.report = extract_zhou.build()

    def rows(self, name):
        with (extract_zhou.OUT / (name + ".csv")).open(newline="") as stream:
            return list(csv.DictReader(stream))

    def test_all_extracted_measurements_match_original_cells(self):
        for name in self.report["datasets"]:
            rows = self.rows(name)
            book = extract_zhou.openpyxl.load_workbook(extract_zhou.RAW / rows[0]["source_file"], data_only=False)
            for row in rows:
                value = book[row["source_sheet"]][row["value_cell"]].value
                self.assertEqual(row["value"], "" if value is None else str(value), row["value_cell"])

    def test_growth_grid_and_reported_negative_and_zero_values(self):
        rows = self.rows("s1c_growth")
        counts = Counter((r["strain"], r["treatment"], r["replicate"]) for r in rows)
        self.assertEqual(len(counts), 12)
        self.assertEqual(set(counts.values()), {4})
        self.assertEqual({r["hours"] for r in rows}, {"12", "24", "36", "48"})
        self.assertEqual(sum(float(r["value"]) < 0 for r in rows), 4)
        self.assertEqual(sum(float(r["value"]) == 0 for r in rows), 4)
        self.assertEqual(self.report["datasets"]["s1c_growth"]["blank"], 0)

    def test_excel_errors_fail_structure_without_becoming_missing(self):
        rows = self.rows("s9b_community")
        errors = [r for r in rows if r["value"] == "#NUM!"]
        self.assertEqual([r["value_cell"] for r in errors], ["V131", "V132", "V133", "V134", "V135", "V136", "V138", "V139"])
        self.assertEqual({r["organism"] for r in errors}, {"Ruminococcus torques"})
        self.assertFalse(self.report["datasets"]["s9b_community"]["structural_ok"])
        self.assertEqual(sum(r["value"] == "" for r in rows), 17)

    def test_mouse_joins_need_arm_context(self):
        self.assertEqual(self.report["mouse_rows"], 77)
        for report in self.report["datasets"].values():
            self.assertEqual(report["join_data_keys_without_metadata"], [])
        self.assertGreater(self.report["mouse_id_collisions"]["experiment+sex+mouse_id"]["duplicate_key_groups"], 0)
        self.assertEqual(self.report["mouse_id_collisions"]["experiment+community+diet+sex+mouse_id"]["duplicate_key_groups"], 0)

    def test_missing_samples_are_retained_with_correct_metadata(self):
        blanks = self.report["entirely_blank_samples"]
        self.assertEqual(len(blanks), 2)
        self.assertEqual(self.report["sample_slots"], 308)
        self.assertEqual(self.report["sample_slots_with_any_value"], 306)
        self.assertEqual((blanks[0]["sex"], blanks[0]["mouse_id"], blanks[0]["diet"]), ("Male", "14", "Mirpur-18 + Glucomannan"))
        self.assertEqual((blanks[1]["sex"], blanks[1]["mouse_id"], blanks[1]["diet"]), ("Female", "16", "Mirpur-18"))
        self.assertIn("I35", blanks[0]["cells"])
        self.assertIn("H398", blanks[1]["cells"])

    def test_roster_taxon_join_and_gavage_groups(self):
        self.assertEqual(self.report["gavage_groups"], {"Arm Dependent": 1, "Early": 5, "Late": 12})
        self.assertEqual(self.report["unmatched_roster_taxa"], [])
        self.assertEqual(self.report["unmatched_observed_taxa"], [])
        self.assertEqual(len(self.rows("s7_roster")), 18)

    def test_two_supplied_abundance_scales_agree_without_transforming_exports(self):
        self.assertLess(self.report["s9a_log10_max_abs_difference"], 1e-12)
        self.assertEqual(self.report["s9a_paired_missing_scales"], 1)
        self.assertEqual(len(self.rows("s9a_copri_fecal")), 117)

    def test_extraction_preserves_original_bytes(self):
        after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in extract_zhou.RAW.iterdir() if p.is_file()}
        self.assertEqual(after, self.before)


if __name__ == "__main__":
    unittest.main()
