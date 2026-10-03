"""Source retention and diagnostic checks; synthetic cases only test parser behavior."""
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import openpyxl

from research.data_readiness import extract_zhou as extract
from research.data_readiness import prepare_observations as prep


class PreparationEdgeTests(unittest.TestCase):
    def test_missing_error_nonfinite_and_negative_values_stay_distinct(self):
        cases = {"": (None, "missing_source_cell"), "#NUM!": (None, "source_excel_error"),
                 "NaN": (None, "nonfinite"), "inf": (None, "nonfinite"),
                 "unknown": (None, "nonnumeric"), "-0.02": (-0.02, "observed_finite"),
                 "0": (0.0, "observed_finite")}
        for raw, expected in cases.items():
            self.assertEqual(prep.parse_value(raw), expected)

    def test_identifiers_with_delimiters_do_not_collide(self):
        self.assertNotEqual(prep.stable_id(["a|b", "c"])[0], prep.stable_id(["a", "b|c"])[0])
        self.assertEqual(json.loads(prep.stable_id(["5L", "10.1"])[1]), ["5L", "10.1"])

    def test_flagged_group_is_not_silently_filtered_for_formula_check(self):
        for value in ("", "#NUM!", "0", "NaN", "unknown"):
            result = prep.check_summary([{"value": "9", "value_cell": "A1"}, {"value": value, "value_cell": "A2"}], "9 ± 0")
            self.assertEqual(result["diagnostic_status"], "not_tested_flagged_group")
            self.assertIsNone(result["mean_log_sd_log_matches"])
            self.assertIsNone(result["log_linear_mean_upper_sd_matches"])
            self.assertEqual(result["n_slots"], 2)

    def test_sample_sd_and_backtransformation_are_distinguished(self):
        rows = [{"value": "1", "value_cell": "A1"}, {"value": "3", "value_cell": "A2"}]
        result = prep.check_summary(rows, "2.00 ± 1.41")
        self.assertTrue(result["mean_log_sd_log_matches"])
        self.assertFalse(result["log_linear_mean_upper_sd_matches"])
        result = prep.check_summary(rows, "2.70 ± 0.38")
        self.assertFalse(result["mean_log_sd_log_matches"])
        self.assertTrue(result["log_linear_mean_upper_sd_matches"])

    def test_extreme_log_values_do_not_overflow(self):
        rows = [{"value": "1000", "value_cell": "A1"}, {"value": "1002", "value_cell": "A2"}]
        result = prep.check_summary(rows, "1001.70 ± 0.38")
        self.assertTrue(result["log_linear_mean_upper_sd_matches"])

    def test_singletons_have_no_invented_sd(self):
        row = {"numeric_value": 3.0, "measurement_state": "observed_finite", "source_cell": "A1"}
        self.assertIsNone(prep.describe([row])["sd_of_available_reported_values"])
        self.assertEqual(prep.check_summary([{"value": "3", "value_cell": "A1"}], "3 ± 0")["diagnostic_status"], "not_tested_singleton")

    def test_formula_cells_are_rejected_and_metadata_does_not_leak_sections(self):
        sheet = openpyxl.Workbook().active
        sheet["A1"] = "=1+1"
        with self.assertRaisesRegex(ValueError, "formula requires reviewed handling"):
            extract.numeric_value(sheet["A1"])
        sheet["A1"], sheet["B1"], sheet["B2"] = "Iteration 1", "5L", "10.1"
        records = list(extract.metadata(sheet, 1, 2, {1: "experiment", 2: "mouse_id"}))
        self.assertEqual(records[1][1]["experiment"], "Iteration 1")
        self.assertEqual(records[1][2]["experiment"], "A1")
        with self.assertRaisesRegex(ValueError, "no anchor"):
            list(extract.metadata(sheet, 2, 2, {1: "experiment", 2: "mouse_id"}))
        sheet["B2"] = None
        with self.assertRaisesRegex(ValueError, "missing identifier"):
            list(extract.metadata(sheet, 1, 2, {1: "experiment", 2: "mouse_id"}))


class PublishedPreparationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.base = Path(cls.temp.name)
        original = extract.OUT
        cls.addClassCleanup(setattr, extract, "OUT", original)
        extract.OUT = cls.base / "extracted"
        cls.before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in extract.RAW.iterdir() if p.is_file()}
        extract.build()
        cls.output = cls.base / "prepared"
        cls.report = prep.build(extract.OUT, cls.output)

    def rows(self, name):
        return prep.read_csv(self.output / (name + ".csv"))

    def test_all_measurements_retained_with_unique_keys_and_distinct_states(self):
        rows = self.rows("observations")
        self.assertEqual(len(rows), 5440)
        self.assertEqual(len({r["observation_id"] for r in rows}), 5440)
        self.assertEqual(self.report["states"], {"observed_finite": 5397, "missing_source_cell": 35, "source_excel_error": 8})
        self.assertEqual(sum("log_table_zero_semantics_unknown" in r["flags"] for r in rows), 512)
        self.assertEqual(sum("reported_negative_normalized_OD" in r["flags"] for r in rows), 4)

    def test_provided_log_cells_and_alternate_linear_values_retained(self):
        book = openpyxl.load_workbook(extract.RAW / "pnas.2521522122.sd09.xlsx", data_only=False)
        sheet = book["Dataset S9"]
        rows = [r for r in self.rows("observations") if r["source_section"] == "s9a_copri_fecal"]
        for row in rows:
            for field, cell_field in [("raw_value", "source_cell"), ("alternate_reported_linear_value", "alternate_linear_cell")]:
                value = sheet[row[cell_field]].value
                self.assertEqual(row[field], "" if value is None else str(value))
        self.assertEqual([r["source_cell"] for r in rows if r["measurement_state"] == "missing_source_cell"], ["M35"])

    def test_published_summary_formulas_differ_between_specimens(self):
        counts = self.report["published_summary_checks"]
        for section, tested, formula in [("b", 68, "log_linear_mean_upper_sd_matches"), ("d", 65, "log_linear_mean_upper_sd_matches"), ("c", 27, "mean_log_sd_log_matches"), ("e", 26, "mean_log_sd_log_matches")]:
            self.assertEqual(counts[f"s9{section}_community"]["tested_complete_zero_free"], tested)
            self.assertEqual(counts[f"s9{section}_community"][formula], tested)
        diagnostic_rows = self.rows("published_summary_diagnostics")
        self.assertEqual(len(diagnostic_rows), 274)
        self.assertEqual(sum(r["diagnostic_status"] == "not_tested_flagged_group" for r in diagnostic_rows), 88)

    def test_missing_replicate_is_disclosed_in_descriptive_contrast(self):
        rows = self.rows("descriptive_contrasts")
        target = next(r for r in rows if r["experiment"] == "Iteration 2" and r["phase_label"] == "Weaning + 4 weeks")
        self.assertEqual((target["treatment_n_slots"], target["treatment_n_finite"]), ("10", "9"))
        self.assertIn("M35", json.loads(target["treatment_cells"]))
        self.assertEqual(target["measurement_scale"], "published_log10_genome_units_per_g")
        self.assertEqual(len(rows), 14)
        self.assertFalse(any("success" in k or "p_value" in k for k in target))

    def test_edited_derived_csv_is_rejected(self):
        path = extract.OUT / "s1c_growth.csv"
        original = path.read_bytes()
        try:
            path.write_bytes(original + b"\n")
            with self.assertRaisesRegex(ValueError, "extracted CSV checksum mismatch"):
                prep.build(extract.OUT, self.base / "reject")
        finally:
            path.write_bytes(original)

    def test_repeat_outputs_are_identical_and_originals_unchanged(self):
        second = self.base / "second"
        prep.build(extract.OUT, second)
        self.assertEqual({p.name: p.read_bytes() for p in self.output.iterdir()}, {p.name: p.read_bytes() for p in second.iterdir()})
        self.assertEqual(self.before, {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in extract.RAW.iterdir() if p.is_file()})


if __name__ == "__main__":
    unittest.main()
