"""Research parameter inventory guards; production loader behavior stays intact."""

import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

from audit_species_parameters import audit_species_parameters, main


class ParameterAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "profiles.csv"

    def write_profile(self, row):
        with self.path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)

    def test_missing_blank_zero_and_alias_origins_are_distinct(self):
        self.write_profile({
            "name": "test", "init": "2", "aff_calorie": "0",
            "temp_min": "", "aggression": "0.25", "defense": "0.75",
        })
        original = self.path.read_bytes()
        report = audit_species_parameters(self.path)
        parameters = {row["field"]: row for row in report["profiles"][0]["parameters"]}
        self.assertEqual(parameters["aff_calorie"]["value_origin"], "csv_explicit")
        self.assertEqual(parameters["aff_calorie"]["constructor_value"], 0)
        self.assertEqual(parameters["temp_min"]["value_origin"], "compatibility_default")
        self.assertEqual(parameters["temp_min"]["loader_value"], "4.0")
        self.assertEqual(parameters["t6ss_activity"]["value_origin"], "legacy_alias")
        self.assertEqual(parameters["t6ss_activity"]["source_field"], "aggression")
        self.assertEqual(parameters["t6ss_activity"]["constructor_value"], 0.25)
        self.assertEqual(parameters["immunity_factor"]["constructor_value"], 0.75)
        self.assertEqual(self.path.read_bytes(), original)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(["--species-data", str(self.path), "--require-explicit"]), 1)

    def test_fixed_temperature_is_not_unit_or_parameter_independent(self):
        self.write_profile({"name": "test", "init": "1", "temp_min": "4", "temp_max": "47"})
        baseline = audit_species_parameters(self.path, 37)["profiles"][0]["temperature_factor"]
        self.assertGreater(baseline, 0.98)
        self.assertLess(baseline, 0.999)
        self.write_profile({"name": "test", "init": "1", "temp_min": "4", "temp_max": "42"})
        shifted = audit_species_parameters(self.path, 37)["profiles"][0]["temperature_factor"]
        self.assertLess(shifted, 0.7)
        self.assertLess(shifted, baseline)

    def test_constructor_clamping_is_visible(self):
        self.write_profile({"name": "test", "init": "2.9", "baranyi_q0": "0"})
        report = audit_species_parameters(self.path)
        parameters = {row["field"]: row for row in report["profiles"][0]["parameters"]}
        self.assertEqual(report["total_founders"], 2)
        self.assertTrue(parameters["init"]["constructor_changed_value"])
        self.assertTrue(parameters["baranyi_q0"]["constructor_changed_value"])
        self.assertGreater(parameters["baranyi_q0"]["constructor_value"], 0)
        self.assertEqual(parameters["baranyi_q0"]["empirical_provenance"], "not_recorded_per_value")

    def test_invalid_numbers_are_flagged_without_creating_false_totals(self):
        for bad_value in ("nonsense", "nan", "inf", "-inf"):
            with self.subTest(bad_value=bad_value):
                self.write_profile({"name": "test", "init": "1", "rep": bad_value})
                report = audit_species_parameters(self.path)
                self.assertTrue(report["numeric_and_structural_issues"])
                self.assertFalse(report["all_profiles_evaluated"])
                self.assertIsNone(report["total_founders"])
                json.dumps(report, allow_nan=False)

    def test_interactive_app_init_parsing_is_not_confused_with_constructor(self):
        for text in ("2.9", "100.0", "1e2"):
            with self.subTest(init=text):
                self.write_profile({"name": "test", "init": text})
                report = audit_species_parameters(self.path)
                profile = report["profiles"][0]
                self.assertEqual(profile["founder_count"], int(float(text)))
                self.assertFalse(profile["interactive_app_init_compatible"])
                self.assertIsNone(profile["interactive_app_init_count"])
                self.assertTrue(any("interactive app" in issue for issue in profile["issues"]))
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(main(["--species-data", str(self.path)]), 1)
        self.write_profile({"name": "test", "init": "2"})
        profile = audit_species_parameters(self.path)["profiles"][0]
        self.assertTrue(profile["interactive_app_init_compatible"])
        self.assertEqual(profile["interactive_app_init_count"], 2)

    def test_csv_locations_account_for_blank_lines_and_multiline_records(self):
        self.path.write_text('name,init\n\n"line one\nline two",2\nnext,3\n', encoding="utf-8")
        report = audit_species_parameters(self.path)
        self.assertEqual([p["csv_row"] for p in report["profiles"]], [2, 3])
        self.assertEqual([p["csv_line_end"] for p in report["profiles"]], [4, 5])

    def test_duplicate_headers_are_rejected(self):
        self.path.write_text("name,init,init\na,2,3\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duplicate CSV headers"):
            audit_species_parameters(self.path)

    def test_output_never_overwrites_prior_audit(self):
        self.write_profile({"name": "test", "init": "1"})
        output = Path(self.temp.name) / "audit.json"
        args = ["--species-data", str(self.path), "--output", str(output)]
        self.assertEqual(main(args), 0)
        original = output.read_bytes()
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main(args), 2)
        self.assertEqual(output.read_bytes(), original)

    def test_shipped_csv_inventory(self):
        path = Path(__file__).resolve().parents[1] / "species_data.csv"
        report = audit_species_parameters(path)
        self.assertEqual(report["loader_field_count"], 33)
        self.assertEqual(report["numeric_field_count"], 32)
        self.assertEqual(report["trait_field_count_excluding_name_and_init"], 31)
        self.assertEqual(report["profile_count"], 7)
        self.assertEqual(report["inoculated_profile_count"], 2)
        self.assertEqual(report["total_founders"], 200)
        self.assertEqual(report["numeric_and_structural_issues"], [])
        self.assertEqual(report["value_origin_counts"]["compatibility_default"], 0)


if __name__ == "__main__":
    unittest.main()
