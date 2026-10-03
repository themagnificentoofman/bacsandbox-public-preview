"""Synthetic parser fixtures only; these are not measurements or benchmark data."""
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path
import unittest

from audit_benchmark_table import audit, read_table


class BenchmarkTableTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "table.csv"
        self.contract = {"required_columns": ["experiment", "sample", "value"],
                         "key_columns": ["experiment", "sample"], "value_columns": ["value"],
                         "missing_tokens": ["NA"], "units": {"value": None},
                         "time_mapping": None, "source_locator": "synthetic parser fixture"}

    def write(self, text):
        self.path.write_text(text)

    def test_literal_missing_zero_negative_and_input_preservation(self):
        self.write("experiment,sample,value\ne1,s1,0\ne1,s2,NA\ne1,s3,-1\ne1,s4,\n")
        before = self.path.read_bytes()
        result = audit(self.path, self.contract)
        self.assertTrue(result["structural_ok"])
        self.assertEqual(result["columns"]["value"], dict(missing=2, finite=2, zero=1, negative=1, nonnumeric=0, nonfinite=0))
        self.assertEqual(result["unknown_units"], ["value"])
        self.assertEqual(result["scientific_readiness"], "not_assessed")
        self.assertEqual(result["source_sha256"], hashlib.sha256(before).hexdigest())
        self.assertEqual(before, self.path.read_bytes())

    def test_bad_headers_and_ragged_rows(self):
        for text in ("", "x,x\n1,2\n", "x,\n1,2\n", "x,y\n1\n", "x,y\n1,2,3\n"):
            with self.subTest(text=text):
                self.write(text)
                with self.assertRaises(ValueError):
                    read_table(self.path)

    def test_duplicate_keys_are_scoped_to_experiment(self):
        self.write("experiment,sample,value\ne1,s1,1\ne2,s1,2\ne1,s1,3\n")
        result = audit(self.path, self.contract)
        self.assertFalse(result["structural_ok"])
        self.assertEqual(result["duplicate_keys"], [{"key": ["e1", "s1"], "first_line": 2, "line": 4}])

    def test_missing_keys(self):
        self.write("experiment,sample,value\ne1,,1\ne1,NA,2\n")
        self.assertEqual(audit(self.path, self.contract)["missing_key_lines"], [2, 3])

    def test_unknown_tokens_not_silently_missing(self):
        self.write("experiment,sample,value\ne1,s1,<LOD\ne1,s2,NaN\ne1,s3,inf\n")
        result = audit(self.path, self.contract)
        self.assertFalse(result["structural_ok"])
        self.assertEqual(result["columns"]["value"]["nonfinite"], 2)
        self.assertEqual(result["invalid_values"][0]["raw"], "<LOD")

    def test_bom_and_multiline_physical_line(self):
        self.write('\ufeffexperiment,sample,value\ne1,"s\n1",bad\n')
        self.assertEqual(audit(self.path, self.contract)["invalid_values"][0]["line"], 3)

    def test_missing_columns_and_ambiguous_contracts(self):
        self.write("experiment,sample,value\ne1,s1,1\n")
        with self.assertRaises(ValueError):
            audit(self.path, [])
        for change in ({"key_columns": []}, {"key_columns": ["absent"]}, {"value_columns": ["sample"]}, {"missing_tokens": "NA"}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                audit(self.path, self.contract | change)

    def test_empty_table_not_ready(self):
        self.write("experiment,sample,value\n")
        self.assertFalse(audit(self.path, self.contract)["structural_ok"])

    def test_join_reports_unmatched_and_duplicate_metadata(self):
        self.write("experiment,sample,value\ne1,s1,1\ne2,s1,2\n")
        meta = Path(self.tmp.name) / "meta.csv"
        meta.write_text("experiment,sample\ne1,s1\ne1,s1\ne3,s2\n")
        result = audit(self.path, self.contract | {"join_columns": ["experiment", "sample"]}, meta)
        self.assertFalse(result["structural_ok"])
        self.assertEqual(result["metadata_join"]["data_keys_without_metadata"], [["e2", "s1"]])
        self.assertEqual(result["metadata_join"]["metadata_keys_without_data"], [["e3", "s2"]])
        self.assertEqual(len(result["metadata_join"]["duplicate_metadata_keys"]), 1)

    def test_many_measurements_can_join_one_sample(self):
        self.write("experiment,sample,organism,value\ne1,s1,o1,1\ne1,s1,o2,2\n")
        meta = Path(self.tmp.name) / "meta.csv"
        meta.write_text("experiment,sample\ne1,s1\n")
        contract = self.contract | {"key_columns": ["experiment", "sample", "organism"], "join_columns": ["experiment", "sample"]}
        self.assertTrue(audit(self.path, contract, meta)["structural_ok"])

    def cli(self, output):
        contract = Path(self.tmp.name) / "contract.json"
        contract.write_text(json.dumps(self.contract))
        script = Path(__file__).resolve().parents[1] / "audit_benchmark_table.py"
        return subprocess.run([sys.executable, str(script), str(self.path), "--contract", str(contract), "--output", str(output)], capture_output=True, text=True)

    def test_cli_refuses_overwrite(self):
        self.write("experiment,sample,value\ne1,s1,1\n")
        output = Path(self.tmp.name) / "result.json"
        output.write_text("preserve this")
        self.assertEqual(self.cli(output).returncode, 2)
        self.assertEqual(output.read_text(), "preserve this")

    def test_cli_rejects_nonfinite_contract_without_partial_output(self):
        self.write("experiment,sample,value\ne1,s1,1\n")
        for i, value in enumerate((float("nan"), float("inf"), float("-inf"))):
            with self.subTest(value=value):
                self.contract["time_mapping"] = {"days": value}
                output = Path(self.tmp.name) / f"result-{i}.json"
                result = self.cli(output)
                self.assertEqual(result.returncode, 2)
                self.assertFalse(output.exists())
                self.assertNotIn("Traceback", result.stderr)

    def test_cli_writes_failure_evidence(self):
        self.write("experiment,sample,value\ne1,s1,<LOD\n")
        output = Path(self.tmp.name) / "result.json"
        self.assertEqual(self.cli(output).returncode, 2)
        self.assertFalse(json.loads(output.read_text())["structural_ok"])

    def test_cli_success_is_structural_only(self):
        self.write("experiment,sample,value\ne1,s1,1\n")
        output = Path(self.tmp.name) / "result.json"
        self.assertEqual(self.cli(output).returncode, 0)
        self.assertEqual(json.loads(output.read_text())["scientific_readiness"], "not_assessed")


if __name__ == "__main__":
    unittest.main()
