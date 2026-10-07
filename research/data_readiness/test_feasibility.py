"""Tests for source-cohort integrity and bounded configuration diagnostics."""
import math
from types import SimpleNamespace
from unittest.mock import patch
import unittest

from research.data_readiness import audit_feasibility as audit, map_benchmark


class RuntimeDependencyTests(unittest.TestCase):
    def provenance(self, *, is_ce, installed, module_version):
        module = SimpleNamespace(version=SimpleNamespace(ver=module_version))
        if is_ce is not None:
            module.IS_CE = is_ce
        versions = {"numpy": "2.3.5", "openpyxl": "3.1.5", **installed}

        def version(name):
            if name not in versions:
                raise audit.importlib.metadata.PackageNotFoundError(name)
            return versions[name]

        with patch.object(audit, "pygame", module), patch.object(
                audit.importlib.metadata, "version", side_effect=version):
            return audit.runtime_dependencies()

    def test_pygame_ce_without_legacy_distribution(self):
        dependencies, provenance = self.provenance(
            is_ce=True, installed={"pygame-ce": "2.5.8"}, module_version="2.5.8")
        self.assertEqual(dependencies, {"numpy": "2.3.5", "pygame-ce": "2.5.8", "openpyxl": "3.1.5"})
        self.assertEqual(provenance["distribution"], "pygame-ce")
        self.assertTrue(provenance["module_matches_distribution_version"])
        self.assertFalse(provenance["coinstalled_distributions"])

    def test_legacy_pygame_without_ce_marker(self):
        dependencies, provenance = self.provenance(
            is_ce=None, installed={"pygame": "2.6.1"}, module_version="2.6.1")
        self.assertEqual(dependencies["pygame"], "2.6.1")
        self.assertNotIn("pygame-ce", dependencies)
        self.assertEqual(provenance["distribution"], "pygame")
        self.assertFalse(provenance["coinstalled_distributions"])

    def test_coinstallation_reports_loaded_provider_and_both_metadata_records(self):
        installed = {"pygame": "2.6.1", "pygame-ce": "2.5.8"}
        for is_ce, provider in ((True, "pygame-ce"), (False, "pygame")):
            with self.subTest(provider=provider):
                dependencies, provenance = self.provenance(
                    is_ce=is_ce, installed=installed, module_version=installed[provider])
                self.assertEqual(provenance["distribution"], provider)
                self.assertEqual(provenance["installed_distributions"], installed)
                self.assertTrue(provenance["coinstalled_distributions"])
                self.assertEqual(dependencies[provider], installed[provider])
                self.assertEqual(len(set(installed).intersection(dependencies)), 1)

    def test_missing_actual_provider_metadata_is_not_replaced_with_other_provider(self):
        with self.assertRaises(audit.importlib.metadata.PackageNotFoundError):
            self.provenance(is_ce=True, installed={"pygame": "2.6.1"}, module_version="2.5.8")

    def test_module_metadata_version_disagreement_is_visible(self):
        _, provenance = self.provenance(
            is_ce=True, installed={"pygame-ce": "2.5.8"}, module_version="2.5.7")
        self.assertEqual(provenance["module_version"], "2.5.7")
        self.assertEqual(provenance["distribution_version"], "2.5.8")
        self.assertFalse(provenance["module_matches_distribution_version"])


class FeasibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mapping = map_benchmark.build(map_benchmark.HERE / "raw/supplements", audit.ROOT / "species_data.csv")
        cls.payloads = audit.build(map_benchmark.HERE / "raw/supplements")

    def test_incomplete_vectors_are_not_silently_renormalized(self):
        result = self.payloads["feasibility_audit.json"]
        self.assertEqual(result["complete_baseline_animals"], 8)
        self.assertEqual(result["complete_baseline_iterations"], {"Iteration 1": 8})
        self.assertEqual(len(result["excluded_incomplete_baselines"]), 31)
        self.assertEqual(sum(len(r["unresolved_cells"]) for r in result["excluded_incomplete_baselines"]), 49)
        self.assertFalse(result["benchmark_execution_ready"])
        self.assertFalse(result["ecological_simulation_executed"])
        self.assertFalse(result["parameters_fitted"])

    def test_resolution_uses_whole_complete_vector_and_validates_inputs(self):
        rows = self.payloads["founder_resolution.csv"]
        self.assertEqual(len(rows), 24)
        minimum_budgets = {r["budget_for_minimum_expected_quota_one"] for r in rows}
        self.assertEqual((min(minimum_budgets), max(minimum_budgets)), (277836, 923032))
        self.assertTrue(all(r["rarest_taxon"] == "Ruminococcus torques" for r in rows))
        self.assertTrue(all(r["minimum_expected_quota"] < 1 for r in rows))
        for r in rows:
            n = r["budget_for_minimum_expected_quota_one"]
            self.assertGreaterEqual(n * r["minimum_fraction"], 1)
            self.assertLess((n - 1) * r["minimum_fraction"], 1)
            self.assertFalse(r["runnable_founders"])
        with self.assertRaises(ValueError):
            audit.founder_resolution(self.mapping["baseline_observations.csv"], [False])
        with self.assertRaises(ValueError):
            audit.founder_resolution(self.mapping["baseline_observations.csv"][:-1], [500])

    def test_actual_diet_interfaces_reject_both_unimplemented_arms(self):
        probes = self.payloads["diet_interface_probes.json"]["requested_diet_probes"]
        self.assertEqual(len(probes), 4)
        self.assertEqual({r["diet"] for r in probes}, set(map_benchmark.DIETS))
        self.assertTrue(all(not r["accepted"] for r in probes))
        self.assertTrue(all(r["reason"] for r in probes))

    def test_common_substrate_route_and_local_scale_counterexample(self):
        probes = self.payloads["structural_probes.json"]
        self.assertEqual(probes["substrate_component_groups"]["mdcf_prebiotic"],
                         ["beta_1_4_galactan", "beta_1_4_mannan"])
        first, second = probes["uptake_scale_counterexample"]
        self.assertNotEqual(first["concentration"], second["concentration"])
        self.assertTrue(math.isclose(first["uptake_bound"], second["uptake_bound"], rel_tol=1e-14, abs_tol=0))
        self.assertTrue(all(not p["starved"] for p in probes["growth_only_policy_probes"]))

    def test_profiles_are_configuration_not_empirical_distinctions(self):
        profiles = self.payloads["profile_provenance.csv"]
        target, = [r for r in profiles if r["name"] == "Segatella copri BgF5_2 (cultured representative)"]
        self.assertEqual(target["non_name_non_init_fields"], 31)
        self.assertEqual(target["fields_equal_to_synthetic_a"], 29)
        self.assertEqual(target["configured_init"], 0)
        self.assertEqual([r["name"] for r in profiles if r["configured_init"] > 0], ["a", "b"])

    def test_physical_day_arithmetic_does_not_select_clock(self):
        rows = self.payloads["time_workload.csv"]
        self.assertEqual(len(rows), 8)
        self.assertEqual({r["ticks_if_literal_seconds"] for r in rows},
                         {8640000, 9504000, 10368000, 24192000, 25056000, 25920000})
        self.assertTrue(all(not r["clock_selected"] and not r["wall_time_measured"] for r in rows))
        self.assertEqual(self.payloads["structural_probes.json"]["default_harness"]["horizon_seconds"], 100.0)


if __name__ == "__main__":
    unittest.main()
