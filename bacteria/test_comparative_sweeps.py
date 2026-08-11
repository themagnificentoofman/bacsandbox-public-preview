import csv
import random
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

import bacteria.core_sim as core_sim
from run_comparative_sweeps import (
    CONDITION_HOMOGENEOUS,
    CONDITION_SPATIAL_NICHE,
    CSV_SCHEMAS,
    DEFAULT_OUTPUT_DIRECTORY,
    DEFAULT_SPECIES_DATA_PATH,
    EXTINCTION_FIELDS,
    HOMOGENEOUS_CONTROL,
    SPATIAL_NICHE_CONTROL,
    ConditionRun,
    SweepConfig,
    _extinction_rows,
    applied_spatial_control,
    jaccard_community_distance,
    run_comparative_sweep,
    write_sweep_outputs,
)


SPECIES_FIELDS = (
    "name",
    "init",
    "rep",
    "aging",
    "so",
    "sf",
    "speed",
    "ph",
    "food_capacity",
    "metabolic_rate",
    "consumption_speed",
    "greediness",
    "t6ss_activity",
    "bacteriocin_production",
    "immunity_factor",
    "rep_thresh",
    "share_thresh",
    "predation_thresh",
    "base_size",
    "sensor_radius",
    "temp_min",
    "temp_max",
    "acetate_secretion",
    "acetate_uptake",
    "baranyi_q0",
    "persister_switch",
    "persister_enter_thresh",
    "persister_exit_thresh",
    "phage_resistance",
    "hgt_competence",
    "aff_calorie",
    "aff_glycan",
    "aff_mdcf",
)


def write_small_species_data(path: Path) -> None:
    base = {field: "0" for field in SPECIES_FIELDS}
    base.update({
        "ph": "7",
        "food_capacity": "100",
        "rep_thresh": "0.8",
        "share_thresh": "0.6",
        "predation_thresh": "0.9",
        "base_size": "1",
        "sensor_radius": "200",
        "temp_min": "4",
        "temp_max": "47",
        "baranyi_q0": "0.05",
        "persister_enter_thresh": "0.2",
        "persister_exit_thresh": "0.5",
        "aff_calorie": "1",
    })
    rare = dict(base, name="Rare survivor", init="2")
    absent = dict(base, name="Not inoculated", init="0")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SPECIES_FIELDS)
        writer.writeheader()
        writer.writerows((rare, absent))


class TestSpeciesAffinityAudit(unittest.TestCase):
    def test_mag_tiers_and_k12_zero_are_explicit(self):
        with DEFAULT_SPECIES_DATA_PATH.open(newline="", encoding="utf-8") as handle:
            rows = {row["name"]: row for row in csv.DictReader(handle)}
        self.assertEqual(
            rows["Prevotella copri MAG Bg0018 (BgF5_2-like model)"]["aff_mdcf"],
            "3",
        )
        self.assertEqual(
            rows["Prevotella copri MAG Bg0019 (BgF5_2-like model)"]["aff_mdcf"],
            "3",
        )
        self.assertEqual(float(rows["Escherichia coli K-12"]["aff_mdcf"]), 0.0)
        self.assertEqual(int(rows["Escherichia coli K-12"]["init"]), 0)


class TestJaccardCommunityDistance(unittest.TestCase):
    def test_presence_absence_definition_and_bounds(self):
        distance = jaccard_community_distance(
            {"a": 4, "b": 1, "absent": 0},
            {"b": 9, "c": 2, "absent": 0},
        )
        self.assertAlmostEqual(distance, 2.0 / 3.0)
        self.assertGreaterEqual(distance, 0.0)
        self.assertLessEqual(distance, 1.0)
        self.assertEqual(jaccard_community_distance({}, {}), 0.0)
        self.assertEqual(jaccard_community_distance({"a": 1}, {}), 1.0)


class TestSpatialControlDefinitions(unittest.TestCase):
    def test_homogeneous_and_spatial_niche_are_not_aliases(self):
        self.assertFalse(HOMOGENEOUS_CONTROL.partition_enabled)
        self.assertEqual(HOMOGENEOUS_CONTROL.mucus_motility_factor, 1.0)
        self.assertEqual(
            HOMOGENEOUS_CONTROL.lumen_food_washout_rate_s_inv,
            HOMOGENEOUS_CONTROL.mucus_food_washout_rate_s_inv,
        )
        self.assertEqual(HOMOGENEOUS_CONTROL.lumen_competition_cull_weight, 1.0)

        self.assertTrue(SPATIAL_NICHE_CONTROL.partition_enabled)
        self.assertLess(SPATIAL_NICHE_CONTROL.mucus_layer_radius_fraction, 1.0)
        self.assertLess(SPATIAL_NICHE_CONTROL.mucus_motility_factor, 1.0)
        self.assertGreater(
            SPATIAL_NICHE_CONTROL.lumen_food_washout_rate_s_inv,
            SPATIAL_NICHE_CONTROL.mucus_food_washout_rate_s_inv,
        )
        self.assertGreater(SPATIAL_NICHE_CONTROL.lumen_competition_cull_weight, 1.0)

    def test_control_context_restores_production_globals(self):
        original = core_sim.MUCUS_LAYER_RADIUS_FRACTION
        with applied_spatial_control(HOMOGENEOUS_CONTROL):
            self.assertEqual(core_sim.MUCUS_LAYER_RADIUS_FRACTION, 1.0)
            self.assertEqual(core_sim.MUCUS_MOTILITY_FACTOR, 1.0)
        self.assertEqual(core_sim.MUCUS_LAYER_RADIUS_FRACTION, original)


class TestComparativeSweep(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary_directory = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary_directory.name)
        cls.species_path = cls.root / "species.csv"
        write_small_species_data(cls.species_path)
        cls.config = SweepConfig(
            species_data_path=cls.species_path,
            output_directory=cls.root / "outputs",
            seeds=(17,),
            diets=("malnutrition",),
            steps=2,
            dt_seconds=0.01,
            world_width_um=120,
            world_height_um=80,
            carrying_capacity_cells=20,
            drop_frequency_s_inv=0.0,
            drop_range_pg=(0.0, 0.0),
            backend="python",
            low_abundance_max_initial_count=2,
            sweep_id="deterministic_test",
        )
        cls.result = run_comparative_sweep(cls.config)

    @classmethod
    def tearDownClass(cls):
        cls.temporary_directory.cleanup()

    def test_default_contract_is_one_thousand_steps_and_scratch_output(self):
        defaults = SweepConfig()
        self.assertEqual(defaults.steps, 1_000)
        self.assertEqual(defaults.output_directory, DEFAULT_OUTPUT_DIRECTORY)

    def test_paired_seed_is_reproducible_and_restores_global_rngs(self):
        random.seed(801)
        np.random.seed(802)
        python_state = repr(random.getstate())
        numpy_state = repr(np.random.get_state())

        repeated = run_comparative_sweep(self.config)

        self.assertEqual(repeated, self.result)
        self.assertEqual(repr(random.getstate()), python_state)
        self.assertEqual(repr(np.random.get_state()), numpy_state)

        metadata = {row["condition"]: row for row in repeated.metadata_rows}
        self.assertEqual(
            metadata[CONDITION_HOMOGENEOUS]["founder_state_sha256"],
            metadata[CONDITION_SPATIAL_NICHE]["founder_state_sha256"],
        )
        self.assertEqual(
            metadata[CONDITION_HOMOGENEOUS]["post_founder_rng_state_sha256"],
            metadata[CONDITION_SPATIAL_NICHE]["post_founder_rng_state_sha256"],
        )
        self.assertEqual(metadata[CONDITION_HOMOGENEOUS]["step_count"], 2)
        self.assertEqual(metadata[CONDITION_HOMOGENEOUS]["timepoint_count"], 3)
        self.assertFalse(metadata[CONDITION_HOMOGENEOUS]["reduced_model_active"])
        self.assertFalse(
            metadata[CONDITION_HOMOGENEOUS]["fba_extinction_dynamics_validated"]
        )
        self.assertEqual(
            metadata[CONDITION_HOMOGENEOUS]["ecm_backend_used"],
            "not_exercised",
        )
        self.assertEqual(
            metadata[CONDITION_SPATIAL_NICHE]["ecm_backend_used"],
            "not_exercised",
        )

    def test_controls_are_active_during_each_condition_update(self):
        observed = []
        original_update = core_sim.Simulation.update

        def observe_update(simulation, *args, **kwargs):
            observed.append((
                core_sim.MUCUS_LAYER_RADIUS_FRACTION,
                core_sim.MUCUS_MOTILITY_FACTOR,
                core_sim.LUMEN_FOOD_WASHOUT_RATE_PER_SEC,
                core_sim.MUCUS_FOOD_WASHOUT_RATE_PER_SEC,
                core_sim.MUCUS_DROP_AMOUNT_FACTOR,
                core_sim.LUMEN_COMPETITION_CULL_WEIGHT,
            ))
            return original_update(simulation, *args, **kwargs)

        with mock.patch.object(core_sim.Simulation, 'update', new=observe_update):
            run_comparative_sweep(self.config)

        homogeneous_expected = (
            HOMOGENEOUS_CONTROL.mucus_layer_radius_fraction,
            HOMOGENEOUS_CONTROL.mucus_motility_factor,
            HOMOGENEOUS_CONTROL.lumen_food_washout_rate_s_inv,
            HOMOGENEOUS_CONTROL.mucus_food_washout_rate_s_inv,
            HOMOGENEOUS_CONTROL.mucus_drop_amount_factor,
            HOMOGENEOUS_CONTROL.lumen_competition_cull_weight,
        )
        spatial_expected = (
            SPATIAL_NICHE_CONTROL.mucus_layer_radius_fraction,
            SPATIAL_NICHE_CONTROL.mucus_motility_factor,
            SPATIAL_NICHE_CONTROL.lumen_food_washout_rate_s_inv,
            SPATIAL_NICHE_CONTROL.mucus_food_washout_rate_s_inv,
            SPATIAL_NICHE_CONTROL.mucus_drop_amount_factor,
            SPATIAL_NICHE_CONTROL.lumen_competition_cull_weight,
        )
        self.assertEqual(observed[:self.config.steps], [homogeneous_expected] * self.config.steps)
        self.assertEqual(observed[self.config.steps:], [spatial_expected] * self.config.steps)

    def test_surviving_low_abundance_strain_is_right_censored(self):
        self.assertEqual(len(self.result.extinction_rows), 2)
        for row in self.result.extinction_rows:
            self.assertEqual(tuple(row), EXTINCTION_FIELDS)
            self.assertEqual(row["species"], "Rare survivor")
            self.assertFalse(row["event_observed"])
            self.assertEqual(row["status"], "right_censored")
            self.assertAlmostEqual(row["time_to_event_seconds"], 0.02)
            self.assertAlmostEqual(row["censor_time_seconds"], 0.02)

    def test_observed_extinction_uses_first_crossing_time(self):
        run = ConditionRun(
            condition=CONDITION_HOMOGENEOUS,
            seed=17,
            diet="malnutrition",
            founder_state_sha256="founder",
            post_founder_rng_state_sha256="rng",
            elapsed_seconds=(0.0, 0.5, 1.0),
            population_counts=(
                {"Rare": 1},
                {"Rare": 0},
                {"Rare": 0},
            ),
            initial_counts={"Rare": 1},
            extinction_times={"Rare": 0.5},
            physics_backend_used="python",
            ecm_backend_used="python",
        )
        rows = list(_extinction_rows(self.config, "pair", run))
        self.assertEqual(len(rows), 1)
        self.assertTrue(rows[0]["event_observed"])
        self.assertEqual(rows[0]["status"], "extinct")
        self.assertEqual(rows[0]["time_to_event_seconds"], 0.5)
        self.assertEqual(rows[0]["censor_time_seconds"], 1.0)

    def test_csv_schemas_and_non_overwrite_behavior(self):
        first_paths = write_sweep_outputs(
            self.result,
            self.config.output_directory,
            directory_name="fixed_run",
        )
        first_bytes = {name: path.read_bytes() for name, path in first_paths.items()}
        second_paths = write_sweep_outputs(
            self.result,
            self.config.output_directory,
            directory_name="fixed_run",
        )

        self.assertNotEqual(
            next(iter(first_paths.values())).parent,
            next(iter(second_paths.values())).parent,
        )
        for filename, fields in CSV_SCHEMAS.items():
            with first_paths[filename].open(newline="", encoding="utf-8") as handle:
                header = tuple(next(csv.reader(handle)))
            self.assertEqual(header, fields)
            self.assertEqual(first_paths[filename].read_bytes(), first_bytes[filename])
            self.assertTrue(second_paths[filename].exists())


if __name__ == "__main__":
    unittest.main()