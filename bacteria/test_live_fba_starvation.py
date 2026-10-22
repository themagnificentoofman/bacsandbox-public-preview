"""Synthetic lifecycle checks; these do not calibrate or validate microbial biology."""

from dataclasses import replace
import logging
import math
import random
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from bacteria.core_sim import Bacterium, BacteriaSpecies, Simulation
from bacteria.diet_translator import DietTranslationEngine
from bacteria.ecm import SubstrateLayerSpec
from bacteria.metabolic_integration import (
    BiomassMaintenanceFbaLifecyclePolicy, DietBolus, FbaAgentState,
    FbaConservationError, FbaIntegrationAdapter, FbaSpatialConfig,
    FbaTickDisposition, GrowthOnlyFbaLifecyclePolicy,
    FBA_LIFECYCLE_STAT_FIELDS,
)
from bacteria.metabolism import (
    MetabolicOrchestrator, SUBSTRATE_ORDER, build_neutral_phenotype_config,
    load_metabolic_model, load_metabolic_priors,
)


def make_fixture(*, policy=None, backend="python"):
    """Explicit synthetic fast kinetics keep lifecycle checks small and observable."""
    species = BacteriaSpecies("synthetic", 0, 0.1, 0.01, 0, 1, 0, 7, (1, 2, 3))
    agent = Bacterium(species, 50, 50, size=1)
    agent.speed = 0
    species.members = [agent]
    if policy is None:
        policy = BiomassMaintenanceFbaLifecyclePolicy(
            initial_biomass_gdw=1e-9, maintenance_rate_s=0.1,
            starvation_biomass_fraction=0.5,
        )
    config = build_neutral_phenotype_config(
        np.array([[0.5, 0.0, 0.0]], dtype=np.float64),
        np.full((1, 3), 1e-4, dtype=np.float64),
    )
    logger = logging.getLogger("test-live-fba-starvation")
    logger.disabled = True
    events = []
    adapter = FbaIntegrationAdapter(
        MetabolicOrchestrator(mode="fba", backend=backend, phenotype_config=config,
                              model=load_metabolic_model(), priors=load_metabolic_priors()),
        DietTranslationEngine(),
        FbaSpatialConfig(
            environment_volume_l=1e-3,
            substrate_layers=tuple(
                SubstrateLayerSpec(name, 0, 0, 0, "Own (synthetic test)")
                for name in SUBSTRATE_ORDER
            ),
            cache_enabled=False,
        ),
        {species: 0}, policy, event_sink=events.append, logger=logger,
    )
    sim = Simulation(
        [species], dict(ph=7, o2=0.21, co2=0, temperature=37,
                        drop_frequency=0, drop_range=(0, 0)),
        cap=100, seed=2401, ecm_backend="python", metabolic_adapter=adapter,
    )
    sim._populated = True
    sim.running = True
    sim._ensure_petri_geometry(100, 100)
    agent.x, agent.y = sim.petri_center
    adapter.state_store.reconcile((agent,))
    return sim, species, agent, adapter, events


def feed(sim, agent):
    sim.queue_diet_bolus(DietBolus(
        "test-feed", "mdcf_2", 1.0,
        ((sim.ecm.world_to_cell(agent.x, agent.y), 1.0),),
    ))


class TestMaintenancePolicy(unittest.TestCase):
    def policy(self):
        return BiomassMaintenanceFbaLifecyclePolicy(
            initial_biomass_gdw=1e-9, maintenance_rate_s=0.1,
            starvation_biomass_fraction=0.5,
        )

    def test_parameters_are_explicit_finite_and_non_degenerate(self):
        with self.assertRaises(TypeError):
            BiomassMaintenanceFbaLifecyclePolicy(initial_biomass_gdw=1e-9)
        for field, values in {
            "maintenance_rate_s": (0, -1, math.inf, math.nan, True),
            "starvation_biomass_fraction": (0, 1, -1, math.inf, math.nan, True),
            "initial_biomass_gdw": (0, math.inf, math.nan, True),
        }.items():
            for value in values:
                with self.subTest(field=field, value=value):
                    with self.assertRaises((ValueError, TypeError)):
                        replace(self.policy(), **{field: value})

    def test_exact_no_feed_decay_and_separate_gross_net_ledger(self):
        policy = self.policy()
        state = FbaAgentState(biomass_gdw=1.2e-9, biomass_delta_gdw=0.2e-9,
                              substrate_consumed_mol=(3e-12, 0, 0))
        after = policy.advance_metabolic_state(state, dt_seconds=2)
        self.assertAlmostEqual(after.biomass_gdw, 1.2e-9 * math.exp(-0.2), delta=1e-24)
        self.assertAlmostEqual(after.biomass_gdw + after.maintenance_loss_gdw,
                               state.biomass_gdw, delta=1e-24)
        self.assertEqual(after.biomass_delta_gdw, state.biomass_delta_gdw)
        self.assertEqual(after.substrate_consumed_mol, state.substrate_consumed_mol)
        self.assertAlmostEqual(after.biomass_gdw - 1e-9,
                               after.net_biomass_delta_gdw, delta=1e-24)
        zero = policy.advance_metabolic_state(state, dt_seconds=0)
        self.assertEqual(zero.biomass_gdw, state.biomass_gdw)
        self.assertEqual(zero.maintenance_loss_gdw, 0)

    def test_division_partitions_all_accounted_amounts(self):
        state = FbaAgentState(2.1e-9, 0.3, 0.2e-9, (1e-12, 2e-12, 3e-12), 4e-12, 0.1e-9)
        parent, child = self.policy().partition_for_fission(object(), object(), state)
        for name in ("biomass_gdw", "biomass_delta_gdw", "maintenance_loss_gdw",
                     "acetate_secreted_mol", "net_biomass_delta_gdw"):
            self.assertEqual(getattr(parent, name) + getattr(child, name), getattr(state, name))
        self.assertEqual(tuple(a + b for a, b in zip(parent.substrate_consumed_mol,
                                                   child.substrate_consumed_mol)),
                         state.substrate_consumed_mol)

    def test_large_loss_does_not_cancel_positive_surviving_biomass(self):
        policy = replace(self.policy(), maintenance_rate_s=40,
                         starvation_biomass_fraction=1e-30)
        state = policy.advance_metabolic_state(FbaAgentState(1e-9), dt_seconds=1)
        self.assertEqual(state.biomass_gdw, 1e-9 * math.exp(-40))
        physiology = policy.evaluate(object(), state, dt_seconds=1,
                                     environmental_fitness=1, carrying_capacity_factor=1)
        self.assertFalse(physiology.starved)


class TestLiveStarvation(unittest.TestCase):
    def test_no_feed_extinction_is_first_completed_threshold_tick(self):
        crossing = math.log(2) / 0.1
        for dt in (0.2, 0.1):
            sim, species, agent, adapter, events = make_fixture()
            total_loss = 0.0
            last_before_death = None
            for _ in range(math.ceil(crossing / dt) + 1):
                if not species.members:
                    break
                last_before_death = adapter.state_store.state_for(agent).biomass_gdw
                result = sim.update(100, 100, dt)
                total_loss += result.maintenance_loss_gdw
                self.assertEqual(result.substrate_consumed_mol, (0, 0, 0))
                self.assertEqual(result.biomass_delta_gdw, 0)
            self.assertEqual(len(species.members), 0)
            self.assertEqual(species.total_deaths, 1)
            self.assertEqual(species.total_births, 0)
            self.assertEqual(sum(sim.stats['fba_starvation_deaths']), 1)
            self.assertEqual(sum(sim.stats['fba_capacity_deaths']), 0)
            self.assertGreaterEqual(species.extinction_time, crossing)
            self.assertLessEqual(species.extinction_time, crossing + dt + 1e-12)
            self.assertGreater(last_before_death, adapter.lifecycle_policy.starvation_biomass_gdw)
            retained_at_death = 1e-9 * math.exp(-0.1 * species.extinction_time)
            self.assertAlmostEqual(sum(sim.stats['fba_starvation_removed_gdw']),
                                   retained_at_death, delta=1e-23)
            self.assertAlmostEqual(total_loss + retained_at_death, 1e-9, delta=1e-23)
            self.assertEqual(adapter.state_store.snapshot(), {})
            empty = sim.update(100, 100, dt)
            self.assertEqual(empty.disposition, FbaTickDisposition.NO_AGENTS)
            self.assertEqual(empty.maintenance_loss_gdw, 0)
            self.assertTrue(all(e.disposition != FbaTickDisposition.ROLLED_BACK for e in events))

    def test_well_fed_cell_recovers_then_really_divides(self):
        sim, species, agent, adapter, _ = make_fixture()
        sim.update(100, 100, 1)
        depleted = adapter.state_store.state_for(agent).biomass_gdw
        self.assertLess(depleted, 1e-9)
        feed(sim, agent)
        gross = loss = 0.0
        for _ in range(60):
            result = sim.update(100, 100, 0.1)
            gross += result.biomass_delta_gdw
            loss += result.maintenance_loss_gdw
            if species.total_births:
                break
        self.assertEqual(species.total_births, 1)  # No vacuous daughter check.
        self.assertEqual(species.total_deaths, 0)
        self.assertEqual(len(species.members), 2)
        total = sum(adapter.state_store.state_for(a).biomass_gdw for a in species.members)
        self.assertAlmostEqual(total, depleted + gross - loss, delta=1e-23)
        self.assertTrue(all(adapter.state_store.state_for(a).biomass_gdw >= 1e-9
                            for a in species.members))

    def test_growth_only_policy_retains_zero_food_survival(self):
        sim, species, agent, adapter, _ = make_fixture(
            policy=GrowthOnlyFbaLifecyclePolicy(initial_biomass_gdw=1e-9))
        result = sim.update(100, 100, 100)
        self.assertEqual(species.members, [agent])
        self.assertEqual(adapter.state_store.state_for(agent).biomass_gdw, 1e-9)
        self.assertEqual(result.maintenance_loss_gdw, 0)

    def test_consumption_and_growth_remain_recorded_when_cell_dies_same_tick(self):
        policy = BiomassMaintenanceFbaLifecyclePolicy(
            initial_biomass_gdw=1e-9, maintenance_rate_s=10,
            starvation_biomass_fraction=0.5,
        )
        sim, species, agent, adapter, _ = make_fixture(policy=policy)
        feed(sim, agent)
        result = sim.update(100, 100, 0.1)
        self.assertEqual(species.members, [])
        self.assertEqual(species.total_deaths, 1)
        self.assertEqual(adapter.state_store.snapshot(), {})
        self.assertGreater(sum(result.substrate_consumed_mol), 0)
        self.assertGreater(result.biomass_delta_gdw, 0)
        self.assertGreater(result.maintenance_loss_gdw, result.biomass_delta_gdw)
        self.assertLess(result.net_biomass_delta_gdw, 0)
        self.assertEqual(sim.species_step_stats[-1]['fba_events'][0]['event'], 'starvation')
        self.assertAlmostEqual(
            1e-9 + result.biomass_delta_gdw - result.maintenance_loss_gdw,
            sim.stats['fba_starvation_removed_gdw'][-1], delta=1e-24,
        )

    def test_capacity_removal_is_separate_and_csv_retains_ledger_and_events(self):
        sim, species, agent, adapter, _ = make_fixture()
        species.members.append(Bacterium(species, agent.x, agent.y, size=1))
        sim.cap = 1
        result = sim.update(100, 100, 0.1)
        self.assertEqual(len(species.members), 1)
        self.assertEqual(sim.stats['fba_starvation_deaths'], [0])
        self.assertEqual(sim.stats['fba_capacity_deaths'], [1])
        self.assertGreater(sim.stats['fba_capacity_removed_gdw'][0], 0)
        self.assertAlmostEqual(
            sim.stats['fba_biomass_gdw'][0] + sim.stats['fba_capacity_removed_gdw'][0]
            + result.maintenance_loss_gdw, 2e-9, delta=1e-24,
        )
        self.assertEqual(sim.species_step_stats[0]['fba_events'][0]['event'], 'capacity')
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'stats.csv'
            with patch('bacteria.core_sim.EXPORT_FILE', str(path)):
                sim.export()
            restored, _, _, _, _ = make_fixture()
            restored.import_stats(str(path))
            for name in FBA_LIFECYCLE_STAT_FIELDS:
                self.assertEqual(restored.stats[name], sim.stats[name])
            self.assertEqual(restored.species_step_stats[0]['fba_events'],
                             sim.species_step_stats[0]['fba_events'])

    def test_partition_failure_cannot_create_biomass_or_leave_partial_events(self):
        sim, species, agent, adapter, events = make_fixture()
        adapter.state_store.restore({agent: FbaAgentState(2.2e-9)})
        before = adapter.state_store.snapshot()
        with patch.object(BiomassMaintenanceFbaLifecyclePolicy, 'partition_for_fission',
                          side_effect=lambda parent, child, state: (state, state)):
            with self.assertRaises(FbaConservationError):
                sim.update(100, 100, 0.1)
        self.assertEqual(adapter.state_store.snapshot(), before)
        self.assertEqual(sim.species_step_stats, [])
        self.assertEqual(species.members, [agent])
        self.assertEqual(sim.step_idx, 0)
        self.assertEqual(events[-1].disposition, FbaTickDisposition.ROLLED_BACK)

    def test_bad_lifecycle_sink_rolls_back_feed_solver_rng_and_tick(self):
        sim, species, agent, adapter, events = make_fixture()
        feed(sim, agent)
        before_state = adapter.state_store.snapshot()
        before_grid = sim.ecm.snapshot_substrate_layers()
        before_rng = random.getstate()
        before_queue = tuple(sim._pending_diet_boluses)
        with patch.object(BiomassMaintenanceFbaLifecyclePolicy, "advance_metabolic_state",
                          side_effect=lambda state, **kw: replace(state, biomass_gdw=state.biomass_gdw * 2)):
            with self.assertRaises(FbaConservationError):
                sim.update(100, 100, 0.1)
        self.assertEqual(adapter.state_store.snapshot(), before_state)
        self.assertEqual(sim.ecm.snapshot_substrate_layers(), before_grid)
        self.assertEqual(random.getstate(), before_rng)
        self.assertEqual(tuple(sim._pending_diet_boluses), before_queue)
        self.assertEqual(sim.step_idx, 0)
        self.assertEqual(species.members, [agent])
        self.assertEqual(species.total_deaths, 0)
        self.assertEqual(events[-1].disposition, FbaTickDisposition.ROLLED_BACK)
        result = sim.update(100, 100, 0.1)
        self.assertEqual(result.disposition, FbaTickDisposition.COMMITTED)
        self.assertEqual(sim.step_idx, 1)


if __name__ == "__main__":
    unittest.main()
