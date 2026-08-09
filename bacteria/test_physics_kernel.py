import copy
import gc
import math
import os
import random
import statistics
import time
import tracemalloc
import unittest
import weakref
from dataclasses import replace
from unittest.mock import patch

import numpy as np

from bacteria.constants import NEUTRAL_PH
from bacteria.core_sim import Bacterium, BacteriaSpecies
from bacteria.physics import (
    EVENT_CONTACT,
    EVENT_HGT_TRANSFER,
    EVENT_T6SS_LETHAL,
    ContactParameters,
    PhysicsBuffers,
    PhysicsOrchestrator,
)

try:
    import bacteria_c_kernel
except ImportError:
    bacteria_c_kernel = None


NATIVE_PHYSICS_AVAILABLE = (
    bacteria_c_kernel is not None
    and hasattr(bacteria_c_kernel, 'process_interactions')
)

# A 30 Hz-equivalent step keeps all contact rates in real elapsed seconds while making
# multi-tick branch accumulation observable; see DOCUMENTATION.md §4.12 and §4.18.
CONTACT_DT_SECONDS = 1.0 / 30.0
FLOAT_ABSOLUTE_TOLERANCE = 1.0e-5
MEDIUM_RANDOM_SEED = 18473
# A half-megabyte ceiling allows allocator bookkeeping noise while rejecting retained per-tick
# contact structures after warmup; see DOCUMENTATION.md §10.2.
ALLOCATION_GROWTH_TOLERANCE_BYTES = 512 * 1024
NATIVE_PHYSICS_STRESS_DEFAULT_ITERATIONS = 100
RUN_PHYSICS_SCALE = os.environ.get('BACSANDBOX_RUN_PHYSICS_SCALE') == '1'
RUN_PHYSICS_BENCHMARK = os.environ.get('BACSANDBOX_RUN_PHYSICS_BENCHMARK') == '1'
SCALE_RANDOM_SEED = 910247
SCALE_AGENTS_PER_SPECIES = 200
SCALE_TICK_COUNT = 500
SCALE_TURNOVER_INTERVAL = 10
# A near-zero positive reserve guarantees that a positive T6SS dose crosses the lethal boundary
# without changing the production damage equation; see DOCUMENTATION.md §4.12.
SCALE_LETHAL_PROBE_FOOD = 1.0e-9
# Seventy percent fullness keeps non-probe cells away from lethal depletion during one contact tick;
# see DOCUMENTATION.md §4.12.
SCALE_STEADY_FULLNESS = 0.7
BENCHMARK_REPEAT_COUNT = 3
BENCHMARK_TICK_COUNT = 8

SCALE_SPECIES_NAMES = (
    'Bacteroides thetaiotaomicron VPI-5482',
    'Escherichia coli K-12',
    'Faecalibacterium prausnitzii A2-165',
    'Bifidobacterium longum NCC 2705',
    'Akkermansia muciniphila MucT',
)

FLOAT_AGENT_ARGUMENTS = (
    'x',
    'y',
    'interaction_radius',
    'base_radius',
    'food_level',
    'food_capacity',
    'metabolic_burden_factor',
)
UINT8_AGENT_ARGUMENTS = ('is_persister', 'phage_resistant', 'event_flags')
FLOAT_SPECIES_ARGUMENTS = (
    't6ss_activity',
    'immunity_factor',
    'hgt_competence',
    'share_threshold',
    'predation_threshold',
)
MUTABLE_AGENT_ARGUMENTS = (
    'food_level',
    'metabolic_burden_factor',
    'phage_resistant',
    'event_flags',
)


def synthetic_species(
    name,
    *,
    t6ss_activity=0.0,
    immunity_factor=0.0,
    hgt_competence=0.0,
    share_threshold=1.0,
    predation_threshold=0.0,
    base_size=1.0,
):
    """Build a contact-only species whose traits are explicitly synthetic test inputs."""
    return BacteriaSpecies(
        name,
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        NEUTRAL_PH,
        (90, 140, 190),
        # A fixed synthetic capacity makes sharing, damage, and predation state easy to align;
        # see DOCUMENTATION.md §4.12.
        food_capacity=10.0,
        metabolic_rate=0.0,
        consumption_speed=0.0,
        t6ss_activity=t6ss_activity,
        bacteriocin_production=0.0,
        immunity_factor=immunity_factor,
        share_thresh=share_threshold,
        predation_thresh=predation_threshold,
        base_size=base_size,
        persister_switch=0.0,
        phage_resistance=0.0,
        hgt_competence=hgt_competence,
    )


def synthetic_cell(
    species,
    x,
    y,
    food_level,
    *,
    base_size=None,
    interaction_radius=None,
    resistant=False,
    persister=False,
):
    """Build a deterministic real Bacterium with contact-owned state set explicitly."""
    radius = species.base_size if base_size is None else base_size
    with patch('bacteria.core_sim.random.gauss', return_value=1.0), \
         patch('bacteria.core_sim.random.random', return_value=0.5), \
         patch('bacteria.core_sim.random.uniform', return_value=0.0):
        cell = Bacterium(
            species,
            x,
            y,
            size=radius,
            initial_food=food_level,
        )
    cell.size = radius if interaction_radius is None else interaction_radius
    cell.phage_resistant = resistant
    cell.is_persister = persister
    cell.metabolic_burden_factor = 1.0
    return cell


def clone_species_state(species):
    """Clone species, member order, and each member's species back-reference exactly."""
    return copy.deepcopy(species)


def ordered_cells(species):
    return [cell for current_species in species for cell in current_species.members]


def object_state(species):
    cells = ordered_cells(species)
    float_state = np.asarray(
        [
            (
                cell.x,
                cell.y,
                cell.size,
                cell.base_size,
                cell.food_level,
                cell.food_capacity,
                cell.metabolic_burden_factor,
            )
            for cell in cells
        ],
        dtype=np.float64,
    ).reshape((-1, 7))
    discrete_state = tuple(
        (cell.phage_resistant, cell.is_persister)
        for cell in cells
    )
    return float_state, discrete_state


def valid_native_arguments(agent_count=2, species_count=2):
    """Return independent contiguous arrays satisfying the complete native ABI contract."""
    if agent_count and not species_count:
        species_ids = np.zeros(agent_count, dtype=np.int32)
    elif species_count:
        # Cycling IDs exercises all species rows without introducing a biological distribution;
        # see DOCUMENTATION.md §10.2.
        species_ids = np.ascontiguousarray(
            np.arange(agent_count, dtype=np.int32) % species_count,
            dtype=np.int32,
        )
    else:
        species_ids = np.empty(0, dtype=np.int32)

    parameters = ContactParameters()
    return {
        # Unit-spaced coordinates form deterministic contact candidates for ABI mutation checks;
        # see DOCUMENTATION.md §4.12.
        'x': np.arange(agent_count, dtype=np.float64),
        'y': np.zeros(agent_count, dtype=np.float64),
        'interaction_radius': np.ones(agent_count, dtype=np.float64),
        'base_radius': np.ones(agent_count, dtype=np.float64),
        'food_level': np.full(agent_count, 5.0, dtype=np.float64),
        'food_capacity': np.full(agent_count, 10.0, dtype=np.float64),
        'metabolic_burden_factor': np.ones(agent_count, dtype=np.float64),
        'species_id': species_ids,
        'is_persister': np.zeros(agent_count, dtype=np.uint8),
        'phage_resistant': np.asarray(
            [1 if index == 0 else 0 for index in range(agent_count)],
            dtype=np.uint8,
        ),
        'event_flags': np.zeros(agent_count, dtype=np.uint8),
        't6ss_activity': np.zeros(species_count, dtype=np.float64),
        'immunity_factor': np.zeros(species_count, dtype=np.float64),
        'hgt_competence': np.ones(species_count, dtype=np.float64),
        'share_threshold': np.ones(species_count, dtype=np.float64),
        'predation_threshold': np.zeros(species_count, dtype=np.float64),
        'dt': CONTACT_DT_SECONDS,
        'hash_cell_size': parameters.hash_cell_size,
        'axis_cutoff': parameters.axis_cutoff,
        'contact_margin': parameters.contact_margin,
        'hgt_transfer_rate': parameters.hgt_transfer_rate,
        'hgt_metabolic_burden': parameters.hgt_metabolic_burden,
        't6ss_damage_scaling': parameters.t6ss_damage_scaling,
        't6ss_max_immunity_defense': parameters.t6ss_max_immunity_defense,
        't6ss_immunity_size_epsilon': parameters.t6ss_immunity_size_epsilon,
        'persister_resistance_factor': parameters.persister_resistance_factor,
        'phagocytosis_size_ratio': parameters.phagocytosis_size_ratio,
        'cell_wet_mass': parameters.cell_wet_mass,
        'cell_dry_mass_fraction': parameters.cell_dry_mass_fraction,
        'labile_necromass_fraction': parameters.labile_necromass_fraction,
        'cytosol_recovery_fraction': parameters.cytosol_recovery_fraction,
        'consumed_prey_sentinel': parameters.consumed_prey_sentinel,
        'random_draw': lambda: 0.5,
    }


def medium_fixture():
    """Construct signed, cross-boundary contact groups covering every interaction family."""
    kin = synthetic_species(
        'Synthetic kin donors',
        t6ss_activity=0.15,
        immunity_factor=0.25,
        hgt_competence=1.0,
        share_threshold=0.4,
    )
    attackers = synthetic_species(
        'Synthetic attackers',
        t6ss_activity=0.8,
        immunity_factor=0.1,
        hgt_competence=0.9,
    )
    persisters = synthetic_species(
        'Synthetic persister recipients',
        t6ss_activity=0.1,
        immunity_factor=0.7,
        hgt_competence=1.0,
    )
    predators = synthetic_species(
        'Synthetic predators',
        t6ss_activity=0.05,
        immunity_factor=0.4,
        hgt_competence=0.8,
        predation_threshold=1.0,
        base_size=3.0,
    )
    prey = synthetic_species(
        'Synthetic prey',
        immunity_factor=0.2,
        hgt_competence=0.7,
    )

    # Coordinates straddle signed 20 um hash boundaries while remaining inside the documented
    # strict contact disks; see DOCUMENTATION.md §4.12 and §10.2.
    kin.members = [
        synthetic_cell(kin, -20.1, -0.2, 9.5, resistant=True),
        synthetic_cell(kin, -19.8, 0.1, 2.0),
        synthetic_cell(kin, -60.0, -40.0, 7.0),
    ]
    attackers.members = [
        synthetic_cell(attackers, 19.9, 0.0, 8.0, resistant=True),
        synthetic_cell(attackers, 20.3, 0.2, 7.5),
    ]
    persisters.members = [
        synthetic_cell(persisters, 20.1, 0.1, 0.35, persister=True),
        synthetic_cell(persisters, 21.0, -0.2, 6.0),
    ]
    predators.members = [
        synthetic_cell(
            predators,
            100.0,
            0.0,
            0.1,
            base_size=3.0,
            interaction_radius=3.0,
        ),
    ]
    prey.members = [
        synthetic_cell(prey, 101.0, 0.0, 6.0),
        synthetic_cell(prey, -19.9, -0.1, 5.0),
    ]
    return [kin, attackers, persisters, predators, prey]


def resize_medium_fixture(species, tick):
    """Drive grow/shrink/regrow transitions without changing the shared object-order contract."""
    if tick == 10:
        # Crossing the initial 16-agent capacity forces geometric buffer growth; see future
        # DOCUMENTATION.md §10.2.
        for clone_index in range(12):
            owner = species[clone_index % len(species)]
            clone = copy.copy(owner.members[0])
            clone.species = owner
            clone.x = 300.0 + 30.0 * clone_index
            clone.y = -300.0
            clone.food_level = 4.0 + 0.1 * clone_index
            owner.members.append(clone)
    elif tick == 20:
        # A five-agent active slice leaves a deliberately large inactive capacity tail; see future
        # DOCUMENTATION.md §10.2.
        for owner in species:
            owner.members[:] = owner.members[:1]
    elif tick == 30:
        # Regrowing past the prior high-water mark exercises replacement storage and active slices;
        # see DOCUMENTATION.md §10.2.
        for clone_index in range(35):
            owner = species[clone_index % len(species)]
            clone = copy.copy(owner.members[0])
            clone.species = owner
            clone.x = -600.0 + 25.0 * clone_index
            clone.y = 350.0
            clone.food_level = 5.0 + 0.05 * clone_index
            clone.phage_resistant = clone_index % 4 == 0
            clone.metabolic_burden_factor = 1.0
            owner.members.append(clone)
    elif tick == 40:
        # The final shrink verifies that no stale high-capacity state re-enters active results;
        # see DOCUMENTATION.md §10.2.
        for species_index, owner in enumerate(species):
            owner.members[:] = owner.members[:2 + (species_index % 2)]


@unittest.skipUnless(NATIVE_PHYSICS_AVAILABLE, 'native physics kernel is not built')
class TestPhysicsMediumParity(unittest.TestCase):
    def test_multispecies_fifty_tick_parity_with_population_resizes(self):
        python_species = medium_fixture()
        cpp_species = clone_species_state(python_species)
        # The elevated synthetic HGT hazard makes conditional callback ordering observable over a
        # short deterministic run; the exponential hazard itself is documented in §4.18.
        parameters = replace(ContactParameters(), hgt_transfer_rate=8.0)
        python_orchestrator = PhysicsOrchestrator(
            backend='python',
            parameters=parameters,
        )
        cpp_orchestrator = PhysicsOrchestrator(
            backend='cpp',
            parameters=parameters,
        )
        branch_totals = {
            'contact_pairs': 0,
            't6ss_hits': 0,
            'hgt_transfers': 0,
            'share_events': 0,
            'predation_events': 0,
        }
        external_random_state = random.getstate()

        try:
            for tick in range(50):
                resize_medium_fixture(python_species, tick)
                resize_medium_fixture(cpp_species, tick)
                self.assertEqual(
                    [len(owner.members) for owner in cpp_species],
                    [len(owner.members) for owner in python_species],
                )

                # Restoring the same seeded stream before each backend proves equal conditional HGT
                # draw consumption as well as equal outcomes; see DOCUMENTATION.md §4.18.
                random.seed(MEDIUM_RANDOM_SEED + tick)
                python_summary = python_orchestrator.process(
                    python_species,
                    CONTACT_DT_SECONDS,
                )
                python_random_state = random.getstate()
                random.seed(MEDIUM_RANDOM_SEED + tick)
                cpp_summary = cpp_orchestrator.process(
                    cpp_species,
                    CONTACT_DT_SECONDS,
                )
                cpp_random_state = random.getstate()

                self.assertEqual(cpp_random_state, python_random_state, msg=f'tick={tick}')
                self.assertEqual(cpp_summary, python_summary, msg=f'tick={tick}')
                for field_name in branch_totals:
                    branch_totals[field_name] += getattr(python_summary, field_name)
                python_float_state, python_discrete_state = object_state(python_species)
                cpp_float_state, cpp_discrete_state = object_state(cpp_species)
                np.testing.assert_allclose(
                    cpp_float_state,
                    python_float_state,
                    rtol=0.0,
                    atol=FLOAT_ABSOLUTE_TOLERANCE,
                    err_msg=f'tick={tick}',
                )
                self.assertEqual(
                    cpp_discrete_state,
                    python_discrete_state,
                    msg=f'tick={tick}',
                )
                active_count = len(ordered_cells(python_species))
                np.testing.assert_array_equal(
                    cpp_orchestrator.buffers._agent_arrays['event_flags'][:active_count],
                    python_orchestrator.buffers._agent_arrays['event_flags'][:active_count],
                    err_msg=f'tick={tick}',
                )

            self.assertEqual(python_orchestrator.backend_used, 'python')
            self.assertEqual(cpp_orchestrator.backend_used, 'cpp')
            self.assertGreaterEqual(python_orchestrator.buffers.agent_capacity, 40)
            self.assertGreaterEqual(cpp_orchestrator.buffers.agent_capacity, 40)
            for field_name, total in branch_totals.items():
                self.assertGreater(total, 0, msg=f'{field_name} branch was not exercised')
        finally:
            random.setstate(external_random_state)


@unittest.skipUnless(NATIVE_PHYSICS_AVAILABLE, 'native physics kernel is not built')
class TestNativePhysicsABI(unittest.TestCase):
    def assert_rejected_without_writes(
        self,
        arguments,
        exception_type,
        message_pattern,
    ):
        snapshots = {
            name: arguments[name].copy()
            for name in MUTABLE_AGENT_ARGUMENTS
        }
        with self.assertRaisesRegex(exception_type, message_pattern):
            bacteria_c_kernel.process_interactions(**arguments)
        for name, expected in snapshots.items():
            np.testing.assert_equal(
                arguments[name],
                expected,
                err_msg=f'{name} changed after rejected native call',
            )

    def test_zero_agents_zero_species_and_single_agent(self):
        empty_arguments = valid_native_arguments(agent_count=0, species_count=0)
        self.assertEqual(
            bacteria_c_kernel.process_interactions(**empty_arguments),
            (0, 0, 0, 0, 0, 0),
        )

        single_arguments = valid_native_arguments(agent_count=1, species_count=1)
        mutable_before = {
            name: single_arguments[name].copy()
            for name in MUTABLE_AGENT_ARGUMENTS
        }
        self.assertEqual(
            bacteria_c_kernel.process_interactions(**single_arguments),
            (0, 0, 0, 0, 0, 0),
        )
        for name, expected in mutable_before.items():
            np.testing.assert_array_equal(single_arguments[name], expected)

    def test_signed_coordinate_conversion_boundaries_and_overflow(self):
        for sign in (-1.0, 1.0):
            with self.subTest(case='representable', sign=sign):
                arguments = valid_native_arguments(agent_count=1, species_count=1)
                # A 2^62 cell leaves room for both adjacent int64 keys while probing signed
                # floor conversion far beyond ordinary worlds; see DOCUMENTATION.md §10.2.
                arguments['x'][0] = (
                    sign * float(1 << 62) * arguments['hash_cell_size']
                )
                self.assertEqual(
                    bacteria_c_kernel.process_interactions(**arguments),
                    (0, 0, 0, 0, 0, 0),
                )

            with self.subTest(case='int64 overflow', sign=sign):
                arguments = valid_native_arguments(agent_count=1, species_count=1)
                # A signed 2^63 cell cannot preserve the required neighboring-key range and must
                # fail before resolution; see DOCUMENTATION.md §4.12 and §10.2.
                arguments['x'][0] = (
                    sign * float(1 << 63) * arguments['hash_cell_size']
                )
                self.assert_rejected_without_writes(
                    arguments,
                    OverflowError,
                    'cell coordinate.*adjacent int64',
                )

        arguments = valid_native_arguments(agent_count=1, species_count=1)
        arguments['x'][0] = np.finfo(np.float64).max
        self.assert_rejected_without_writes(
            arguments,
            OverflowError,
            'cell coordinate.*adjacent int64',
        )

    def test_every_mutable_buffer_must_be_writable(self):
        for argument_name in MUTABLE_AGENT_ARGUMENTS:
            with self.subTest(argument=argument_name):
                arguments = valid_native_arguments()
                arguments[argument_name].setflags(write=False)
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    f'{argument_name}.*writable',
                )

    def test_every_array_rejects_wrong_dimensionality_and_dtype(self):
        array_names = (
            FLOAT_AGENT_ARGUMENTS
            + ('species_id',)
            + UINT8_AGENT_ARGUMENTS
            + FLOAT_SPECIES_ARGUMENTS
        )
        for argument_name in array_names:
            with self.subTest(argument=argument_name, defect='dimensions'):
                arguments = valid_native_arguments()
                arguments[argument_name] = arguments[argument_name].reshape((1, -1))
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    f'{argument_name}.*shape',
                )

            with self.subTest(argument=argument_name, defect='dtype'):
                arguments = valid_native_arguments()
                if argument_name in FLOAT_AGENT_ARGUMENTS + FLOAT_SPECIES_ARGUMENTS:
                    invalid_dtype = np.float32
                elif argument_name == 'species_id':
                    invalid_dtype = np.int64
                else:
                    invalid_dtype = np.int8
                arguments[argument_name] = arguments[argument_name].astype(invalid_dtype)
                self.assert_rejected_without_writes(
                    arguments,
                    TypeError,
                    f'{argument_name}.*invalid dtype',
                )

    def test_all_agent_and_species_length_mismatch_families(self):
        agent_dependents = (
            'y',
            'interaction_radius',
            'base_radius',
            'food_level',
            'food_capacity',
            'metabolic_burden_factor',
            'species_id',
            'is_persister',
            'phage_resistant',
            'event_flags',
        )
        for argument_name in agent_dependents:
            with self.subTest(extent='N', argument=argument_name):
                arguments = valid_native_arguments()
                arguments[argument_name] = arguments[argument_name][:-1]
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    rf'{argument_name}.*\(N,\)',
                )

        arguments = valid_native_arguments()
        arguments['x'] = arguments['x'][:-1]
        self.assert_rejected_without_writes(arguments, ValueError, r'y.*\(N,\)')

        species_dependents = (
            'immunity_factor',
            'hgt_competence',
            'share_threshold',
            'predation_threshold',
        )
        for argument_name in species_dependents:
            with self.subTest(extent='S', argument=argument_name):
                arguments = valid_native_arguments()
                arguments[argument_name] = arguments[argument_name][:-1]
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    rf'{argument_name}.*\(S,\)',
                )

        arguments = valid_native_arguments()
        arguments['t6ss_activity'] = arguments['t6ss_activity'][:-1]
        self.assert_rejected_without_writes(
            arguments,
            ValueError,
            r'immunity_factor.*\(S,\)',
        )

    def test_nonfinite_values_rejected_for_every_array_and_scalar_family(self):
        for argument_name in FLOAT_AGENT_ARGUMENTS + FLOAT_SPECIES_ARGUMENTS:
            with self.subTest(kind='array', argument=argument_name):
                arguments = valid_native_arguments()
                arguments[argument_name][0] = math.nan
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    f'{argument_name}.*finite',
                )

        scalar_names = (
            'dt',
            'hash_cell_size',
            'axis_cutoff',
            'contact_margin',
            'hgt_transfer_rate',
            'hgt_metabolic_burden',
            't6ss_damage_scaling',
            't6ss_max_immunity_defense',
            't6ss_immunity_size_epsilon',
            'persister_resistance_factor',
            'phagocytosis_size_ratio',
            'cell_wet_mass',
            'cell_dry_mass_fraction',
            'labile_necromass_fraction',
            'cytosol_recovery_fraction',
            'consumed_prey_sentinel',
        )
        for argument_name in scalar_names:
            with self.subTest(kind='scalar', argument=argument_name):
                arguments = valid_native_arguments()
                arguments[argument_name] = math.nan
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    f'{argument_name}.*finite',
                )

    def test_out_of_range_agent_status_trait_and_parameter_values(self):
        nonnegative_agent_arrays = (
            'interaction_radius',
            'base_radius',
            'food_capacity',
            'metabolic_burden_factor',
        )
        for argument_name in nonnegative_agent_arrays:
            with self.subTest(kind='nonnegative array', argument=argument_name):
                arguments = valid_native_arguments()
                arguments[argument_name][0] = -0.1
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    f'{argument_name}.*nonnegative',
                )

        for argument_name in ('is_persister', 'phage_resistant'):
            with self.subTest(kind='status', argument=argument_name):
                arguments = valid_native_arguments()
                arguments[argument_name][0] = 2
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    f'{argument_name}.*0 or 1',
                )

        for species_id in (-1, 2):
            with self.subTest(kind='species id', value=species_id):
                arguments = valid_native_arguments()
                arguments['species_id'][0] = species_id
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    'species_id.*species arrays',
                )

        arguments = valid_native_arguments()
        arguments['event_flags'][0] = 1 << 6
        self.assert_rejected_without_writes(
            arguments,
            ValueError,
            'event_flags.*unknown event bits',
        )

        for argument_name in FLOAT_SPECIES_ARGUMENTS:
            for invalid_value in (-0.1, 1.1):
                with self.subTest(
                    kind='species trait',
                    argument=argument_name,
                    value=invalid_value,
                ):
                    arguments = valid_native_arguments()
                    arguments[argument_name][0] = invalid_value
                    self.assert_rejected_without_writes(
                        arguments,
                        ValueError,
                        f'{argument_name}.*\[0, 1\]',
                    )

        nonnegative_scalars = (
            'dt',
            'axis_cutoff',
            'contact_margin',
            'hgt_transfer_rate',
            'hgt_metabolic_burden',
            't6ss_damage_scaling',
            'phagocytosis_size_ratio',
            'cell_wet_mass',
        )
        for argument_name in nonnegative_scalars:
            with self.subTest(kind='nonnegative scalar', argument=argument_name):
                arguments = valid_native_arguments()
                arguments[argument_name] = -0.1
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    f'{argument_name}.*nonnegative',
                )

        unit_interval_scalars = (
            't6ss_max_immunity_defense',
            'persister_resistance_factor',
            'cell_dry_mass_fraction',
            'labile_necromass_fraction',
            'cytosol_recovery_fraction',
        )
        for argument_name in unit_interval_scalars:
            with self.subTest(kind='unit scalar', argument=argument_name):
                arguments = valid_native_arguments()
                arguments[argument_name] = 1.1
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    f'{argument_name}.*\[0, 1\]',
                )

        positive_scalars = ('hash_cell_size', 't6ss_immunity_size_epsilon')
        for argument_name in positive_scalars:
            with self.subTest(kind='positive scalar', argument=argument_name):
                arguments = valid_native_arguments()
                arguments[argument_name] = 0.0
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    f'{argument_name}.*positive',
                )

        arguments = valid_native_arguments()
        arguments['consumed_prey_sentinel'] = 0.0
        self.assert_rejected_without_writes(
            arguments,
            ValueError,
            'consumed_prey_sentinel.*negative',
        )

    def test_invalid_random_draws_reject_without_partial_writes(self):
        for invalid_draw in (-0.1, 1.0, math.nan):
            with self.subTest(draw=invalid_draw):
                arguments = valid_native_arguments()
                arguments['random_draw'] = lambda value=invalid_draw: value
                self.assert_rejected_without_writes(
                    arguments,
                    ValueError,
                    r'random_draw.*\[0, 1\)',
                )

        arguments = valid_native_arguments(agent_count=3, species_count=1)
        delayed_draws = iter((0.0, math.nan))
        arguments['random_draw'] = lambda: next(delayed_draws)
        self.assert_rejected_without_writes(
            arguments,
            ValueError,
            r'random_draw.*\[0, 1\)',
        )


class TestPhysicsBufferOwnership(unittest.TestCase):
    @staticmethod
    def _isolated_cells(cell_count):
        species = synthetic_species('Synthetic ownership cells')
        cells = []
        for cell_index in range(cell_count):
            # A 100 um stride keeps cells outside neighboring contact disks so any write-back
            # mismatch is attributable to buffer slicing, not biology; see DOCUMENTATION.md §4.12.
            cells.append(synthetic_cell(
                species,
                100.0 * cell_index,
                25.0,
                2.0 + 0.05 * cell_index,
                resistant=cell_index % 3 == 0,
                persister=cell_index % 4 == 0,
            ))
        return species, cells

    def test_grow_shrink_regrow_writeback_tail_and_flag_isolation(self):
        species, cells = self._isolated_cells(40)
        orchestrator = PhysicsOrchestrator(backend='python')

        species.members = cells[:3]
        # The first pair starts in contact solely to seed event bits that the next flatten must
        # clear; see DOCUMENTATION.md §4.12 and §10.2.
        species.members[0].x = 0.0
        species.members[1].x = 0.5
        first_summary = orchestrator.process([species], CONTACT_DT_SECONDS)
        self.assertGreater(first_summary.contact_pairs, 0)
        self.assertTrue(
            orchestrator.buffers._agent_arrays['event_flags'][0] & EVENT_CONTACT
        )

        for active_count in (24, 5):
            species.members = cells[:active_count]
            for cell_index, cell in enumerate(species.members):
                cell.x = 100.0 * cell_index
                cell.food_level = 3.0 + 0.01 * cell_index
                cell.metabolic_burden_factor = 1.0 + 0.001 * cell_index
            summary = orchestrator.process([species], CONTACT_DT_SECONDS)
            self.assertEqual(summary.contact_pairs, 0)
            np.testing.assert_array_equal(
                orchestrator.buffers._agent_arrays['event_flags'][:active_count],
                np.zeros(active_count, dtype=np.uint8),
            )

        capacity = orchestrator.buffers.agent_capacity
        self.assertGreater(capacity, 5)
        # Poison every inactive field that can affect write-back or event interpretation; active
        # flattening must overwrite it before reuse, see DOCUMENTATION.md §10.2.
        orchestrator.buffers._agent_arrays['food_level'][5:capacity] = -999.0
        orchestrator.buffers._agent_arrays['metabolic_burden_factor'][5:capacity] = 99.0
        orchestrator.buffers._agent_arrays['phage_resistant'][5:capacity] = 1
        orchestrator.buffers._agent_arrays['event_flags'][5:capacity] = 0xFF

        for active_count in (20, 40, 7):
            species.members = cells[:active_count]
            expected_food = []
            expected_burden = []
            expected_resistance = []
            for cell_index, cell in enumerate(species.members):
                cell.x = 100.0 * cell_index
                cell.food_level = 4.0 + 0.02 * cell_index
                cell.metabolic_burden_factor = 1.0 + 0.002 * cell_index
                cell.phage_resistant = cell_index % 5 == 0
                expected_food.append(cell.food_level)
                expected_burden.append(cell.metabolic_burden_factor)
                expected_resistance.append(cell.phage_resistant)

            summary = orchestrator.process([species], CONTACT_DT_SECONDS)
            self.assertEqual(summary.contact_pairs, 0)
            self.assertEqual(len(orchestrator.buffers._agent_objects), active_count)
            for cell_index, cell in enumerate(species.members):
                self.assertIs(orchestrator.buffers._agent_objects[cell_index], cell)
                self.assertEqual(cell.food_level, expected_food[cell_index])
                self.assertEqual(
                    cell.metabolic_burden_factor,
                    expected_burden[cell_index],
                )
                self.assertEqual(
                    cell.phage_resistant,
                    expected_resistance[cell_index],
                )
            np.testing.assert_array_equal(
                orchestrator.buffers._agent_arrays['event_flags'][:active_count],
                np.zeros(active_count, dtype=np.uint8),
            )

    def test_frame_views_and_temporary_native_callback_are_collectible(self):
        species, cells = self._isolated_cells(4)
        species.members = cells
        buffers = PhysicsBuffers()
        frame = buffers.flatten([species])
        frame_reference = weakref.ref(frame)
        view_reference = weakref.ref(frame.agents.food_level)
        del frame
        gc.collect()
        self.assertIsNone(frame_reference())
        self.assertIsNone(view_reference())

        if not NATIVE_PHYSICS_AVAILABLE:
            self.skipTest('native physics kernel is not built')

        class TemporaryDraw:
            def __call__(self):
                return 0.5

        callback = TemporaryDraw()
        callback_reference = weakref.ref(callback)
        arguments = valid_native_arguments()
        arguments['random_draw'] = callback
        result = bacteria_c_kernel.process_interactions(**arguments)
        self.assertEqual(len(result), 6)
        del result
        del arguments
        del callback
        gc.collect()
        self.assertIsNone(callback_reference())

    def test_python_retained_allocation_is_bounded_after_warmup(self):
        species, cells = self._isolated_cells(48)
        orchestrator = PhysicsOrchestrator(backend='python')
        population_pattern = (8, 32, 16, 48, 12)

        for iteration in range(25):
            species.members = cells[:population_pattern[iteration % len(population_pattern)]]
            orchestrator.process([species], CONTACT_DT_SECONDS)

        gc.collect()
        tracemalloc.start()
        try:
            for iteration in range(25):
                species.members = cells[:population_pattern[iteration % len(population_pattern)]]
                orchestrator.process([species], CONTACT_DT_SECONDS)
            gc.collect()
            baseline_bytes = tracemalloc.get_traced_memory()[0]

            samples = []
            for iteration in range(200):
                species.members = cells[:population_pattern[iteration % len(population_pattern)]]
                orchestrator.process([species], CONTACT_DT_SECONDS)
                if (iteration + 1) % 20 == 0:
                    gc.collect()
                    samples.append(tracemalloc.get_traced_memory()[0])

            retained_growth = samples[-1] - baseline_bytes
            self.assertLessEqual(
                retained_growth,
                ALLOCATION_GROWTH_TOLERANCE_BYTES,
                msg=f'retained Python allocation grew by {retained_growth} bytes',
            )
            monotonically_increasing = all(
                later > earlier
                for earlier, later in zip(samples, samples[1:])
            )
            self.assertFalse(
                monotonically_increasing
                and retained_growth > ALLOCATION_GROWTH_TOLERANCE_BYTES // 4,
                msg=f'suspicious retained-allocation samples: {samples}',
            )
        finally:
            tracemalloc.stop()


@unittest.skipUnless(NATIVE_PHYSICS_AVAILABLE, 'native physics kernel is not built')
class TestNativePhysicsStress(unittest.TestCase):
    """Valgrind workload target; Python assertions alone do not prove native leak freedom."""

    def test_repeated_preallocated_calls(self):
        iteration_count = int(os.environ.get(
            'BACSANDBOX_NATIVE_PHYSICS_STRESS_ITERATIONS',
            str(NATIVE_PHYSICS_STRESS_DEFAULT_ITERATIONS),
        ))
        self.assertGreater(iteration_count, 0)
        arguments = valid_native_arguments(agent_count=32, species_count=4)
        # Synthetic positive activity drives mutable T6SS state without invoking unrelated
        # simulation systems; see DOCUMENTATION.md §4.12.
        arguments['t6ss_activity'][:] = (0.1, 0.2, 0.3, 0.4)
        arguments['hgt_competence'].fill(0.0)
        array_addresses = {
            name: arguments[name].ctypes.data
            for name in (
                FLOAT_AGENT_ARGUMENTS
                + ('species_id',)
                + UINT8_AGENT_ARGUMENTS
                + FLOAT_SPECIES_ARGUMENTS
            )
        }

        for iteration in range(iteration_count):
            # Resetting mutable state isolates repeated native ownership while a small signed shift
            # changes hash contents without reallocating Python inputs; see DOCUMENTATION.md §10.2.
            arguments['x'][0] = -0.25 * (iteration % 5)
            arguments['food_level'].fill(5.0 + 0.01 * (iteration % 3))
            arguments['metabolic_burden_factor'].fill(1.0)
            arguments['phage_resistant'].fill(0)
            arguments['event_flags'].fill(0)
            summary = bacteria_c_kernel.process_interactions(**arguments)

            self.assertIsInstance(summary, tuple)
            self.assertEqual(len(summary), 6)
            self.assertTrue(all(isinstance(value, int) for value in summary))
            self.assertTrue(all(value >= 0 for value in summary))
            self.assertTrue(np.all(np.isfinite(arguments['food_level'])))
            self.assertTrue(np.all(np.isfinite(arguments['metabolic_burden_factor'])))

        for name, address in array_addresses.items():
            self.assertEqual(arguments[name].ctypes.data, address, msg=name)


def scale_fixture():
    """Build five named, explicitly synthetic 200-cell contact populations."""
    # These branch-coverage traits are synthetic and do not represent measurements for the named
    # strains; they span T6SS, immunity, HGT, and persistence behavior documented in §4.12/§4.18.
    trait_rows = (
        (0.20, 0.15, 1.00, 1.00),
        (0.35, 0.30, 0.90, 1.05),
        (0.50, 0.45, 0.80, 0.95),
        (0.65, 0.60, 0.70, 1.10),
        (0.80, 0.75, 0.60, 1.00),
    )
    species = []
    for name, (t6ss, immunity, competence, base_size) in zip(
        SCALE_SPECIES_NAMES,
        trait_rows,
    ):
        species.append(synthetic_species(
            name,
            t6ss_activity=t6ss,
            immunity_factor=immunity,
            hgt_competence=competence,
            share_threshold=1.0,
            predation_threshold=0.0,
            base_size=base_size,
        ))

    fixture_random = random.Random(SCALE_RANDOM_SEED)
    for species_index, owner in enumerate(species):
        for cohort_index in range(SCALE_AGENTS_PER_SPECIES):
            # Eighty-micrometer cohort spacing isolates five-cell contact groups, keeping the
            # 1,000x500 fixture bounded while still dispatching real contacts; see DOCUMENTATION.md §4.12 and §10.2.
            cohort_x = 80.0 * (cohort_index % 20)
            cohort_y = 80.0 * (cohort_index // 20)
            # Sub-micrometer offsets keep all five species in strict disk contact without exact
            # coordinate ties; see DOCUMENTATION.md §4.12.
            x = cohort_x + 0.32 * species_index + fixture_random.uniform(-0.01, 0.01)
            y = cohort_y + 0.07 * species_index + fixture_random.uniform(-0.01, 0.01)
            # Seeded donors only in the first population make cross-species plasmid propagation
            # measurable while remaining a synthetic HGT fixture; see DOCUMENTATION.md §4.18.
            resistant = species_index == 0 and fixture_random.random() < 0.25
            persister = (cohort_index + species_index) % 11 == 0
            owner.members.append(synthetic_cell(
                owner,
                x,
                y,
                owner.food_capacity * SCALE_STEADY_FULLNESS,
                base_size=owner.base_size,
                interaction_radius=owner.base_size,
                resistant=resistant,
                persister=persister,
            ))
    return species


def prepare_scale_tick(species, tick):
    """Reset unrelated energy drift and select deterministic T6SS probes."""
    for species_index, owner in enumerate(species):
        for cell in owner.members:
            cell.food_level = cell.food_capacity * SCALE_STEADY_FULLNESS
        if tick % SCALE_TURNOVER_INTERVAL == 0:
            # Coprime index strides rotate probes through each population without stochastic
            # selection, preserving backend identity order; see DOCUMENTATION.md §4.12 and §10.2.
            victim_index = (
                7 * (tick // SCALE_TURNOVER_INTERVAL) + 13 * species_index
            ) % len(owner.members)
            owner.members[victim_index].food_level = SCALE_LETHAL_PROBE_FOOD


def apply_scale_turnover(species, event_flags, tick):
    """Remove T6SS-lethal cells and add deterministic replacement births."""
    lethal_by_species = []
    start = 0
    for species_index, owner in enumerate(species):
        stop = start + len(owner.members)
        owner_flags = event_flags[start:stop]
        lethal_indices = [
            member_index
            for member_index, flags in enumerate(owner_flags)
            if flags & EVENT_T6SS_LETHAL
        ]
        lethal_by_species.append(len(lethal_indices))
        lethal_index_set = set(lethal_indices)
        dead_cells = [owner.members[index] for index in lethal_indices]
        owner.members[:] = [
            cell
            for member_index, cell in enumerate(owner.members)
            if member_index not in lethal_index_set
        ]
        owner.total_deaths += len(dead_cells)

        for death_index, dead_cell in enumerate(dead_cells):
            # Replacement births preserve the dead cell's cohort while inheriting resistance and
            # burden from a deterministic survivor, matching §4.18 without invoking other systems.
            parent_index = (tick + species_index + death_index) % len(owner.members)
            parent = owner.members[parent_index]
            child = copy.copy(parent)
            child.species = owner
            child.x = dead_cell.x
            child.y = dead_cell.y
            child.food_level = child.food_capacity * SCALE_STEADY_FULLNESS
            child.is_persister = False
            owner.members.append(child)
        owner.total_births += len(dead_cells)
        start = stop
    return lethal_by_species


def assert_runtime_species_equal(test_case, expected, actual, tick):
    test_case.assertEqual(
        [len(owner.members) for owner in actual],
        [len(owner.members) for owner in expected],
        msg=f'population tick={tick}',
    )
    test_case.assertEqual(
        [owner.total_births for owner in actual],
        [owner.total_births for owner in expected],
        msg=f'birth totals tick={tick}',
    )
    test_case.assertEqual(
        [owner.total_deaths for owner in actual],
        [owner.total_deaths for owner in expected],
        msg=f'death totals tick={tick}',
    )
    expected_float, expected_discrete = object_state(expected)
    actual_float, actual_discrete = object_state(actual)
    np.testing.assert_allclose(
        actual_float,
        expected_float,
        rtol=0.0,
        atol=FLOAT_ABSOLUTE_TOLERANCE,
        err_msg=f'aligned floating state tick={tick}',
    )
    test_case.assertEqual(actual_discrete, expected_discrete, msg=f'tick={tick}')


@unittest.skipUnless(RUN_PHYSICS_SCALE, 'set BACSANDBOX_RUN_PHYSICS_SCALE=1 to run')
@unittest.skipUnless(NATIVE_PHYSICS_AVAILABLE, 'native physics kernel is not built')
class TestPhysicsScale(unittest.TestCase):
    def test_thousand_agents_five_hundred_ticks(self):
        # This contact-only harness disables movement, metabolism, and environmental stochasticity;
        # deterministic turnover isolates §4.12/§4.18 dispatch parity at the required scale.
        python_species = scale_fixture()
        cpp_species = clone_species_state(python_species)
        self.assertEqual(
            [len(owner.members) for owner in python_species],
            [SCALE_AGENTS_PER_SPECIES] * len(SCALE_SPECIES_NAMES),
        )
        self.assertEqual(len(ordered_cells(python_species)), 1000)

        # A synthetic elevated rate produces enough transfers for deterministic cumulative checks;
        # probability still uses the production exponential hazard in DOCUMENTATION.md §4.18.
        parameters = replace(ContactParameters(), hgt_transfer_rate=6.0)
        python_orchestrator = PhysicsOrchestrator('python', parameters=parameters)
        cpp_orchestrator = PhysicsOrchestrator('cpp', parameters=parameters)
        python_duration = 0.0
        cpp_duration = 0.0
        python_t6ss_lethal = 0
        cpp_t6ss_lethal = 0
        python_hgt_transfers = 0
        cpp_hgt_transfers = 0
        python_lethal_by_species = [0] * len(SCALE_SPECIES_NAMES)
        cpp_lethal_by_species = [0] * len(SCALE_SPECIES_NAMES)
        external_random_state = random.getstate()
        scale_start = time.perf_counter()

        try:
            for tick in range(SCALE_TICK_COUNT):
                prepare_scale_tick(python_species, tick)
                prepare_scale_tick(cpp_species, tick)

                random.seed(SCALE_RANDOM_SEED + tick)
                start = time.perf_counter()
                python_summary = python_orchestrator.process(
                    python_species,
                    CONTACT_DT_SECONDS,
                )
                python_duration += time.perf_counter() - start
                python_random_state = random.getstate()

                random.seed(SCALE_RANDOM_SEED + tick)
                start = time.perf_counter()
                cpp_summary = cpp_orchestrator.process(
                    cpp_species,
                    CONTACT_DT_SECONDS,
                )
                cpp_duration += time.perf_counter() - start
                cpp_random_state = random.getstate()

                self.assertEqual(cpp_random_state, python_random_state, msg=f'tick={tick}')
                self.assertEqual(cpp_summary, python_summary, msg=f'tick={tick}')
                active_count = len(ordered_cells(python_species))
                python_flags = python_orchestrator.buffers._agent_arrays[
                    'event_flags'
                ][:active_count].copy()
                cpp_flags = cpp_orchestrator.buffers._agent_arrays[
                    'event_flags'
                ][:active_count].copy()
                np.testing.assert_array_equal(cpp_flags, python_flags, err_msg=f'tick={tick}')

                python_t6ss_lethal += python_summary.t6ss_deaths
                cpp_t6ss_lethal += cpp_summary.t6ss_deaths
                python_hgt_transfers += python_summary.hgt_transfers
                cpp_hgt_transfers += cpp_summary.hgt_transfers
                python_tick_lethal = apply_scale_turnover(
                    python_species,
                    python_flags,
                    tick,
                )
                cpp_tick_lethal = apply_scale_turnover(cpp_species, cpp_flags, tick)
                python_lethal_by_species = [
                    total + count
                    for total, count in zip(
                        python_lethal_by_species,
                        python_tick_lethal,
                    )
                ]
                cpp_lethal_by_species = [
                    total + count
                    for total, count in zip(cpp_lethal_by_species, cpp_tick_lethal)
                ]

                self.assertEqual(cpp_tick_lethal, python_tick_lethal, msg=f'tick={tick}')
                self.assertEqual(cpp_t6ss_lethal, python_t6ss_lethal, msg=f'tick={tick}')
                self.assertEqual(cpp_hgt_transfers, python_hgt_transfers, msg=f'tick={tick}')
                self.assertEqual(
                    cpp_lethal_by_species,
                    python_lethal_by_species,
                    msg=f'tick={tick}',
                )
                assert_runtime_species_equal(
                    self,
                    python_species,
                    cpp_species,
                    tick,
                )
        finally:
            random.setstate(external_random_state)

        total_duration = time.perf_counter() - scale_start
        expected_turnovers = SCALE_TICK_COUNT // SCALE_TURNOVER_INTERVAL
        self.assertEqual(
            [len(owner.members) for owner in python_species],
            [SCALE_AGENTS_PER_SPECIES] * len(SCALE_SPECIES_NAMES),
        )
        self.assertEqual(
            [owner.total_births for owner in python_species],
            [expected_turnovers] * len(SCALE_SPECIES_NAMES),
        )
        self.assertEqual(
            [owner.total_deaths for owner in python_species],
            [expected_turnovers] * len(SCALE_SPECIES_NAMES),
        )
        self.assertEqual(
            python_lethal_by_species,
            [expected_turnovers] * len(SCALE_SPECIES_NAMES),
        )
        self.assertEqual(python_t6ss_lethal, sum(python_lethal_by_species))
        self.assertGreater(python_hgt_transfers, 0)
        print(
            'Physics Scale Outcomes: '
            f'populations={[len(owner.members) for owner in python_species]}; '
            f'births={[owner.total_births for owner in python_species]}; '
            f'deaths={[owner.total_deaths for owner in python_species]}; '
            f'T6SS lethal={python_t6ss_lethal}; HGT transfers={python_hgt_transfers}; '
            f'Python={python_duration:.6f}s; C++={cpp_duration:.6f}s; '
            f'total={total_duration:.6f}s'
        )


def benchmark_fixture():
    """Build a contact-only 1,000-agent mix of dense cohorts and isolated cells."""
    species = [
        synthetic_species(f'Synthetic benchmark species {index}')
        for index in range(5)
    ]
    for species_index, owner in enumerate(species):
        for member_index in range(SCALE_AGENTS_PER_SPECIES):
            if member_index < 100:
                # Twenty 25-cell cohorts create bounded dense contact work; see DOCUMENTATION.md §4.12 and §10.2.
                group_index = member_index // 5
                within_species_slot = member_index % 5
                x = 100.0 * group_index + 0.15 * species_index
                y = 0.04 * within_species_slot
            else:
                # The other 500 cells are separated by 100 um to exercise sparse hashing without
                # accepted contacts; see DOCUMENTATION.md §4.12 and §10.2.
                sparse_index = species_index * 100 + (member_index - 100)
                x = 5000.0 + 100.0 * sparse_index
                y = 5000.0
            owner.members.append(synthetic_cell(owner, x, y, 5.0))
    return species


def timed_contact_workload(backend, baseline, seed):
    species = clone_species_state(baseline)
    orchestrator = PhysicsOrchestrator(backend)
    random.seed(seed)
    summaries = []
    start = time.perf_counter()
    for _ in range(BENCHMARK_TICK_COUNT):
        summaries.append(orchestrator.process(species, CONTACT_DT_SECONDS))
    duration = time.perf_counter() - start
    float_state, discrete_state = object_state(species)
    active_count = len(ordered_cells(species))
    flags = orchestrator.buffers._agent_arrays['event_flags'][:active_count].copy()
    return (
        duration,
        summaries,
        float_state,
        discrete_state,
        flags,
        random.getstate(),
    )


@unittest.skipUnless(
    RUN_PHYSICS_BENCHMARK,
    'set BACSANDBOX_RUN_PHYSICS_BENCHMARK=1 to run',
)
@unittest.skipUnless(NATIVE_PHYSICS_AVAILABLE, 'native physics kernel is not built')
class TestPhysicsBenchmark(unittest.TestCase):
    """Contact-kernel benchmark excluding movement, metabolism, and environmental work."""

    def test_contact_kernel_speedup(self):
        baseline = benchmark_fixture()
        self.assertEqual(len(ordered_cells(baseline)), 1000)
        external_random_state = random.getstate()
        python_durations = []
        cpp_durations = []

        try:
            timed_contact_workload('python', baseline, SCALE_RANDOM_SEED)
            timed_contact_workload('cpp', baseline, SCALE_RANDOM_SEED)

            for sample_index in range(BENCHMARK_REPEAT_COUNT):
                seed = SCALE_RANDOM_SEED + sample_index
                python_result = timed_contact_workload('python', baseline, seed)
                cpp_result = timed_contact_workload('cpp', baseline, seed)
                python_durations.append(python_result[0])
                cpp_durations.append(cpp_result[0])

                self.assertEqual(cpp_result[1], python_result[1])
                np.testing.assert_allclose(
                    cpp_result[2],
                    python_result[2],
                    rtol=0.0,
                    atol=FLOAT_ABSOLUTE_TOLERANCE,
                )
                self.assertEqual(cpp_result[3], python_result[3])
                np.testing.assert_array_equal(cpp_result[4], python_result[4])
                self.assertEqual(cpp_result[5], python_result[5])
        finally:
            random.setstate(external_random_state)

        python_median = statistics.median(python_durations)
        cpp_median = statistics.median(cpp_durations)
        speedup = python_median / cpp_median
        print(
            f'Physics Kernel Speedup Factor: {speedup:.3f}x '
            f'(Python median: {python_median:.6f}s; C++ median: {cpp_median:.6f}s)'
        )


if __name__ == '__main__':
    unittest.main()