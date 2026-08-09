import importlib
import math
import random
from dataclasses import dataclass

import numpy as np

from .constants import (
    BACTERIA_CONTACT_MARGIN_UM,
    BACTERIA_INTERACTION_AXIS_CUTOFF_UM,
    BACTERIA_INTERACTION_HASH_CELL_SIZE_UM,
    CELL_DRY_MASS_FRACTION,
    CELL_WET_MASS_PG_PER_UM3,
    CONSUMED_PREY_SENTINEL_FOOD_LEVEL,
    CYTOSOL_RECOVERY_FRACTION,
    HGT_METABOLIC_BURDEN,
    HGT_TRANSFER_RATE_PER_SEC,
    LABILE_NECROMASS_FRACTION,
    PERSISTER_RESISTANCE_FACTOR,
    PHAGOCYTOSIS_SIZE_RATIO_THRESHOLD,
    T6SS_DAMAGE_SCALING_FACTOR,
    T6SS_IMMUNITY_SIZE_EPSILON_UM,
    T6SS_MAX_IMMUNITY_DEFENSE,
)


_BACKENDS = frozenset(('auto', 'python', 'cpp'))
_NATIVE_MODULE_NAME = 'bacteria_c_kernel'
# The three signed offsets cover the source bucket and every adjacent bucket; see DOCUMENTATION.md §4.12.
_NEIGHBOR_OFFSETS = (-1, 0, 1)
# Geometric doubling starts with a small allocation to amortize resize costs without retaining a large empty arena;
# see DOCUMENTATION.md §10.2.
_MIN_BUFFER_CAPACITY = 16

EVENT_CONTACT = 1 << 0
EVENT_T6SS_HIT = 1 << 1
EVENT_T6SS_LETHAL = 1 << 2
EVENT_HGT_TRANSFER = 1 << 3
EVENT_SHARE = 1 << 4
EVENT_PREDATION = 1 << 5


class _NativePhysicsUnavailable(ImportError):
    pass


@dataclass(frozen=True)
class AgentArrays:
    x: np.ndarray
    y: np.ndarray
    interaction_radius: np.ndarray
    base_radius: np.ndarray
    food_level: np.ndarray
    food_capacity: np.ndarray
    metabolic_burden_factor: np.ndarray
    species_id: np.ndarray
    is_persister: np.ndarray
    phage_resistant: np.ndarray
    event_flags: np.ndarray


@dataclass(frozen=True)
class SpeciesArrays:
    t6ss_activity: np.ndarray
    immunity_factor: np.ndarray
    hgt_competence: np.ndarray
    share_threshold: np.ndarray
    predation_threshold: np.ndarray


@dataclass(frozen=True)
class InteractionFrame:
    agents: AgentArrays
    species: SpeciesArrays
    agent_objects: list


@dataclass(frozen=True)
class InteractionSummary:
    contact_pairs: int = 0
    t6ss_hits: int = 0
    t6ss_deaths: int = 0
    hgt_transfers: int = 0
    share_events: int = 0
    predation_events: int = 0

    @classmethod
    def from_native(cls, result):
        if not isinstance(result, tuple) or len(result) != 6:
            raise TypeError('bacteria_c_kernel.process_interactions must return a six-integer tuple')
        values = []
        for value in result:
            if not isinstance(value, (int, np.integer)):
                raise TypeError('native interaction summary values must be integers')
            if value < 0:
                raise ValueError('native interaction summary values must be non-negative')
            values.append(int(value))
        return cls(*values)


@dataclass(frozen=True)
class ContactParameters:
    hash_cell_size: float = BACTERIA_INTERACTION_HASH_CELL_SIZE_UM
    axis_cutoff: float = BACTERIA_INTERACTION_AXIS_CUTOFF_UM
    contact_margin: float = BACTERIA_CONTACT_MARGIN_UM
    hgt_transfer_rate: float = HGT_TRANSFER_RATE_PER_SEC
    hgt_metabolic_burden: float = HGT_METABOLIC_BURDEN
    t6ss_damage_scaling: float = T6SS_DAMAGE_SCALING_FACTOR
    t6ss_max_immunity_defense: float = T6SS_MAX_IMMUNITY_DEFENSE
    t6ss_immunity_size_epsilon: float = T6SS_IMMUNITY_SIZE_EPSILON_UM
    persister_resistance_factor: float = PERSISTER_RESISTANCE_FACTOR
    phagocytosis_size_ratio: float = PHAGOCYTOSIS_SIZE_RATIO_THRESHOLD
    cell_wet_mass: float = CELL_WET_MASS_PG_PER_UM3
    cell_dry_mass_fraction: float = CELL_DRY_MASS_FRACTION
    labile_necromass_fraction: float = LABILE_NECROMASS_FRACTION
    cytosol_recovery_fraction: float = CYTOSOL_RECOVERY_FRACTION
    consumed_prey_sentinel: float = CONSUMED_PREY_SENTINEL_FOOD_LEVEL


class PhysicsBuffers:
    """Reusable structure-of-arrays storage for the contact boundary."""

    _FLOAT_AGENT_FIELDS = (
        'x',
        'y',
        'interaction_radius',
        'base_radius',
        'food_level',
        'food_capacity',
        'metabolic_burden_factor',
    )
    _UINT8_AGENT_FIELDS = ('is_persister', 'phage_resistant', 'event_flags')
    _FLOAT_SPECIES_FIELDS = (
        't6ss_activity',
        'immunity_factor',
        'hgt_competence',
        'share_threshold',
        'predation_threshold',
    )

    def __init__(self):
        self._agent_capacity = 0
        self._species_capacity = 0
        self._agent_arrays = {
            name: np.empty(0, dtype=np.float64)
            for name in self._FLOAT_AGENT_FIELDS
        }
        self._agent_arrays['species_id'] = np.empty(0, dtype=np.int32)
        for name in self._UINT8_AGENT_FIELDS:
            self._agent_arrays[name] = np.empty(0, dtype=np.uint8)
        self._species_arrays = {
            name: np.empty(0, dtype=np.float64)
            for name in self._FLOAT_SPECIES_FIELDS
        }
        self._agent_objects = []

    @property
    def agent_capacity(self):
        return self._agent_capacity

    @property
    def species_capacity(self):
        return self._species_capacity

    @staticmethod
    def _grown_capacity(current, required):
        if required <= current:
            return current
        capacity = max(current, _MIN_BUFFER_CAPACITY)
        # Doubling bounds reallocations logarithmically as populations grow; see DOCUMENTATION.md §10.2.
        while capacity < required:
            capacity *= 2
        return capacity

    def _ensure_agent_capacity(self, required):
        capacity = self._grown_capacity(self._agent_capacity, required)
        if capacity == self._agent_capacity:
            return
        self._agent_arrays = {
            name: np.empty(capacity, dtype=np.float64)
            for name in self._FLOAT_AGENT_FIELDS
        }
        self._agent_arrays['species_id'] = np.empty(capacity, dtype=np.int32)
        for name in self._UINT8_AGENT_FIELDS:
            self._agent_arrays[name] = np.empty(capacity, dtype=np.uint8)
        self._agent_capacity = capacity

    def _ensure_species_capacity(self, required):
        capacity = self._grown_capacity(self._species_capacity, required)
        if capacity == self._species_capacity:
            return
        self._species_arrays = {
            name: np.empty(capacity, dtype=np.float64)
            for name in self._FLOAT_SPECIES_FIELDS
        }
        self._species_capacity = capacity

    def flatten(self, species):
        agent_count = sum(len(sp.members) for sp in species)
        self._ensure_agent_capacity(agent_count)
        self._ensure_species_capacity(len(species))
        self._agent_objects.clear()

        agent_index = 0
        for species_id, sp in enumerate(species):
            self._species_arrays['t6ss_activity'][species_id] = sp.t6ss_activity
            self._species_arrays['immunity_factor'][species_id] = sp.immunity_factor
            self._species_arrays['hgt_competence'][species_id] = sp.hgt_competence
            self._species_arrays['share_threshold'][species_id] = sp.share_thresh
            self._species_arrays['predation_threshold'][species_id] = sp.predation_thresh
            # Species order followed by member order is the deterministic contact identity contract;
            # see DOCUMENTATION.md §10.2.
            for bacterium in sp.members:
                self._agent_arrays['x'][agent_index] = bacterium.x
                self._agent_arrays['y'][agent_index] = bacterium.y
                self._agent_arrays['interaction_radius'][agent_index] = bacterium.size
                self._agent_arrays['base_radius'][agent_index] = bacterium.base_size
                self._agent_arrays['food_level'][agent_index] = bacterium.food_level
                self._agent_arrays['food_capacity'][agent_index] = bacterium.food_capacity
                self._agent_arrays['metabolic_burden_factor'][agent_index] = bacterium.metabolic_burden_factor
                self._agent_arrays['species_id'][agent_index] = species_id
                self._agent_arrays['is_persister'][agent_index] = bacterium.is_persister
                self._agent_arrays['phage_resistant'][agent_index] = bacterium.phage_resistant
                self._agent_arrays['event_flags'][agent_index] = 0
                self._agent_objects.append(bacterium)
                agent_index += 1

        agent_views = {
            name: array[:agent_count]
            for name, array in self._agent_arrays.items()
        }
        species_views = {
            name: array[:len(species)]
            for name, array in self._species_arrays.items()
        }
        return InteractionFrame(
            agents=AgentArrays(**agent_views),
            species=SpeciesArrays(**species_views),
            agent_objects=self._agent_objects,
        )


def _build_spatial_hash(x, y, cell_size):
    spatial_hash = {}
    for agent_index in range(x.shape[0]):
        # Signed floor division keeps negative coordinates in the geometrically correct hash cell;
        # see DOCUMENTATION.md §4.12.
        cell_x = math.floor(x[agent_index] / cell_size)
        cell_y = math.floor(y[agent_index] / cell_size)
        spatial_hash.setdefault((cell_x, cell_y), []).append(agent_index)
    return spatial_hash


def process_interactions_python(agents, species, dt, parameters=None, random_draw=None):
    """Mutate contact-owned arrays in deterministic historical traversal order."""
    parameters = parameters or ContactParameters()
    random_draw = random_draw or random.random
    agent_count = agents.x.shape[0]
    spatial_hash = _build_spatial_hash(
        agents.x,
        agents.y,
        parameters.hash_cell_size,
    )

    contact_pairs = 0
    t6ss_hits = 0
    t6ss_deaths = 0
    hgt_transfers = 0
    share_events = 0
    predation_events = 0

    for source_index in range(agent_count):
        if agents.food_level[source_index] < 0.0:
            continue
        source_cell_x = math.floor(agents.x[source_index] / parameters.hash_cell_size)
        source_cell_y = math.floor(agents.y[source_index] / parameters.hash_cell_size)
        source_species_id = int(agents.species_id[source_index])

        # The nx-then-ny bucket walk and insertion-ordered members preserve immediate mutation cascades;
        # see DOCUMENTATION.md §4.12 and §4.18.
        for offset_x in _NEIGHBOR_OFFSETS:
            neighbor_x = source_cell_x + offset_x
            for offset_y in _NEIGHBOR_OFFSETS:
                neighbor_y = source_cell_y + offset_y
                for target_index in spatial_hash.get((neighbor_x, neighbor_y), ()):
                    if target_index == source_index or agents.food_level[target_index] < 0.0:
                        continue

                    delta_x = agents.x[source_index] - agents.x[target_index]
                    delta_y = agents.y[source_index] - agents.y[target_index]
                    # The axis gate cheaply rejects cells beyond the documented near-contact neighborhood;
                    # see DOCUMENTATION.md §4.12.
                    if abs(delta_x) > parameters.axis_cutoff or abs(delta_y) > parameters.axis_cutoff:
                        continue

                    # Disk contact includes the soft margin and remains strict at exact tangency;
                    # see DOCUMENTATION.md §4.12.
                    distance_squared = delta_x ** 2 + delta_y ** 2
                    interaction_distance = (
                        agents.interaction_radius[source_index]
                        + agents.interaction_radius[target_index]
                        + parameters.contact_margin
                    )
                    if distance_squared >= interaction_distance ** 2:
                        continue

                    contact_pairs += 1
                    agents.event_flags[source_index] |= EVENT_CONTACT
                    target_species_id = int(agents.species_id[target_index])

                    # HGT precedes kin/competitor behavior so a successful recipient can donate later in
                    # this same ordered pass; P = 1 - exp(-k*c_d*c_r*dt), see DOCUMENTATION.md §4.18.
                    if agents.phage_resistant[source_index] and not agents.phage_resistant[target_index]:
                        donor_competence = species.hgt_competence[source_species_id]
                        recipient_competence = species.hgt_competence[target_species_id]
                        if donor_competence > 0.0 and recipient_competence > 0.0:
                            transfer_probability = 1.0 - math.exp(
                                -parameters.hgt_transfer_rate
                                * donor_competence
                                * recipient_competence
                                * dt
                            )
                            if random_draw() < transfer_probability:
                                agents.phage_resistant[target_index] = 1
                                agents.metabolic_burden_factor[target_index] = (
                                    1.0 + parameters.hgt_metabolic_burden
                                )
                                agents.event_flags[target_index] |= EVENT_HGT_TRANSFER
                                hgt_transfers += 1

                    if source_species_id == target_species_id:
                        # Kin sharing transfers diff*dt without exceeding recipient storage capacity;
                        # see DOCUMENTATION.md §4.12.
                        if agents.food_level[source_index] > (
                            agents.food_capacity[source_index]
                            * species.share_threshold[source_species_id]
                        ):
                            food_difference = (
                                agents.food_level[source_index] - agents.food_level[target_index]
                            )
                            if food_difference > 0.0:
                                transfer_amount = food_difference * dt
                                remaining_capacity = (
                                    agents.food_capacity[target_index]
                                    - agents.food_level[target_index]
                                )
                                actual_transfer = min(transfer_amount, remaining_capacity)
                                if actual_transfer > 0.0:
                                    agents.food_level[source_index] -= actual_transfer
                                    agents.food_level[target_index] += actual_transfer
                                    agents.event_flags[source_index] |= EVENT_SHARE
                                    agents.event_flags[target_index] |= EVENT_SHARE
                                    share_events += 1
                        continue

                    # T6SS dose scales with attacker area and elapsed seconds, then target immunity and
                    # persistence attenuate it; see DOCUMENTATION.md §4.12.
                    target_food_before_t6ss = agents.food_level[target_index]
                    base_damage = (
                        species.t6ss_activity[source_species_id]
                        * agents.interaction_radius[source_index] ** 2
                        * parameters.t6ss_damage_scaling
                        * dt
                    )
                    immunity_defense = min(
                        parameters.t6ss_max_immunity_defense,
                        (
                            species.immunity_factor[target_species_id]
                            * agents.interaction_radius[target_index]
                        ) / (
                            agents.interaction_radius[source_index]
                            + parameters.t6ss_immunity_size_epsilon
                        ),
                    )
                    actual_damage = base_damage * (1.0 - immunity_defense)
                    if agents.is_persister[target_index]:
                        actual_damage *= parameters.persister_resistance_factor
                    agents.food_level[target_index] -= actual_damage
                    if actual_damage > 0.0:
                        agents.event_flags[target_index] |= EVENT_T6SS_HIT
                        t6ss_hits += 1
                        if target_food_before_t6ss > 0.0 and agents.food_level[target_index] <= 0.0:
                            agents.event_flags[target_index] |= EVENT_T6SS_LETHAL
                            t6ss_deaths += 1

                    # A sufficiently larger hungry predator recovers remaining cytosol plus labile dry
                    # biomass, capped by its own storage; see DOCUMENTATION.md §4.12.
                    if agents.food_level[source_index] < (
                        agents.food_capacity[source_index]
                        * species.predation_threshold[source_species_id]
                    ):
                        if agents.interaction_radius[source_index] > (
                            agents.interaction_radius[target_index]
                            * parameters.phagocytosis_size_ratio
                        ):
                            prey_volume = agents.base_radius[target_index] ** 3
                            prey_necromass = (
                                prey_volume
                                * parameters.cell_wet_mass
                                * parameters.cell_dry_mass_fraction
                                * parameters.labile_necromass_fraction
                            )
                            digestible_gain = (
                                max(0.0, agents.food_level[target_index])
                                * parameters.cytosol_recovery_fraction
                                + prey_necromass
                            )
                            agents.food_level[source_index] = min(
                                agents.food_capacity[source_index],
                                agents.food_level[source_index] + digestible_gain,
                            )
                            agents.food_level[target_index] = parameters.consumed_prey_sentinel
                            agents.event_flags[source_index] |= EVENT_PREDATION
                            agents.event_flags[target_index] |= EVENT_PREDATION
                            predation_events += 1

    return InteractionSummary(
        contact_pairs=contact_pairs,
        t6ss_hits=t6ss_hits,
        t6ss_deaths=t6ss_deaths,
        hgt_transfers=hgt_transfers,
        share_events=share_events,
        predation_events=predation_events,
    )


class PhysicsOrchestrator:
    """Flatten, dispatch, validate, and commit contact-owned bacterium state."""

    def __init__(self, backend='auto', buffers=None, parameters=None):
        if backend not in _BACKENDS:
            raise ValueError("backend must be one of 'auto', 'python', or 'cpp'")
        self._backend = backend
        self._buffers = buffers or PhysicsBuffers()
        self._parameters = parameters or ContactParameters()
        self._backend_used = None
        self._native_import_attempted = False
        self._native_module = None
        self._native_callable = None
        self._native_load_error = None
        self._last_summary = InteractionSummary()

    @property
    def buffers(self):
        return self._buffers

    @property
    def backend_used(self):
        return self._backend_used

    @property
    def native_available(self):
        if not self._native_import_attempted:
            return None
        return self._native_callable is not None

    @property
    def last_summary(self):
        return self._last_summary

    def process(self, species, dt):
        frame = self._buffers.flatten(species)
        if self._backend == 'python':
            summary = process_interactions_python(
                frame.agents,
                frame.species,
                dt,
                self._parameters,
                random.random,
            )
            backend_used = 'python'
        else:
            try:
                native_callable = self._load_native_backend()
            except _NativePhysicsUnavailable:
                if self._backend == 'cpp':
                    raise
                summary = process_interactions_python(
                    frame.agents,
                    frame.species,
                    dt,
                    self._parameters,
                    random.random,
                )
                backend_used = 'python'
            else:
                summary = self._process_cpp(native_callable, frame, dt)
                backend_used = 'cpp'

        # Objects remain untouched until dispatch and native-result validation both succeed;
        # see DOCUMENTATION.md §10.2.
        self._write_back(frame)
        self._backend_used = backend_used
        self._last_summary = summary
        return summary

    def _load_native_backend(self):
        if self._native_import_attempted:
            if self._native_load_error is not None:
                raise self._native_load_error
            return self._native_callable

        self._native_import_attempted = True
        try:
            self._native_module = importlib.import_module(_NATIVE_MODULE_NAME)
        except ModuleNotFoundError as exc:
            if exc.name != _NATIVE_MODULE_NAME:
                self._native_load_error = exc
                raise
            error = _NativePhysicsUnavailable(
                "The 'bacteria_c_kernel' extension is unavailable; build it or use physics_backend='python'."
            )
            self._native_load_error = error
            raise error from exc
        except Exception as exc:
            self._native_load_error = exc
            raise

        try:
            native_callable = getattr(self._native_module, 'process_interactions')
        except AttributeError as exc:
            error = _NativePhysicsUnavailable(
                "The installed 'bacteria_c_kernel' has no 'process_interactions' symbol; "
                "build the Phase B extension or use physics_backend='python'."
            )
            self._native_load_error = error
            raise error from exc
        if not callable(native_callable):
            error = TypeError("bacteria_c_kernel.process_interactions must be callable")
            self._native_load_error = error
            raise error

        self._native_callable = native_callable
        return native_callable

    def _process_cpp(self, native_callable, frame, dt):
        agents = frame.agents
        species = frame.species
        parameters = self._parameters
        # The callback lets native code draw only for HGT-eligible ordered contacts, preserving
        # Python's conditional global RNG sequence and immediate cascades; see DOCUMENTATION.md §4.18.
        result = native_callable(
            agents.x,
            agents.y,
            agents.interaction_radius,
            agents.base_radius,
            agents.food_level,
            agents.food_capacity,
            agents.metabolic_burden_factor,
            agents.species_id,
            agents.is_persister,
            agents.phage_resistant,
            agents.event_flags,
            species.t6ss_activity,
            species.immunity_factor,
            species.hgt_competence,
            species.share_threshold,
            species.predation_threshold,
            dt,
            parameters.hash_cell_size,
            parameters.axis_cutoff,
            parameters.contact_margin,
            parameters.hgt_transfer_rate,
            parameters.hgt_metabolic_burden,
            parameters.t6ss_damage_scaling,
            parameters.t6ss_max_immunity_defense,
            parameters.t6ss_immunity_size_epsilon,
            parameters.persister_resistance_factor,
            parameters.phagocytosis_size_ratio,
            parameters.cell_wet_mass,
            parameters.cell_dry_mass_fraction,
            parameters.labile_necromass_fraction,
            parameters.cytosol_recovery_fraction,
            parameters.consumed_prey_sentinel,
            random.random,
        )
        return InteractionSummary.from_native(result)

    @staticmethod
    def _write_back(frame):
        for agent_index, bacterium in enumerate(frame.agent_objects):
            bacterium.food_level = float(frame.agents.food_level[agent_index])
            bacterium.phage_resistant = bool(frame.agents.phage_resistant[agent_index])
            bacterium.metabolic_burden_factor = float(
                frame.agents.metabolic_burden_factor[agent_index]
            )
