#include "physics_kernel.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#ifndef BACSANDBOX_OPENMP_ENABLED
#define BACSANDBOX_OPENMP_ENABLED 0
#endif

namespace bacsandbox::physics {

namespace {

struct SortableMember {
    CellKey key;
    std::size_t original_index;
};

bool is_finite_float64(double value) noexcept {
    static_assert(sizeof(double) == sizeof(std::uint64_t), "Physics buffers require IEEE-754 float64");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    constexpr std::uint64_t exponent_mask = 0x7ff0000000000000ULL;
    return (bits & exponent_mask) != exponent_mask;
}

bool cell_key_less(const CellKey& left, const CellKey& right) noexcept {
    return left.x < right.x || (left.x == right.x && left.y < right.y);
}

bool cell_key_equal(const CellKey& left, const CellKey& right) noexcept {
    return left.x == right.x && left.y == right.y;
}

CellKey floor_cell(double x, double y, double cell_size) {
    // Mathematical floor preserves signed geometry, and the strict conversion bounds reserve room
    // for each adjacent-cell offset; see DOCUMENTATION.md §4.12 and §10.2.
    const double scaled_x = x / cell_size;
    const double scaled_y = y / cell_size;
    if (!is_finite_float64(scaled_x) || !is_finite_float64(scaled_y)) {
        throw std::overflow_error("physics coordinate/cell conversion is outside the finite float64 range");
    }
    const double floored_x = std::floor(scaled_x);
    const double floored_y = std::floor(scaled_y);
    const double minimum_cell = static_cast<double>(std::numeric_limits<std::int64_t>::min());
    const double maximum_cell_exclusive = -minimum_cell;
    if (floored_x <= minimum_cell || floored_x >= maximum_cell_exclusive ||
        floored_y <= minimum_cell || floored_y >= maximum_cell_exclusive) {
        throw std::overflow_error("physics cell coordinate cannot represent all adjacent int64 keys");
    }
    return {
        static_cast<std::int64_t>(floored_x),
        static_cast<std::int64_t>(floored_y),
    };
}

template <typename Visitor>
void for_each_geometric_candidate(
    std::size_t source_index,
    const AgentBuffers& agents,
    const InteractionParameters& parameters,
    const SpatialGrid& grid,
    Visitor&& visitor) noexcept {
    const CellKey source_cell = grid.source_cells[source_index];

    // The nx-then-ny adjacent-key scan and original-order buckets preserve Python pair order;
    // see DOCUMENTATION.md §4.12 and §10.2.
    for (std::int64_t offset_x = -1; offset_x <= 1; ++offset_x) {
        for (std::int64_t offset_y = -1; offset_y <= 1; ++offset_y) {
            const CellKey neighbor_key{
                source_cell.x + offset_x,
                source_cell.y + offset_y,
            };
            const auto key_iterator = std::lower_bound(
                grid.keys.begin(),
                grid.keys.end(),
                neighbor_key,
                cell_key_less);
            if (key_iterator == grid.keys.end() || !cell_key_equal(*key_iterator, neighbor_key)) {
                continue;
            }

            const std::size_t key_index = static_cast<std::size_t>(
                key_iterator - grid.keys.begin());
            for (std::size_t member_offset = grid.offsets[key_index];
                 member_offset < grid.offsets[key_index + 1];
                 ++member_offset) {
                const std::size_t target_index = grid.members[member_offset];
                if (target_index == source_index) {
                    continue;
                }

                const double delta_x = agents.x[source_index] - agents.x[target_index];
                const double delta_y = agents.y[source_index] - agents.y[target_index];
                // The axis gate is inclusive at its boundary and only rejects values above it;
                // see DOCUMENTATION.md §4.12.
                if (std::fabs(delta_x) > parameters.axis_cutoff ||
                    std::fabs(delta_y) > parameters.axis_cutoff) {
                    continue;
                }

                // Contact uses a strict disk inequality with the Python-supplied soft margin;
                // see DOCUMENTATION.md §4.12.
                const double distance_squared = delta_x * delta_x + delta_y * delta_y;
                const double interaction_distance =
                    agents.interaction_radius[source_index] +
                    agents.interaction_radius[target_index] +
                    parameters.contact_margin;
                if (distance_squared >= interaction_distance * interaction_distance) {
                    continue;
                }
                visitor(target_index);
            }
        }
    }
}

void set_event(MutableBufferView<std::uint8_t> event_flags, std::size_t index, EventBit event) noexcept {
    event_flags[index] = static_cast<std::uint8_t>(event_flags[index] | event);
}

}  // namespace

SpatialGrid build_spatial_grid(
    const AgentBuffers& agents,
    const InteractionParameters& parameters) {
    if (agents.x.size != agents.y.size) {
        throw std::invalid_argument("physics x and y buffer sizes differ");
    }

    const std::size_t agent_count = agents.x.size;
    SpatialGrid grid;
    grid.source_cells.resize(agent_count);
    std::vector<SortableMember> sorted_members;
    sorted_members.reserve(agent_count);
    for (std::size_t agent_index = 0; agent_index < agent_count; ++agent_index) {
        const CellKey cell = floor_cell(
            agents.x[agent_index],
            agents.y[agent_index],
            parameters.hash_cell_size);
        grid.source_cells[agent_index] = cell;
        sorted_members.push_back({cell, agent_index});
    }

    // Stable tuple order yields sorted keys while retaining flattened source order within buckets;
    // see DOCUMENTATION.md §10.2.
    std::stable_sort(
        sorted_members.begin(),
        sorted_members.end(),
        [](const SortableMember& left, const SortableMember& right) noexcept {
            if (cell_key_less(left.key, right.key)) {
                return true;
            }
            if (cell_key_less(right.key, left.key)) {
                return false;
            }
            return left.original_index < right.original_index;
        });

    grid.keys.reserve(agent_count);
    grid.offsets.reserve(agent_count + 1);
    grid.members.reserve(agent_count);
    for (const SortableMember& member : sorted_members) {
        if (grid.keys.empty() || !cell_key_equal(grid.keys.back(), member.key)) {
            grid.keys.push_back(member.key);
            grid.offsets.push_back(grid.members.size());
        }
        grid.members.push_back(member.original_index);
    }
    grid.offsets.push_back(grid.members.size());
    return grid;
}

ContactPairs discover_contact_pairs(
    const AgentBuffers& agents,
    const InteractionParameters& parameters) {
    const SpatialGrid grid = build_spatial_grid(agents, parameters);
    const std::size_t agent_count = agents.x.size;
    if (agent_count > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error("physics agent count exceeds the OpenMP loop range");
    }

    std::vector<std::size_t> source_counts(agent_count, 0);
#if BACSANDBOX_OPENMP_ENABLED
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t source = 0; source < static_cast<std::int64_t>(agent_count); ++source) {
        const std::size_t source_index = static_cast<std::size_t>(source);
        std::size_t count = 0;
        for_each_geometric_candidate(
            source_index,
            agents,
            parameters,
            grid,
            [&count](std::size_t) noexcept {
                ++count;
            });
        source_counts[source_index] = count;
    }

    ContactPairs contacts;
    contacts.source_offsets.resize(agent_count + 1, 0);
    for (std::size_t source_index = 0; source_index < agent_count; ++source_index) {
        const std::size_t preceding = contacts.source_offsets[source_index];
        if (source_counts[source_index] > std::numeric_limits<std::size_t>::max() - preceding) {
            throw std::overflow_error("physics ordered-pair count exceeds the native size range");
        }
        contacts.source_offsets[source_index + 1] = preceding + source_counts[source_index];
    }
    contacts.pairs.resize(contacts.source_offsets.back());

#if BACSANDBOX_OPENMP_ENABLED
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t source = 0; source < static_cast<std::int64_t>(agent_count); ++source) {
        const std::size_t source_index = static_cast<std::size_t>(source);
        std::size_t write_offset = contacts.source_offsets[source_index];
        for_each_geometric_candidate(
            source_index,
            agents,
            parameters,
            grid,
            [&contacts, &write_offset, source_index](std::size_t target_index) noexcept {
                contacts.pairs[write_offset] = {source_index, target_index};
                ++write_offset;
            });
    }
    return contacts;
}

InteractionSummary resolve_interactions(
    const AgentBuffers& agents,
    const SpeciesBuffers& species,
    const InteractionParameters& parameters,
    const ContactPairs& contacts,
    const RandomDraw& random_draw) {
    InteractionSummary summary;
    for (const OrderedPair& pair : contacts.pairs) {
        const std::size_t source_index = pair.source;
        const std::size_t target_index = pair.target;
        // Mutable liveness is rechecked in pair order so consumed prey stop participating immediately;
        // see DOCUMENTATION.md §4.12 and §10.2.
        if (agents.food_level[source_index] < 0.0 || agents.food_level[target_index] < 0.0) {
            continue;
        }

        ++summary.contact_pairs;
        set_event(agents.event_flags, source_index, event_contact);
        const std::size_t source_species = static_cast<std::size_t>(agents.species_id[source_index]);
        const std::size_t target_species = static_cast<std::size_t>(agents.species_id[target_index]);

        // HGT resolves first so a newly resistant recipient can donate later in the same ordered pass;
        // P = 1 - exp(-rate * donor competence * recipient competence * dt), see DOCUMENTATION.md §4.18.
        if (agents.phage_resistant[source_index] != 0 &&
            agents.phage_resistant[target_index] == 0) {
            const double donor_competence = species.hgt_competence[source_species];
            const double recipient_competence = species.hgt_competence[target_species];
            if (donor_competence > 0.0 && recipient_competence > 0.0) {
                const double transfer_probability = 1.0 - std::exp(
                    -parameters.hgt_transfer_rate * donor_competence *
                    recipient_competence * parameters.dt);
                const double draw = random_draw();
                if (!is_finite_float64(draw) || draw < 0.0 || draw >= 1.0) {
                    throw std::invalid_argument("random_draw must return a finite value in [0, 1)");
                }
                if (draw < transfer_probability) {
                    agents.phage_resistant[target_index] = 1;
                    agents.metabolic_burden_factor[target_index] =
                        1.0 + parameters.hgt_metabolic_burden;
                    set_event(agents.event_flags, target_index, event_hgt_transfer);
                    ++summary.hgt_transfers;
                }
            }
        }

        if (source_species == target_species) {
            // Kin sharing transfers the food difference times elapsed seconds without exceeding
            // recipient capacity; see DOCUMENTATION.md §4.12.
            if (agents.food_level[source_index] >
                agents.food_capacity[source_index] * species.share_threshold[source_species]) {
                const double food_difference =
                    agents.food_level[source_index] - agents.food_level[target_index];
                if (food_difference > 0.0) {
                    const double transfer_amount = food_difference * parameters.dt;
                    const double remaining_capacity =
                        agents.food_capacity[target_index] - agents.food_level[target_index];
                    const double actual_transfer = std::min(transfer_amount, remaining_capacity);
                    if (actual_transfer > 0.0) {
                        agents.food_level[source_index] -= actual_transfer;
                        agents.food_level[target_index] += actual_transfer;
                        set_event(agents.event_flags, source_index, event_share);
                        set_event(agents.event_flags, target_index, event_share);
                        ++summary.share_events;
                    }
                }
            }
            continue;
        }

        // T6SS dose scales with attacker area and elapsed seconds, then target size-immunity and
        // persister state attenuate it; see DOCUMENTATION.md §4.12.
        const double target_food_before_t6ss = agents.food_level[target_index];
        const double source_radius = agents.interaction_radius[source_index];
        const double target_radius = agents.interaction_radius[target_index];
        const double base_damage =
            species.t6ss_activity[source_species] * source_radius * source_radius *
            parameters.t6ss_damage_scaling * parameters.dt;
        const double immunity_defense = std::min(
            parameters.t6ss_max_immunity_defense,
            species.immunity_factor[target_species] * target_radius /
                (source_radius + parameters.t6ss_immunity_size_epsilon));
        double actual_damage = base_damage * (1.0 - immunity_defense);
        if (agents.is_persister[target_index] != 0) {
            actual_damage *= parameters.persister_resistance_factor;
        }
        agents.food_level[target_index] -= actual_damage;
        if (actual_damage > 0.0) {
            set_event(agents.event_flags, target_index, event_t6ss_hit);
            ++summary.t6ss_hits;
            if (target_food_before_t6ss > 0.0 && agents.food_level[target_index] <= 0.0) {
                set_event(agents.event_flags, target_index, event_t6ss_lethal);
                ++summary.t6ss_deaths;
            }
        }

        // Hungry larger predators recover remaining cytosol and labile dry biomass, capped by
        // predator capacity; see DOCUMENTATION.md §4.12.
        if (agents.food_level[source_index] <
            agents.food_capacity[source_index] * species.predation_threshold[source_species]) {
            if (source_radius > target_radius * parameters.phagocytosis_size_ratio) {
                const double prey_base_radius = agents.base_radius[target_index];
                const double prey_volume =
                    prey_base_radius * prey_base_radius * prey_base_radius;
                const double prey_necromass =
                    prey_volume * parameters.cell_wet_mass *
                    parameters.cell_dry_mass_fraction *
                    parameters.labile_necromass_fraction;
                const double digestible_gain =
                    std::max(0.0, agents.food_level[target_index]) *
                        parameters.cytosol_recovery_fraction +
                    prey_necromass;
                agents.food_level[source_index] = std::min(
                    agents.food_capacity[source_index],
                    agents.food_level[source_index] + digestible_gain);
                agents.food_level[target_index] = parameters.consumed_prey_sentinel;
                set_event(agents.event_flags, source_index, event_predation);
                set_event(agents.event_flags, target_index, event_predation);
                ++summary.predation_events;
            }
        }
    }
    return summary;
}

}  // namespace bacsandbox::physics