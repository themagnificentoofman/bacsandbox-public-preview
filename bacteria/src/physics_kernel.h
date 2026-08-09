#ifndef BACSANDBOX_PHYSICS_KERNEL_H
#define BACSANDBOX_PHYSICS_KERNEL_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace bacsandbox::physics {

template <typename Element>
struct ConstBufferView {
    const Element* data;
    std::size_t size;

    const Element& operator[](std::size_t index) const noexcept {
        return data[index];
    }
};

template <typename Element>
struct MutableBufferView {
    Element* data;
    std::size_t size;

    Element& operator[](std::size_t index) const noexcept {
        return data[index];
    }
};

enum EventBit : std::uint8_t {
    // Marks an ordered geometric contact accepted by the strict disk test; see DOCUMENTATION.md §4.12.
    event_contact = 1U << 0,
    // Marks receipt of a positive T6SS dose; see DOCUMENTATION.md §4.12.
    event_t6ss_hit = 1U << 1,
    // Marks a positive-to-nonpositive food transition caused by T6SS; see DOCUMENTATION.md §4.12.
    event_t6ss_lethal = 1U << 2,
    // Marks acquisition of phage resistance through HGT; see DOCUMENTATION.md §4.18.
    event_hgt_transfer = 1U << 3,
    // Marks participation in a successful same-species nutrient transfer; see DOCUMENTATION.md §4.12.
    event_share = 1U << 4,
    // Marks predator and consumed prey in a successful phagocytosis event; see DOCUMENTATION.md §4.12.
    event_predation = 1U << 5,
};

constexpr std::uint8_t all_event_bits =
    event_contact | event_t6ss_hit | event_t6ss_lethal |
    event_hgt_transfer | event_share | event_predation;

struct AgentBuffers {
    // Positions locate agents in the signed floor spatial grid; see DOCUMENTATION.md §10.2.
    ConstBufferView<double> x;
    ConstBufferView<double> y;
    // Interaction radius defines the effective contact disk and attacker area; see DOCUMENTATION.md §4.12.
    ConstBufferView<double> interaction_radius;
    // Base radius is cubed to estimate prey biomass independently of transient size; see DOCUMENTATION.md §4.12.
    ConstBufferView<double> base_radius;
    // Food is the mutable energy state used by sharing, damage, and predation; see DOCUMENTATION.md §4.12.
    MutableBufferView<double> food_level;
    // Capacity bounds nutrient sharing and predator recovery; see DOCUMENTATION.md §4.12.
    ConstBufferView<double> food_capacity;
    // Metabolic burden records the resistance cost after successful HGT; see DOCUMENTATION.md §4.18.
    MutableBufferView<double> metabolic_burden_factor;
    // Species IDs select per-species contact traits without object pointers; see DOCUMENTATION.md §10.2.
    ConstBufferView<std::int32_t> species_id;
    // Persister status attenuates T6SS damage for dormant targets; see DOCUMENTATION.md §4.12.
    ConstBufferView<std::uint8_t> is_persister;
    // Resistance status is both the HGT donor gate and mutable recipient state; see DOCUMENTATION.md §4.18.
    MutableBufferView<std::uint8_t> phage_resistant;
    // Event bits retain all successful ordered-contact outcomes for the tick; see DOCUMENTATION.md §10.2.
    MutableBufferView<std::uint8_t> event_flags;
};

struct SpeciesBuffers {
    // T6SS activity scales attacker-area damage per elapsed second; see DOCUMENTATION.md §4.12.
    ConstBufferView<double> t6ss_activity;
    // Immunity factor reduces T6SS dose according to target-to-attacker size; see DOCUMENTATION.md §4.12.
    ConstBufferView<double> immunity_factor;
    // HGT competence gates donor and recipient transfer probability; see DOCUMENTATION.md §4.18.
    ConstBufferView<double> hgt_competence;
    // Sharing threshold gates donation by source fullness; see DOCUMENTATION.md §4.12.
    ConstBufferView<double> share_threshold;
    // Predation threshold gates phagocytosis by source hunger; see DOCUMENTATION.md §4.12.
    ConstBufferView<double> predation_threshold;
};

struct InteractionParameters {
    // Elapsed seconds scale all contact rates independently of frame rate; see DOCUMENTATION.md §4.12 and §4.18.
    double dt;
    // Hash-cell width defines only broadphase membership, not contact distance; see DOCUMENTATION.md §10.2.
    double hash_cell_size;
    // The inclusive axis cutoff rejects distant candidate pairs before the disk test; see DOCUMENTATION.md §4.12.
    double axis_cutoff;
    // The soft contact margin expands the summed interaction radii; see DOCUMENTATION.md §4.12.
    double contact_margin;
    // HGT rate sets the exponential transfer hazard per elapsed second; see DOCUMENTATION.md §4.18.
    double hgt_transfer_rate;
    // HGT burden raises the recipient's metabolic cost after resistance acquisition; see DOCUMENTATION.md §4.18.
    double hgt_metabolic_burden;
    // T6SS scaling converts attacker area and activity into food damage per second; see DOCUMENTATION.md §4.12.
    double t6ss_damage_scaling;
    // Maximum immunity defense prevents attenuation from exceeding the model ceiling; see DOCUMENTATION.md §4.12.
    double t6ss_max_immunity_defense;
    // Size epsilon keeps the immunity ratio defined for a zero-radius attacker; see DOCUMENTATION.md §4.12.
    double t6ss_immunity_size_epsilon;
    // Persister factor attenuates damage to dormant targets; see DOCUMENTATION.md §4.12.
    double persister_resistance_factor;
    // Phagocytosis ratio requires a predator to exceed prey interaction size; see DOCUMENTATION.md §4.12.
    double phagocytosis_size_ratio;
    // Wet mass density converts prey reference volume into cellular mass; see DOCUMENTATION.md §4.12.
    double cell_wet_mass;
    // Dry mass fraction selects structural biomass from wet cellular mass; see DOCUMENTATION.md §4.12.
    double cell_dry_mass_fraction;
    // Labile fraction selects the recoverable share of prey dry biomass; see DOCUMENTATION.md §4.12.
    double labile_necromass_fraction;
    // Cytosol recovery fraction selects the recoverable share of prey food; see DOCUMENTATION.md §4.12.
    double cytosol_recovery_fraction;
    // The negative sentinel makes consumed prey immediately inactive in later ordered pairs; see DOCUMENTATION.md §4.12.
    double consumed_prey_sentinel;
};

struct CellKey {
    // Signed floor-grid coordinates preserve negative-position neighborhoods; see DOCUMENTATION.md §4.12.
    std::int64_t x;
    std::int64_t y;
};

struct SpatialGrid {
    // Source cells permit deterministic adjacent-key lookup without recomputing floor conversion; see DOCUMENTATION.md §10.2.
    std::vector<CellKey> source_cells;
    // Sorted unique keys and offsets form a contiguous sparse grid; see DOCUMENTATION.md §10.2.
    std::vector<CellKey> keys;
    std::vector<std::size_t> offsets;
    // Bucket members retain original flattened-agent order; see DOCUMENTATION.md §10.2.
    std::vector<std::size_t> members;
};

struct OrderedPair {
    std::size_t source;
    std::size_t target;
};

struct ContactPairs {
    // Source offsets partition deterministic source-owned pair slices; see DOCUMENTATION.md §10.2.
    std::vector<std::size_t> source_offsets;
    std::vector<OrderedPair> pairs;
};

struct InteractionSummary {
    // Counts ordered pairs that remain active when serial resolution reaches them; see DOCUMENTATION.md §10.2.
    std::size_t contact_pairs = 0;
    // Counts strictly positive T6SS doses applied to active targets; see DOCUMENTATION.md §4.12.
    std::size_t t6ss_hits = 0;
    // Counts T6SS positive-to-nonpositive target transitions; see DOCUMENTATION.md §4.12.
    std::size_t t6ss_deaths = 0;
    // Counts false-to-true resistance acquisitions through HGT; see DOCUMENTATION.md §4.18.
    std::size_t hgt_transfers = 0;
    // Counts successful capacity-bounded kin nutrient transfers; see DOCUMENTATION.md §4.12.
    std::size_t share_events = 0;
    // Counts successful size- and hunger-gated prey consumptions; see DOCUMENTATION.md §4.12.
    std::size_t predation_events = 0;
};

using RandomDraw = std::function<double()>;

SpatialGrid build_spatial_grid(
    const AgentBuffers& agents,
    const InteractionParameters& parameters);
ContactPairs discover_contact_pairs(
    const AgentBuffers& agents,
    const InteractionParameters& parameters);
InteractionSummary resolve_interactions(
    const AgentBuffers& agents,
    const SpeciesBuffers& species,
    const InteractionParameters& parameters,
    const ContactPairs& contacts,
    const RandomDraw& random_draw);

}  // namespace bacsandbox::physics

#endif