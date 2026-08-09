#include "ecm_kernel.h"
#include "physics_kernel.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace py = pybind11;

#ifndef BACSANDBOX_OPENMP_ENABLED
#define BACSANDBOX_OPENMP_ENABLED 0
#endif

namespace bacsandbox::ecm {

namespace {

constexpr std::size_t layer_count = 4;

struct CoordinateHash {
    std::size_t operator()(const Coordinate& coordinate) const noexcept {
        constexpr unsigned coordinate_bits = std::numeric_limits<std::uint64_t>::digits;
        constexpr unsigned half_coordinate_bits = coordinate_bits / 2;
        const std::uint64_t x_bits = static_cast<std::uint64_t>(coordinate.x);
        const std::uint64_t y_bits = static_cast<std::uint64_t>(coordinate.y);
        const std::uint64_t rotated_y =
            (y_bits << half_coordinate_bits) |
            (y_bits >> (coordinate_bits - half_coordinate_bits));
        return static_cast<std::size_t>(x_bits ^ rotated_y);
    }
};

struct OrderedLayerResult {
    std::vector<Coordinate> coordinates;
    std::vector<double> values;
};

std::int64_t checked_decrement(std::int64_t coordinate) {
    if (coordinate == std::numeric_limits<std::int64_t>::min()) {
        throw std::overflow_error("ECM neighbor coordinate is below the int64 range");
    }
    return coordinate - 1;
}

std::int64_t checked_increment(std::int64_t coordinate) {
    if (coordinate == std::numeric_limits<std::int64_t>::max()) {
        throw std::overflow_error("ECM neighbor coordinate is above the int64 range");
    }
    return coordinate + 1;
}

void add_contribution(
    OrderedLayerResult& result,
    std::unordered_map<Coordinate, std::size_t, CoordinateHash>& coordinate_indices,
    const Coordinate& coordinate,
    double contribution) {
    const auto insertion = coordinate_indices.emplace(coordinate, result.coordinates.size());
    if (insertion.second) {
        result.coordinates.push_back(coordinate);
        result.values.push_back(contribution);
        return;
    }
    result.values[insertion.first->second] += contribution;
}

OrderedLayerResult process_layer(
    const LayerBatch& layers,
    const KernelParameters& parameters,
    std::size_t layer_index) {
    // The capped diffusion fraction keeps each explicit 2D stencil step nonnegative and stable; see DOCUMENTATION.md §4.9.
    const double alpha = std::min(
        parameters.alpha_max,
        std::max(
            0.0,
            parameters.diffusion[layer_index] * parameters.dt /
                (parameters.cell_size * parameters.cell_size)));
    // The clamped first-order decay multiplier prevents a timestep from removing more concentration than exists; see DOCUMENTATION.md §4.9.
    const double decay_multiplier = std::max(
        0.0,
        std::min(1.0, 1.0 - parameters.decay[layer_index] * parameters.dt));

    OrderedLayerResult accumulated;
    std::unordered_map<Coordinate, std::size_t, CoordinateHash> coordinate_indices;
    const std::size_t start = static_cast<std::size_t>(layers.offsets[layer_index]);
    const std::size_t stop = static_cast<std::size_t>(layers.offsets[layer_index + 1]);
    coordinate_indices.reserve(stop - start);

    // Ordered vectors reproduce Python dict first-insertion order while the hash map only accelerates accumulation; see DOCUMENTATION.md §4.9.
    for (std::size_t source_index = start; source_index < stop; ++source_index) {
        const double value = layers.values[source_index];
        if (value < parameters.min_concentration) {
            continue;
        }

        const double decayed = value * decay_multiplier;
        const double retained = decayed * (1.0 - 4.0 * alpha);
        const double neighbor_share = decayed * alpha;
        const Coordinate source = layers.coordinates[source_index];

        if (retained >= parameters.min_concentration) {
            add_contribution(accumulated, coordinate_indices, source, retained);
        }
        if (neighbor_share >= parameters.min_concentration) {
            add_contribution(
                accumulated,
                coordinate_indices,
                {checked_decrement(source.x), source.y},
                neighbor_share);
            add_contribution(
                accumulated,
                coordinate_indices,
                {checked_increment(source.x), source.y},
                neighbor_share);
            add_contribution(
                accumulated,
                coordinate_indices,
                {source.x, checked_decrement(source.y)},
                neighbor_share);
            add_contribution(
                accumulated,
                coordinate_indices,
                {source.x, checked_increment(source.y)},
                neighbor_share);
        }
    }

    OrderedLayerResult filtered;
    filtered.coordinates.reserve(accumulated.coordinates.size());
    filtered.values.reserve(accumulated.values.size());
    for (std::size_t index = 0; index < accumulated.values.size(); ++index) {
        if (accumulated.values[index] >= parameters.min_concentration) {
            filtered.coordinates.push_back(accumulated.coordinates[index]);
            filtered.values.push_back(accumulated.values[index]);
        }
    }
    return filtered;
}

void validate_layer_batch(const LayerBatch& layers) {
    if (layers.coordinates.size() != layers.values.size()) {
        throw std::invalid_argument("ECM coordinate and value counts differ");
    }
    if (layers.coordinates.size() >
        static_cast<std::uintmax_t>(std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error("ECM input row count exceeds the int64 offset range");
    }
    if (layers.offsets.front() != 0 ||
        layers.offsets.back() != static_cast<std::int64_t>(layers.coordinates.size())) {
        throw std::invalid_argument("ECM layer offsets must span all input rows");
    }
    for (std::size_t index = 0; index < layer_count; ++index) {
        if (layers.offsets[index] < 0 || layers.offsets[index] > layers.offsets[index + 1]) {
            throw std::invalid_argument("ECM layer offsets must be nonnegative and monotonic");
        }
    }
}

}  // namespace

LayerResult dispatch_process_layers(const LayerBatch& layers, const KernelParameters& parameters) {
    validate_layer_batch(layers);

    std::array<OrderedLayerResult, layer_count> layer_results;
    std::array<std::exception_ptr, layer_count> layer_errors{};
#if BACSANDBOX_OPENMP_ENABLED
#pragma omp parallel for schedule(static)
#endif
    for (int layer_index = 0; layer_index < static_cast<int>(layer_count); ++layer_index) {
        try {
            layer_results[static_cast<std::size_t>(layer_index)] =
                process_layer(layers, parameters, static_cast<std::size_t>(layer_index));
        } catch (...) {
            layer_errors[static_cast<std::size_t>(layer_index)] = std::current_exception();
        }
    }

    for (const std::exception_ptr& error : layer_errors) {
        if (error) {
            std::rethrow_exception(error);
        }
    }

    LayerResult flattened;
    flattened.offsets[0] = 0;
    std::size_t total_size = 0;
    const std::size_t maximum_offset =
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max());
    for (std::size_t layer_index = 0; layer_index < layer_count; ++layer_index) {
        const std::size_t layer_size = layer_results[layer_index].coordinates.size();
        if (layer_size > maximum_offset - total_size) {
            throw std::overflow_error("ECM output row count exceeds the int64 offset range");
        }
        total_size += layer_size;
        flattened.offsets[layer_index + 1] = static_cast<std::int64_t>(total_size);
    }

    flattened.coordinates.reserve(total_size);
    flattened.values.reserve(total_size);
    for (std::size_t layer_index = 0; layer_index < layer_count; ++layer_index) {
        OrderedLayerResult& layer_result = layer_results[layer_index];
        flattened.coordinates.insert(
            flattened.coordinates.end(),
            layer_result.coordinates.begin(),
            layer_result.coordinates.end());
        flattened.values.insert(
            flattened.values.end(),
            layer_result.values.begin(),
            layer_result.values.end());
    }
    return flattened;
}

bool openmp_enabled() noexcept {
    return BACSANDBOX_OPENMP_ENABLED != 0;
}

}  // namespace bacsandbox::ecm

namespace {

bool is_finite_float64(double value) noexcept {
    static_assert(sizeof(double) == sizeof(std::uint64_t), "ECM buffers require IEEE-754 float64");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    // IEEE-754 reserves an all-ones exponent for infinity/NaN, so this check remains valid under fast-math; see DOCUMENTATION.md §4.9.
    constexpr std::uint64_t exponent_mask = 0x7ff0000000000000ULL;
    return (bits & exponent_mask) != exponent_mask;
}

template <typename Element>
void require_c_array(const py::array& array, const char* name) {
    if (!array.dtype().is(py::dtype::of<Element>())) {
        throw py::type_error(std::string(name) + " has an invalid dtype");
    }
    if ((array.flags() & py::array::c_style) == 0) {
        throw py::value_error(std::string(name) + " must be C-contiguous");
    }
}

void require_finite_nonnegative(double value, const char* name) {
    if (!is_finite_float64(value) || value < 0.0) {
        throw py::value_error(std::string(name) + " must be finite and nonnegative");
    }
}

void require_finite(double value, const char* name) {
    if (!is_finite_float64(value)) {
        throw py::value_error(std::string(name) + " must be finite");
    }
}

void require_unit_interval(double value, const char* name) {
    if (!is_finite_float64(value) || value < 0.0 || value > 1.0) {
        throw py::value_error(std::string(name) + " must be finite and in [0, 1]");
    }
}

template <typename Element>
std::size_t require_vector(const py::array& array, const char* name) {
    require_c_array<Element>(array, name);
    if (array.ndim() != 1) {
        throw py::value_error(std::string(name) + " must have shape (N,)");
    }
    return static_cast<std::size_t>(array.shape(0));
}

void require_length(std::size_t actual, std::size_t expected, const char* name, const char* extent) {
    if (actual != expected) {
        throw py::value_error(
            std::string(name) + " must have shape (" + extent + ",)");
    }
}

void require_writable(const py::array& array, const char* name) {
    if (!array.writeable()) {
        throw py::value_error(std::string(name) + " must be writable");
    }
}

bacsandbox::ecm::LayerBatch read_layers(
    const py::array& coordinates,
    const py::array& values,
    const py::array& offsets) {
    require_c_array<std::int64_t>(coordinates, "coordinates");
    require_c_array<double>(values, "values");
    require_c_array<std::int64_t>(offsets, "offsets");

    if (coordinates.ndim() != 2 || coordinates.shape(1) != 2) {
        throw py::value_error("coordinates must have shape (N, 2)");
    }
    const py::ssize_t row_count = coordinates.shape(0);
    if (values.ndim() != 1 || values.shape(0) != row_count) {
        throw py::value_error("values must have shape (N,)");
    }
    if (offsets.ndim() != 1 || offsets.shape(0) != 5) {
        throw py::value_error("offsets must have shape (5,)");
    }

    const auto* coordinate_data = static_cast<const std::int64_t*>(coordinates.data());
    const auto* value_data = static_cast<const double*>(values.data());
    const auto* offset_data = static_cast<const std::int64_t*>(offsets.data());
    if (offset_data[0] != 0 || offset_data[4] != row_count) {
        throw py::value_error("offsets must start at zero and end at N");
    }
    for (std::size_t index = 0; index < 4; ++index) {
        if (offset_data[index] < 0 || offset_data[index] > offset_data[index + 1]) {
            throw py::value_error("offsets must be nonnegative and monotonic");
        }
    }

    bacsandbox::ecm::LayerBatch layers;
    layers.coordinates.reserve(static_cast<std::size_t>(row_count));
    layers.values.reserve(static_cast<std::size_t>(row_count));
    for (py::ssize_t row = 0; row < row_count; ++row) {
        if (!is_finite_float64(value_data[row])) {
            throw py::value_error("values must contain only finite concentrations");
        }
        layers.coordinates.push_back({coordinate_data[2 * row], coordinate_data[2 * row + 1]});
        layers.values.push_back(value_data[row]);
    }
    std::copy_n(offset_data, layers.offsets.size(), layers.offsets.begin());
    return layers;
}

bacsandbox::ecm::KernelParameters read_parameters(
    const py::array& diffusion,
    const py::array& decay,
    double dt,
    double cell_size,
    double alpha_max,
    double min_concentration) {
    require_c_array<double>(diffusion, "diffusion");
    require_c_array<double>(decay, "decay");
    if (diffusion.ndim() != 1 || diffusion.shape(0) != 4) {
        throw py::value_error("diffusion must have shape (4,)");
    }
    if (decay.ndim() != 1 || decay.shape(0) != 4) {
        throw py::value_error("decay must have shape (4,)");
    }

    const auto* diffusion_data = static_cast<const double*>(diffusion.data());
    const auto* decay_data = static_cast<const double*>(decay.data());
    bacsandbox::ecm::KernelParameters parameters;
    for (std::size_t index = 0; index < parameters.diffusion.size(); ++index) {
        require_finite_nonnegative(diffusion_data[index], "diffusion entries");
        require_finite_nonnegative(decay_data[index], "decay entries");
        parameters.diffusion[index] = diffusion_data[index];
        parameters.decay[index] = decay_data[index];
    }

    require_finite_nonnegative(dt, "dt");
    require_finite_nonnegative(alpha_max, "alpha_max");
    require_finite_nonnegative(min_concentration, "min_concentration");
    if (!is_finite_float64(cell_size) || cell_size <= 0.0) {
        throw py::value_error("cell_size must be finite and positive");
    }
    parameters.dt = dt;
    parameters.cell_size = cell_size;
    parameters.alpha_max = alpha_max;
    parameters.min_concentration = min_concentration;
    return parameters;
}

py::tuple make_python_result(bacsandbox::ecm::LayerResult result) {
    py::array_t<std::int64_t> coordinates({
        static_cast<py::ssize_t>(result.coordinates.size()),
        static_cast<py::ssize_t>(2),
    });
    py::array_t<double> values(static_cast<py::ssize_t>(result.values.size()));
    py::array_t<std::int64_t> offsets(static_cast<py::ssize_t>(result.offsets.size()));

    auto coordinate_output = coordinates.mutable_unchecked<2>();
    for (py::ssize_t row = 0; row < coordinate_output.shape(0); ++row) {
        coordinate_output(row, 0) = result.coordinates[static_cast<std::size_t>(row)].x;
        coordinate_output(row, 1) = result.coordinates[static_cast<std::size_t>(row)].y;
    }
    auto value_output = values.mutable_unchecked<1>();
    for (py::ssize_t row = 0; row < value_output.shape(0); ++row) {
        value_output(row) = result.values[static_cast<std::size_t>(row)];
    }
    auto offset_output = offsets.mutable_unchecked<1>();
    for (py::ssize_t index = 0; index < offset_output.shape(0); ++index) {
        offset_output(index) = result.offsets[static_cast<std::size_t>(index)];
    }

    return py::make_tuple(std::move(coordinates), std::move(values), std::move(offsets));
}

py::tuple process_layers_binding(
    const py::array& coordinates,
    const py::array& values,
    const py::array& offsets,
    const py::array& diffusion,
    const py::array& decay,
    double dt,
    double cell_size,
    double alpha_max,
    double min_concentration) {
    bacsandbox::ecm::LayerBatch layers = read_layers(coordinates, values, offsets);
    bacsandbox::ecm::KernelParameters parameters = read_parameters(
        diffusion,
        decay,
        dt,
        cell_size,
        alpha_max,
        min_concentration);

    bacsandbox::ecm::LayerResult result;
    {
        py::gil_scoped_release release;
        result = bacsandbox::ecm::dispatch_process_layers(layers, parameters);
    }
    return make_python_result(std::move(result));
}

py::tuple process_interactions_binding(
    const py::array& x,
    const py::array& y,
    const py::array& interaction_radius,
    const py::array& base_radius,
    py::array& food_level,
    const py::array& food_capacity,
    py::array& metabolic_burden_factor,
    const py::array& species_id,
    const py::array& is_persister,
    py::array& phage_resistant,
    py::array& event_flags,
    const py::array& t6ss_activity,
    const py::array& immunity_factor,
    const py::array& hgt_competence,
    const py::array& share_threshold,
    const py::array& predation_threshold,
    double dt,
    double hash_cell_size,
    double axis_cutoff,
    double contact_margin,
    double hgt_transfer_rate,
    double hgt_metabolic_burden,
    double t6ss_damage_scaling,
    double t6ss_max_immunity_defense,
    double t6ss_immunity_size_epsilon,
    double persister_resistance_factor,
    double phagocytosis_size_ratio,
    double cell_wet_mass,
    double cell_dry_mass_fraction,
    double labile_necromass_fraction,
    double cytosol_recovery_fraction,
    double consumed_prey_sentinel,
    const py::function& random_draw) {
    const std::size_t agent_count = require_vector<double>(x, "x");
    require_length(require_vector<double>(y, "y"), agent_count, "y", "N");
    require_length(
        require_vector<double>(interaction_radius, "interaction_radius"),
        agent_count,
        "interaction_radius",
        "N");
    require_length(
        require_vector<double>(base_radius, "base_radius"),
        agent_count,
        "base_radius",
        "N");
    require_length(
        require_vector<double>(food_level, "food_level"),
        agent_count,
        "food_level",
        "N");
    require_length(
        require_vector<double>(food_capacity, "food_capacity"),
        agent_count,
        "food_capacity",
        "N");
    require_length(
        require_vector<double>(metabolic_burden_factor, "metabolic_burden_factor"),
        agent_count,
        "metabolic_burden_factor",
        "N");
    require_length(
        require_vector<std::int32_t>(species_id, "species_id"),
        agent_count,
        "species_id",
        "N");
    require_length(
        require_vector<std::uint8_t>(is_persister, "is_persister"),
        agent_count,
        "is_persister",
        "N");
    require_length(
        require_vector<std::uint8_t>(phage_resistant, "phage_resistant"),
        agent_count,
        "phage_resistant",
        "N");
    require_length(
        require_vector<std::uint8_t>(event_flags, "event_flags"),
        agent_count,
        "event_flags",
        "N");

    const std::size_t species_count = require_vector<double>(t6ss_activity, "t6ss_activity");
    require_length(
        require_vector<double>(immunity_factor, "immunity_factor"),
        species_count,
        "immunity_factor",
        "S");
    require_length(
        require_vector<double>(hgt_competence, "hgt_competence"),
        species_count,
        "hgt_competence",
        "S");
    require_length(
        require_vector<double>(share_threshold, "share_threshold"),
        species_count,
        "share_threshold",
        "S");
    require_length(
        require_vector<double>(predation_threshold, "predation_threshold"),
        species_count,
        "predation_threshold",
        "S");

    require_writable(food_level, "food_level");
    require_writable(metabolic_burden_factor, "metabolic_burden_factor");
    require_writable(phage_resistant, "phage_resistant");
    require_writable(event_flags, "event_flags");

    const auto* x_data = static_cast<const double*>(x.data());
    const auto* y_data = static_cast<const double*>(y.data());
    const auto* interaction_radius_data = static_cast<const double*>(interaction_radius.data());
    const auto* base_radius_data = static_cast<const double*>(base_radius.data());
    auto* food_level_data = static_cast<double*>(food_level.mutable_data());
    const auto* food_capacity_data = static_cast<const double*>(food_capacity.data());
    auto* metabolic_burden_data = static_cast<double*>(metabolic_burden_factor.mutable_data());
    const auto* species_id_data = static_cast<const std::int32_t*>(species_id.data());
    const auto* is_persister_data = static_cast<const std::uint8_t*>(is_persister.data());
    auto* phage_resistant_data = static_cast<std::uint8_t*>(phage_resistant.mutable_data());
    auto* event_flags_data = static_cast<std::uint8_t*>(event_flags.mutable_data());

    for (std::size_t index = 0; index < agent_count; ++index) {
        require_finite(x_data[index], "x entries");
        require_finite(y_data[index], "y entries");
        require_finite_nonnegative(interaction_radius_data[index], "interaction_radius entries");
        require_finite_nonnegative(base_radius_data[index], "base_radius entries");
        require_finite(food_level_data[index], "food_level entries");
        require_finite_nonnegative(food_capacity_data[index], "food_capacity entries");
        require_finite_nonnegative(metabolic_burden_data[index], "metabolic_burden_factor entries");
        if (species_id_data[index] < 0 ||
            static_cast<std::size_t>(species_id_data[index]) >= species_count) {
            throw py::value_error("species_id entries must index the species arrays");
        }
        if (is_persister_data[index] > 1) {
            throw py::value_error("is_persister entries must be 0 or 1");
        }
        if (phage_resistant_data[index] > 1) {
            throw py::value_error("phage_resistant entries must be 0 or 1");
        }
        if ((event_flags_data[index] & ~bacsandbox::physics::all_event_bits) != 0) {
            throw py::value_error("event_flags entries contain unknown event bits");
        }
    }

    const auto* t6ss_activity_data = static_cast<const double*>(t6ss_activity.data());
    const auto* immunity_factor_data = static_cast<const double*>(immunity_factor.data());
    const auto* hgt_competence_data = static_cast<const double*>(hgt_competence.data());
    const auto* share_threshold_data = static_cast<const double*>(share_threshold.data());
    const auto* predation_threshold_data = static_cast<const double*>(predation_threshold.data());
    for (std::size_t index = 0; index < species_count; ++index) {
        require_unit_interval(t6ss_activity_data[index], "t6ss_activity entries");
        require_unit_interval(immunity_factor_data[index], "immunity_factor entries");
        require_unit_interval(hgt_competence_data[index], "hgt_competence entries");
        require_unit_interval(share_threshold_data[index], "share_threshold entries");
        require_unit_interval(predation_threshold_data[index], "predation_threshold entries");
    }

    require_finite_nonnegative(dt, "dt");
    if (!is_finite_float64(hash_cell_size) || hash_cell_size <= 0.0) {
        throw py::value_error("hash_cell_size must be finite and positive");
    }
    require_finite_nonnegative(axis_cutoff, "axis_cutoff");
    require_finite_nonnegative(contact_margin, "contact_margin");
    require_finite_nonnegative(hgt_transfer_rate, "hgt_transfer_rate");
    require_finite_nonnegative(hgt_metabolic_burden, "hgt_metabolic_burden");
    require_finite_nonnegative(t6ss_damage_scaling, "t6ss_damage_scaling");
    require_unit_interval(t6ss_max_immunity_defense, "t6ss_max_immunity_defense");
    if (!is_finite_float64(t6ss_immunity_size_epsilon) ||
        t6ss_immunity_size_epsilon <= 0.0) {
        throw py::value_error("t6ss_immunity_size_epsilon must be finite and positive");
    }
    require_unit_interval(persister_resistance_factor, "persister_resistance_factor");
    require_finite_nonnegative(phagocytosis_size_ratio, "phagocytosis_size_ratio");
    require_finite_nonnegative(cell_wet_mass, "cell_wet_mass");
    require_unit_interval(cell_dry_mass_fraction, "cell_dry_mass_fraction");
    require_unit_interval(labile_necromass_fraction, "labile_necromass_fraction");
    require_unit_interval(cytosol_recovery_fraction, "cytosol_recovery_fraction");
    if (!is_finite_float64(consumed_prey_sentinel) || consumed_prey_sentinel >= 0.0) {
        throw py::value_error("consumed_prey_sentinel must be finite and negative");
    }

    const bacsandbox::physics::AgentBuffers agents{
        {x_data, agent_count},
        {y_data, agent_count},
        {interaction_radius_data, agent_count},
        {base_radius_data, agent_count},
        {food_level_data, agent_count},
        {food_capacity_data, agent_count},
        {metabolic_burden_data, agent_count},
        {species_id_data, agent_count},
        {is_persister_data, agent_count},
        {phage_resistant_data, agent_count},
        {event_flags_data, agent_count},
    };
    const bacsandbox::physics::SpeciesBuffers species{
        {t6ss_activity_data, species_count},
        {immunity_factor_data, species_count},
        {hgt_competence_data, species_count},
        {share_threshold_data, species_count},
        {predation_threshold_data, species_count},
    };
    const bacsandbox::physics::InteractionParameters parameters{
        dt,
        hash_cell_size,
        axis_cutoff,
        contact_margin,
        hgt_transfer_rate,
        hgt_metabolic_burden,
        t6ss_damage_scaling,
        t6ss_max_immunity_defense,
        t6ss_immunity_size_epsilon,
        persister_resistance_factor,
        phagocytosis_size_ratio,
        cell_wet_mass,
        cell_dry_mass_fraction,
        labile_necromass_fraction,
        cytosol_recovery_fraction,
        consumed_prey_sentinel,
    };

    bacsandbox::physics::ContactPairs contacts;
    {
        py::gil_scoped_release release;
        contacts = bacsandbox::physics::discover_contact_pairs(agents, parameters);
    }

    // Resolution remains serial with the GIL held because only eligible HGT contacts invoke Python RNG;
    // see DOCUMENTATION.md §4.18 and §10.2.
    const bacsandbox::physics::RandomDraw native_random_draw = [&random_draw]() {
        return random_draw().cast<double>();
    };
    // Caller-owned arrays commit atomically: snapshots are restored only when callback validation
    // or interaction resolution fails; see DOCUMENTATION.md §4.18 and §10.2.
    const std::vector<double> original_food_level(
        food_level_data,
        food_level_data + agent_count);
    const std::vector<double> original_metabolic_burden(
        metabolic_burden_data,
        metabolic_burden_data + agent_count);
    const std::vector<std::uint8_t> original_phage_resistant(
        phage_resistant_data,
        phage_resistant_data + agent_count);
    const std::vector<std::uint8_t> original_event_flags(
        event_flags_data,
        event_flags_data + agent_count);

    bacsandbox::physics::InteractionSummary summary;
    try {
        summary = bacsandbox::physics::resolve_interactions(
            agents,
            species,
            parameters,
            contacts,
            native_random_draw);
    } catch (...) {
        std::copy(original_food_level.begin(), original_food_level.end(), food_level_data);
        std::copy(
            original_metabolic_burden.begin(),
            original_metabolic_burden.end(),
            metabolic_burden_data);
        std::copy(
            original_phage_resistant.begin(),
            original_phage_resistant.end(),
            phage_resistant_data);
        std::copy(original_event_flags.begin(), original_event_flags.end(), event_flags_data);
        throw;
    }
    return py::make_tuple(
        summary.contact_pairs,
        summary.t6ss_hits,
        summary.t6ss_deaths,
        summary.hgt_transfers,
        summary.share_events,
        summary.predation_events);
}

}  // namespace

PYBIND11_MODULE(bacteria_c_kernel, module) {
    module.doc() = "Native BacSandbox computational kernels";
    module.def(
        "process_layers",
        &process_layers_binding,
        py::arg("coordinates"),
        py::arg("values"),
        py::arg("offsets"),
        py::arg("diffusion"),
        py::arg("decay"),
        py::arg("dt"),
        py::arg("cell_size"),
        py::arg("alpha_max"),
        py::arg("min_concentration"));
    module.def(
        "process_interactions",
        &process_interactions_binding,
        py::arg("x").noconvert(),
        py::arg("y").noconvert(),
        py::arg("interaction_radius").noconvert(),
        py::arg("base_radius").noconvert(),
        py::arg("food_level").noconvert(),
        py::arg("food_capacity").noconvert(),
        py::arg("metabolic_burden_factor").noconvert(),
        py::arg("species_id").noconvert(),
        py::arg("is_persister").noconvert(),
        py::arg("phage_resistant").noconvert(),
        py::arg("event_flags").noconvert(),
        py::arg("t6ss_activity").noconvert(),
        py::arg("immunity_factor").noconvert(),
        py::arg("hgt_competence").noconvert(),
        py::arg("share_threshold").noconvert(),
        py::arg("predation_threshold").noconvert(),
        py::arg("dt"),
        py::arg("hash_cell_size"),
        py::arg("axis_cutoff"),
        py::arg("contact_margin"),
        py::arg("hgt_transfer_rate"),
        py::arg("hgt_metabolic_burden"),
        py::arg("t6ss_damage_scaling"),
        py::arg("t6ss_max_immunity_defense"),
        py::arg("t6ss_immunity_size_epsilon"),
        py::arg("persister_resistance_factor"),
        py::arg("phagocytosis_size_ratio"),
        py::arg("cell_wet_mass"),
        py::arg("cell_dry_mass_fraction"),
        py::arg("labile_necromass_fraction"),
        py::arg("cytosol_recovery_fraction"),
        py::arg("consumed_prey_sentinel"),
        py::arg("random_draw"));
    module.def("openmp_enabled", &bacsandbox::ecm::openmp_enabled);
}