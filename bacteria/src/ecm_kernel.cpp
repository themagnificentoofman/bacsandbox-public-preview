#include "ecm_kernel.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <exception>
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

}  // namespace

PYBIND11_MODULE(bacteria_c_kernel, module) {
    module.doc() = "Native BacSandbox ECM kernel boundary";
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
    module.def("openmp_enabled", &bacsandbox::ecm::openmp_enabled);
}