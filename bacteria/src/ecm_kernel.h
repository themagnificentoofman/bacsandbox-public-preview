#ifndef BACSANDBOX_ECM_KERNEL_H
#define BACSANDBOX_ECM_KERNEL_H

#include <array>
#include <cstdint>
#include <vector>

namespace bacsandbox::ecm {

struct Coordinate {
    std::int64_t x;
    std::int64_t y;

    bool operator==(const Coordinate& other) const noexcept {
        return x == other.x && y == other.y;
    }
};

struct LayerBatch {
    std::vector<Coordinate> coordinates;
    std::vector<double> values;
    // Five offsets partition the fixed qs/bac/acetate/phage slices without duplicating model constants; see DOCUMENTATION.md §4.9.
    std::array<std::int64_t, 5> offsets;
};

struct KernelParameters {
    std::array<double, 4> diffusion;
    std::array<double, 4> decay;
    double dt;
    double cell_size;
    double alpha_max;
    double min_concentration;
};

struct LayerResult {
    std::vector<Coordinate> coordinates;
    std::vector<double> values;
    std::array<std::int64_t, 5> offsets;
};

LayerResult dispatch_process_layers(const LayerBatch& layers, const KernelParameters& parameters);
bool openmp_enabled() noexcept;

}  // namespace bacsandbox::ecm

#endif