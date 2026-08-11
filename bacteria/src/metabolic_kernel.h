#ifndef BACSANDBOX_METABOLIC_KERNEL_H
#define BACSANDBOX_METABOLIC_KERNEL_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace bacsandbox::metabolism {

// Three substrate slots preserve the reviewed JSON/NumPy ABI order; see DOCUMENTATION.md §10.3.
constexpr std::size_t substrate_count = 3;
constexpr std::size_t maximum_cache_capacity = 1000000;
constexpr double simplex_tolerance = 1.0e-10;

enum class LpStatus : std::uint8_t {
    optimal,
    infeasible,
    unbounded,
    numerical_error,
};

struct NativeModelData {
    std::vector<std::int64_t> model_reaction_offsets;
    std::vector<std::int64_t> model_metabolite_offsets;
    std::vector<std::int64_t> stoich_row_offsets;
    std::vector<std::int32_t> stoich_column_indices;
    std::vector<double> stoich_values;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> objective_coefficients;
    std::vector<std::int32_t> exchange_reaction_indices;
    std::vector<std::int32_t> biomass_reaction_indices;
    std::vector<std::int32_t> acetate_reaction_indices;
    double cache_tolerance_molar = 0.0;
    std::size_t cache_capacity = 0;
};

struct NativePhenotypeData {
    std::vector<std::int32_t> model_indices;
    std::vector<double> uptake_vmax;
    std::vector<double> uptake_km;
    std::vector<double> inhibition_multiplier;
    std::uint64_t model_revision = 0;
    std::uint64_t prior_revision = 0;
};

struct DenseModel {
    std::size_t reaction_start = 0;
    std::size_t reaction_count = 0;
    std::size_t metabolite_count = 0;
    std::vector<double> stoichiometry;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> objective;
    std::array<std::size_t, substrate_count> exchange_local_indices{};
    std::array<std::int32_t, substrate_count> exchange_global_indices{};
    std::size_t biomass_local_index = 0;
    std::size_t acetate_local_index = 0;
    std::int32_t biomass_global_index = 0;
    std::int32_t acetate_global_index = 0;
};

struct CacheKey {
    std::uint64_t model_revision = 0;
    std::uint64_t prior_revision = 0;
    std::int32_t phenotype_index = 0;
    std::array<std::int64_t, substrate_count> concentration_bins{};
};

struct CacheKeyHash {
    std::size_t operator()(const CacheKey& key) const noexcept;
};

struct CacheKeyEqual {
    bool operator()(const CacheKey& left, const CacheKey& right) const noexcept;
};

struct CacheEntry {
    std::array<double, substrate_count> concentrations{};
    std::array<double, substrate_count> effective_exchange_bounds{};
    std::vector<double> flux_values;
    std::uint64_t insertion_sequence = 0;
    std::uint64_t last_used_sequence = 0;
};

struct CacheStats {
    std::size_t entries = 0;
    std::uint64_t hits = 0;
    std::uint64_t misses = 0;
    std::uint64_t evictions = 0;
    std::uint64_t lp_solves = 0;
};

struct CacheSnapshot {
    std::unordered_map<CacheKey, CacheEntry, CacheKeyHash, CacheKeyEqual> entries;
    CacheStats stats;
    std::uint64_t sequence = 0;
};

struct FluxResult {
    std::vector<double> flux_values;
    bool cache_hit = false;
};

struct BatchSummary {
    std::size_t substeps = 0;
    std::uint64_t lp_solves = 0;
    std::uint64_t cache_hits = 0;
    std::uint64_t cache_misses = 0;
    std::uint64_t cache_evictions = 0;
    double minimum_substep_seconds = 0.0;
};

struct BatchBuffers {
    const std::int32_t* phenotype_id = nullptr;
    const std::int32_t* environment_id = nullptr;
    std::size_t agent_count = 0;
    double* concentrations_m = nullptr;
    const double* environment_volume_l = nullptr;
    std::size_t environment_count = 0;
    double* biomass_gdw = nullptr;
    double* growth_rate_s = nullptr;
    double* biomass_delta_gdw = nullptr;
    double* substrate_consumed_mol = nullptr;
    double* acetate_secreted_mol = nullptr;
};

struct PreparedSolve {
    std::int32_t phenotype_index = 0;
    std::size_t model_index = 0;
    std::array<double, substrate_count> concentrations{};
    std::array<double, substrate_count> effective_exchange_bounds{};
    std::array<std::int64_t, substrate_count> concentration_bins{};
};

struct SimplexWorkspace {
    std::vector<double> constraint_matrix;
    std::vector<double> constraint_bounds;
    std::vector<double> tableau;
    std::vector<double> shifted_solution;
    std::vector<std::int64_t> basis;
    std::vector<std::int64_t> nonbasis;

    void prepare(std::size_t equality_rows, std::size_t reaction_columns);
};

class MetabolicKernel {
public:
    MetabolicKernel(
        std::vector<std::int64_t> model_reaction_offsets,
        std::vector<std::int64_t> model_metabolite_offsets,
        std::vector<std::int64_t> stoich_row_offsets,
        std::vector<std::int32_t> stoich_column_indices,
        std::vector<double> stoich_values,
        std::vector<double> lower_bounds,
        std::vector<double> upper_bounds,
        std::vector<double> objective_coefficients,
        std::vector<std::int32_t> exchange_reaction_indices,
        std::vector<std::int32_t> biomass_reaction_indices,
        std::vector<std::int32_t> acetate_reaction_indices,
        double cache_tolerance_molar,
        std::size_t cache_capacity);

    void configure_phenotypes(
        std::vector<std::int32_t> model_indices,
        std::vector<double> uptake_vmax,
        std::vector<double> uptake_km,
        std::vector<double> inhibition_alpha,
        std::vector<double> inhibition_beta,
        std::vector<double> max_inhibition,
        std::uint64_t model_revision,
        std::uint64_t prior_revision);

    FluxResult solve_flux(
        std::int32_t phenotype_index,
        const std::array<double, substrate_count>& concentrations_m,
        bool cache_enabled = true);

    BatchSummary advance_batch(
        const BatchBuffers& buffers,
        double dt_seconds,
        double low_substrate_threshold_m,
        double depletion_safety,
        std::size_t max_substeps,
        bool cache_enabled);

    void clear_cache() noexcept;
    CacheStats cache_stats() const noexcept;
    CacheSnapshot snapshot_cache() const;
    void restore_cache(const CacheSnapshot& snapshot);

private:
    PreparedSolve prepare_solve(
        std::int32_t phenotype_index,
        const std::array<double, substrate_count>& concentrations_m) const;
    bool lookup_cache(const PreparedSolve& request, double* flux_output);
    void insert_cache(const PreparedSolve& request, const double* flux_values);
    LpStatus solve_prepared(
        const PreparedSolve& request,
        SimplexWorkspace& workspace,
        double* flux_output) const noexcept;
    bool requests_match(
        const PreparedSolve& left,
        const PreparedSolve& right) const noexcept;

    NativeModelData model_;
    NativePhenotypeData phenotypes_;
    std::vector<DenseModel> dense_models_;
    std::unordered_map<CacheKey, CacheEntry, CacheKeyHash, CacheKeyEqual> cache_;
    CacheStats cache_stats_;
    std::uint64_t sequence_ = 0;
};

}  // namespace bacsandbox::metabolism

#endif