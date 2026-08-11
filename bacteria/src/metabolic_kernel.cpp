#include "metabolic_kernel.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#ifndef BACSANDBOX_OPENMP_ENABLED
#define BACSANDBOX_OPENMP_ENABLED 0
#endif

namespace bacsandbox::metabolism {

namespace {

bool is_finite_float64(double value) noexcept {
    static_assert(sizeof(double) == sizeof(std::uint64_t), "metabolic buffers require IEEE-754 float64");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    constexpr std::uint64_t exponent_mask = 0x7ff0000000000000ULL;
    return (bits & exponent_mask) != exponent_mask;
}

std::size_t checked_product(std::size_t left, std::size_t right, const char* name) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::overflow_error(std::string(name) + " exceeds the native size range");
    }
    return left * right;
}

std::size_t checked_extent(std::int64_t value, const char* name) {
    if (value < 0 || static_cast<std::uintmax_t>(value) > std::numeric_limits<std::size_t>::max()) {
        throw std::overflow_error(std::string(name) + " exceeds the native size range");
    }
    return static_cast<std::size_t>(value);
}

void require_finite(double value, const char* name) {
    if (!is_finite_float64(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

void require_finite_nonnegative(double value, const char* name) {
    if (!is_finite_float64(value) || value < 0.0) {
        throw std::invalid_argument(std::string(name) + " must be finite and nonnegative");
    }
}

void require_finite_positive(double value, const char* name) {
    if (!is_finite_float64(value) || value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be finite and positive");
    }
}

void validate_offsets(
    const std::vector<std::int64_t>& offsets,
    std::size_t expected_final,
    const char* name) {
    if (offsets.size() < 2 || offsets.front() != 0) {
        throw std::invalid_argument(std::string(name) + " must start at zero and contain at least two entries");
    }
    for (std::size_t index = 0; index + 1 < offsets.size(); ++index) {
        if (offsets[index] < 0 || offsets[index] > offsets[index + 1]) {
            throw std::invalid_argument(std::string(name) + " must be nonnegative and monotonic");
        }
    }
    if (expected_final > static_cast<std::uintmax_t>(std::numeric_limits<std::int64_t>::max()) ||
        offsets.back() != static_cast<std::int64_t>(expected_final)) {
        throw std::invalid_argument(std::string(name) + " must span its flattened values");
    }
}

void validate_model(const NativeModelData& model) {
    const std::size_t reaction_count = model.lower_bounds.size();
    if (reaction_count == 0 || model.upper_bounds.size() != reaction_count ||
        model.objective_coefficients.size() != reaction_count) {
        throw std::invalid_argument("reaction bounds and objective must have one nonempty shared extent");
    }
    if (reaction_count > static_cast<std::uintmax_t>(std::numeric_limits<std::int32_t>::max())) {
        throw std::overflow_error("reaction count exceeds the int32 stoichiometry-index ABI");
    }
    validate_offsets(model.model_reaction_offsets, reaction_count, "model_reaction_offsets");
    const std::size_t model_count = model.model_reaction_offsets.size() - 1;
    if (model.model_metabolite_offsets.size() != model_count + 1) {
        throw std::invalid_argument("model_metabolite_offsets must have shape (M+1,)");
    }
    const std::size_t metabolite_count = checked_extent(
        model.model_metabolite_offsets.back(), "model_metabolite_offsets");
    if (metabolite_count == 0) {
        throw std::invalid_argument("models must contain at least one metabolite");
    }
    validate_offsets(model.model_metabolite_offsets, metabolite_count, "model_metabolite_offsets");
    if (model.stoich_row_offsets.size() != metabolite_count + 1 ||
        model.stoich_column_indices.size() != model.stoich_values.size()) {
        throw std::invalid_argument("stoichiometry arrays have inconsistent extents");
    }
    validate_offsets(model.stoich_row_offsets, model.stoich_values.size(), "stoich_row_offsets");
    for (std::size_t row = 0; row < metabolite_count; ++row) {
        if (model.stoich_row_offsets[row] == model.stoich_row_offsets[row + 1]) {
            throw std::invalid_argument("every metabolite row must contain stoichiometry");
        }
    }
    for (std::size_t reaction = 0; reaction < reaction_count; ++reaction) {
        require_finite(model.lower_bounds[reaction], "lower_bounds entries");
        require_finite(model.upper_bounds[reaction], "upper_bounds entries");
        require_finite(model.objective_coefficients[reaction], "objective_coefficients entries");
        if (model.lower_bounds[reaction] > model.upper_bounds[reaction]) {
            throw std::invalid_argument("reaction lower bounds must not exceed upper bounds");
        }
    }
    for (double coefficient : model.stoich_values) {
        require_finite(coefficient, "stoich_values entries");
        if (coefficient == 0.0) {
            throw std::invalid_argument("stoich_values entries must be nonzero");
        }
    }
    if (checked_product(model_count, substrate_count, "exchange shape") !=
        model.exchange_reaction_indices.size()) {
        throw std::invalid_argument("exchange_reaction_indices must have shape (M, 3)");
    }
    if (model.biomass_reaction_indices.size() != model_count ||
        model.acetate_reaction_indices.size() != model_count) {
        throw std::invalid_argument("biomass and acetate indices must have shape (M,)");
    }
    for (std::size_t model_index = 0; model_index < model_count; ++model_index) {
        const std::int64_t reaction_start = model.model_reaction_offsets[model_index];
        const std::int64_t reaction_stop = model.model_reaction_offsets[model_index + 1];
        const std::int64_t metabolite_start = model.model_metabolite_offsets[model_index];
        const std::int64_t metabolite_stop = model.model_metabolite_offsets[model_index + 1];
        if (reaction_start == reaction_stop || metabolite_start == metabolite_stop) {
            throw std::invalid_argument("every model must contain reactions and metabolites");
        }
        const auto belongs = [reaction_start, reaction_stop](std::int32_t index) {
            return index >= reaction_start && index < reaction_stop;
        };
        const std::int32_t biomass = model.biomass_reaction_indices[model_index];
        const std::int32_t acetate = model.acetate_reaction_indices[model_index];
        if (!belongs(biomass) || !belongs(acetate) || biomass == acetate ||
            model.objective_coefficients[static_cast<std::size_t>(biomass)] <= 0.0) {
            throw std::invalid_argument("biomass and acetate markers are invalid");
        }
        for (std::int64_t reaction = reaction_start; reaction < reaction_stop; ++reaction) {
            if (reaction != biomass &&
                model.objective_coefficients[static_cast<std::size_t>(reaction)] != 0.0) {
                throw std::invalid_argument("only each model's biomass reaction may carry objective weight");
            }
        }
        std::array<std::int32_t, substrate_count> exchange{};
        for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
            exchange[substrate] =
                model.exchange_reaction_indices[model_index * substrate_count + substrate];
            if (!belongs(exchange[substrate]) || exchange[substrate] == biomass ||
                exchange[substrate] == acetate) {
                throw std::invalid_argument("exchange indices must select non-marker reactions in their model");
            }
        }
        if (exchange[0] == exchange[1] || exchange[0] == exchange[2] ||
            exchange[1] == exchange[2]) {
            throw std::invalid_argument("each model's exchange indices must be unique");
        }
        for (std::int64_t row = metabolite_start; row < metabolite_stop; ++row) {
            const std::size_t begin = checked_extent(
                model.stoich_row_offsets[static_cast<std::size_t>(row)], "stoich row offset");
            const std::size_t end = checked_extent(
                model.stoich_row_offsets[static_cast<std::size_t>(row) + 1], "stoich row offset");
            for (std::size_t entry = begin; entry < end; ++entry) {
                if (!belongs(model.stoich_column_indices[entry])) {
                    throw std::invalid_argument("stoichiometry columns must remain within each model");
                }
            }
        }
    }
    require_finite_positive(model.cache_tolerance_molar, "cache_tolerance_molar");
    if (model.cache_capacity == 0 || model.cache_capacity > maximum_cache_capacity) {
        throw std::invalid_argument("cache_capacity is outside the supported range");
    }
}

std::vector<DenseModel> make_dense_models(const NativeModelData& model) {
    const std::size_t model_count = model.model_reaction_offsets.size() - 1;
    std::vector<DenseModel> result;
    result.reserve(model_count);
    for (std::size_t model_index = 0; model_index < model_count; ++model_index) {
        DenseModel dense;
        dense.reaction_start = checked_extent(
            model.model_reaction_offsets[model_index], "model reaction offset");
        const std::size_t reaction_stop = checked_extent(
            model.model_reaction_offsets[model_index + 1], "model reaction offset");
        dense.reaction_count = reaction_stop - dense.reaction_start;
        const std::size_t metabolite_start = checked_extent(
            model.model_metabolite_offsets[model_index], "model metabolite offset");
        const std::size_t metabolite_stop = checked_extent(
            model.model_metabolite_offsets[model_index + 1], "model metabolite offset");
        dense.metabolite_count = metabolite_stop - metabolite_start;
        dense.stoichiometry.assign(
            checked_product(dense.metabolite_count, dense.reaction_count, "dense stoichiometry"),
            0.0);
        for (std::size_t global_row = metabolite_start; global_row < metabolite_stop; ++global_row) {
            const std::size_t local_row = global_row - metabolite_start;
            const std::size_t begin = checked_extent(model.stoich_row_offsets[global_row], "stoich row offset");
            const std::size_t end = checked_extent(model.stoich_row_offsets[global_row + 1], "stoich row offset");
            for (std::size_t entry = begin; entry < end; ++entry) {
                const std::size_t global_column =
                    static_cast<std::size_t>(model.stoich_column_indices[entry]);
                dense.stoichiometry[local_row * dense.reaction_count +
                                    global_column - dense.reaction_start] = model.stoich_values[entry];
            }
        }
        dense.lower_bounds.assign(
            model.lower_bounds.begin() + static_cast<std::ptrdiff_t>(dense.reaction_start),
            model.lower_bounds.begin() + static_cast<std::ptrdiff_t>(reaction_stop));
        dense.upper_bounds.assign(
            model.upper_bounds.begin() + static_cast<std::ptrdiff_t>(dense.reaction_start),
            model.upper_bounds.begin() + static_cast<std::ptrdiff_t>(reaction_stop));
        dense.objective.assign(
            model.objective_coefficients.begin() + static_cast<std::ptrdiff_t>(dense.reaction_start),
            model.objective_coefficients.begin() + static_cast<std::ptrdiff_t>(reaction_stop));
        for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
            const std::int32_t global =
                model.exchange_reaction_indices[model_index * substrate_count + substrate];
            dense.exchange_global_indices[substrate] = global;
            dense.exchange_local_indices[substrate] =
                static_cast<std::size_t>(global) - dense.reaction_start;
        }
        dense.biomass_global_index = model.biomass_reaction_indices[model_index];
        dense.acetate_global_index = model.acetate_reaction_indices[model_index];
        dense.biomass_local_index =
            static_cast<std::size_t>(dense.biomass_global_index) - dense.reaction_start;
        dense.acetate_local_index =
            static_cast<std::size_t>(dense.acetate_global_index) - dense.reaction_start;
        result.push_back(std::move(dense));
    }
    return result;
}

double& tableau_value(
    SimplexWorkspace& workspace,
    std::size_t tableau_columns,
    std::size_t row,
    std::size_t column) noexcept {
    return workspace.tableau[row * tableau_columns + column];
}

double tableau_value(
    const SimplexWorkspace& workspace,
    std::size_t tableau_columns,
    std::size_t row,
    std::size_t column) noexcept {
    return workspace.tableau[row * tableau_columns + column];
}

void pivot(
    SimplexWorkspace& workspace,
    std::size_t row_count,
    std::size_t column_count,
    std::size_t pivot_row,
    std::size_t pivot_column) noexcept {
    const std::size_t tableau_columns = column_count + 2;
    const double inverse = 1.0 / tableau_value(
        workspace, tableau_columns, pivot_row, pivot_column);
    // Gauss-Jordan elimination preserves the basic dictionary, while stable indices
    // make degenerate pivots deterministic; see DOCUMENTATION.md §10.3.
    for (std::size_t row = 0; row < row_count + 2; ++row) {
        if (row == pivot_row) {
            continue;
        }
        for (std::size_t column = 0; column < column_count + 2; ++column) {
            if (column == pivot_column) {
                continue;
            }
            tableau_value(workspace, tableau_columns, row, column) -=
                tableau_value(workspace, tableau_columns, pivot_row, column) *
                tableau_value(workspace, tableau_columns, row, pivot_column) * inverse;
        }
    }
    for (std::size_t column = 0; column < column_count + 2; ++column) {
        if (column != pivot_column) {
            tableau_value(workspace, tableau_columns, pivot_row, column) *= inverse;
        }
    }
    for (std::size_t row = 0; row < row_count + 2; ++row) {
        if (row != pivot_row) {
            tableau_value(workspace, tableau_columns, row, pivot_column) *= -inverse;
        }
    }
    tableau_value(workspace, tableau_columns, pivot_row, pivot_column) = inverse;
    std::swap(workspace.basis[pivot_row], workspace.nonbasis[pivot_column]);
}

bool simplex_phase(
    SimplexWorkspace& workspace,
    std::size_t row_count,
    std::size_t column_count,
    int phase) noexcept {
    const std::size_t tableau_columns = column_count + 2;
    const std::size_t objective_row = phase == 1 ? row_count + 1 : row_count;
    while (true) {
        std::size_t entering = column_count + 1;
        for (std::size_t column = 0; column < column_count + 1; ++column) {
            if (phase == 2 && workspace.nonbasis[column] == -1) {
                continue;
            }
            if (entering == column_count + 1 ||
                tableau_value(workspace, tableau_columns, objective_row, column) <
                    tableau_value(workspace, tableau_columns, objective_row, entering) - simplex_tolerance ||
                (std::abs(
                     tableau_value(workspace, tableau_columns, objective_row, column) -
                     tableau_value(workspace, tableau_columns, objective_row, entering)) <= simplex_tolerance &&
                 workspace.nonbasis[column] < workspace.nonbasis[entering])) {
                entering = column;
            }
        }
        if (entering == column_count + 1 ||
            tableau_value(workspace, tableau_columns, objective_row, entering) >= -simplex_tolerance) {
            return true;
        }
        std::size_t leaving = row_count;
        for (std::size_t row = 0; row < row_count; ++row) {
            const double coefficient = tableau_value(
                workspace, tableau_columns, row, entering);
            if (coefficient <= simplex_tolerance) {
                continue;
            }
            if (leaving == row_count) {
                leaving = row;
                continue;
            }
            const double candidate_ratio =
                tableau_value(workspace, tableau_columns, row, column_count + 1) / coefficient;
            const double current_ratio =
                tableau_value(workspace, tableau_columns, leaving, column_count + 1) /
                tableau_value(workspace, tableau_columns, leaving, entering);
            if (candidate_ratio < current_ratio - simplex_tolerance ||
                (std::abs(candidate_ratio - current_ratio) <= simplex_tolerance &&
                 workspace.basis[row] < workspace.basis[leaving])) {
                leaving = row;
            }
        }
        if (leaving == row_count) {
            return false;
        }
        pivot(workspace, row_count, column_count, leaving, entering);
    }
}

LpStatus solve_nonnegative_lp(
    SimplexWorkspace& workspace,
    std::size_t row_count,
    std::size_t column_count,
    const double* objective) noexcept {
    const std::size_t tableau_columns = column_count + 2;
    std::fill(workspace.tableau.begin(), workspace.tableau.end(), 0.0);
    for (std::size_t row = 0; row < row_count; ++row) {
        for (std::size_t column = 0; column < column_count; ++column) {
            tableau_value(workspace, tableau_columns, row, column) =
                workspace.constraint_matrix[row * column_count + column];
        }
        workspace.basis[row] = static_cast<std::int64_t>(column_count + row);
        tableau_value(workspace, tableau_columns, row, column_count) = -1.0;
        tableau_value(workspace, tableau_columns, row, column_count + 1) =
            workspace.constraint_bounds[row];
    }
    for (std::size_t column = 0; column < column_count; ++column) {
        workspace.nonbasis[column] = static_cast<std::int64_t>(column);
        tableau_value(workspace, tableau_columns, row_count, column) = -objective[column];
    }
    workspace.nonbasis[column_count] = -1;
    tableau_value(workspace, tableau_columns, row_count + 1, column_count) = 1.0;

    if (row_count != 0) {
        std::size_t initial_row = 0;
        for (std::size_t row = 1; row < row_count; ++row) {
            const double candidate = tableau_value(
                workspace, tableau_columns, row, column_count + 1);
            const double current = tableau_value(
                workspace, tableau_columns, initial_row, column_count + 1);
            if (candidate < current ||
                (candidate == current && workspace.basis[row] < workspace.basis[initial_row])) {
                initial_row = row;
            }
        }
        if (tableau_value(workspace, tableau_columns, initial_row, column_count + 1) <
            -simplex_tolerance) {
            pivot(workspace, row_count, column_count, initial_row, column_count);
            if (!simplex_phase(workspace, row_count, column_count, 1) ||
                std::abs(tableau_value(
                    workspace, tableau_columns, row_count + 1, column_count + 1)) > simplex_tolerance) {
                return LpStatus::infeasible;
            }
            for (std::size_t row = 0; row < row_count; ++row) {
                if (workspace.basis[row] != -1) {
                    continue;
                }
                std::size_t entering = column_count + 1;
                for (std::size_t column = 0; column < column_count + 1; ++column) {
                    if (std::abs(tableau_value(
                            workspace, tableau_columns, row, column)) <= simplex_tolerance) {
                        continue;
                    }
                    if (entering == column_count + 1 ||
                        workspace.nonbasis[column] < workspace.nonbasis[entering]) {
                        entering = column;
                    }
                }
                if (entering != column_count + 1) {
                    pivot(workspace, row_count, column_count, row, entering);
                }
            }
        }
    }
    if (!simplex_phase(workspace, row_count, column_count, 2)) {
        return LpStatus::unbounded;
    }
    std::fill(workspace.shifted_solution.begin(), workspace.shifted_solution.end(), 0.0);
    for (std::size_t row = 0; row < row_count; ++row) {
        if (workspace.basis[row] >= 0 &&
            static_cast<std::size_t>(workspace.basis[row]) < column_count) {
            workspace.shifted_solution[static_cast<std::size_t>(workspace.basis[row])] =
                tableau_value(workspace, tableau_columns, row, column_count + 1);
        }
    }
    return LpStatus::optimal;
}

const char* status_message(LpStatus status) noexcept {
    switch (status) {
        case LpStatus::optimal:
            return "optimal";
        case LpStatus::infeasible:
            return "infeasible";
        case LpStatus::unbounded:
            return "unbounded";
        case LpStatus::numerical_error:
            return "numerical_error";
    }
    return "unknown";
}

}  // namespace

void SimplexWorkspace::prepare(
    std::size_t equality_rows,
    std::size_t reaction_columns) {
    if (equality_rows > (std::numeric_limits<std::size_t>::max() - reaction_columns) / 2) {
        throw std::overflow_error("simplex constraint count exceeds the native size range");
    }
    const std::size_t constraint_rows = 2 * equality_rows + reaction_columns;
    constraint_matrix.resize(checked_product(
        constraint_rows, reaction_columns, "simplex constraint matrix"));
    constraint_bounds.resize(constraint_rows);
    tableau.resize(checked_product(
        constraint_rows + 2, reaction_columns + 2, "simplex tableau"));
    shifted_solution.resize(reaction_columns);
    basis.resize(constraint_rows);
    nonbasis.resize(reaction_columns + 1);
}

std::size_t CacheKeyHash::operator()(const CacheKey& key) const noexcept {
    std::size_t value = std::hash<std::uint64_t>{}(key.model_revision);
    value ^= std::hash<std::uint64_t>{}(key.prior_revision) + (value << 1U);
    value ^= std::hash<std::int32_t>{}(key.phenotype_index) + (value << 1U);
    for (std::int64_t bin : key.concentration_bins) {
        value ^= std::hash<std::int64_t>{}(bin) + (value << 1U);
    }
    return value;
}

bool CacheKeyEqual::operator()(const CacheKey& left, const CacheKey& right) const noexcept {
    return left.model_revision == right.model_revision &&
        left.prior_revision == right.prior_revision &&
        left.phenotype_index == right.phenotype_index &&
        left.concentration_bins == right.concentration_bins;
}

MetabolicKernel::MetabolicKernel(
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
    std::size_t cache_capacity)
    : model_{
          std::move(model_reaction_offsets),
          std::move(model_metabolite_offsets),
          std::move(stoich_row_offsets),
          std::move(stoich_column_indices),
          std::move(stoich_values),
          std::move(lower_bounds),
          std::move(upper_bounds),
          std::move(objective_coefficients),
          std::move(exchange_reaction_indices),
          std::move(biomass_reaction_indices),
          std::move(acetate_reaction_indices),
          cache_tolerance_molar,
          cache_capacity} {
    validate_model(model_);
    dense_models_ = make_dense_models(model_);
    cache_.reserve(model_.cache_capacity);
}

void MetabolicKernel::configure_phenotypes(
    std::vector<std::int32_t> model_indices,
    std::vector<double> uptake_vmax,
    std::vector<double> uptake_km,
    std::vector<double> inhibition_alpha,
    std::vector<double> inhibition_beta,
    std::vector<double> max_inhibition,
    std::uint64_t model_revision,
    std::uint64_t prior_revision) {
    const std::size_t phenotype_count = model_indices.size();
    if (phenotype_count == 0) {
        throw std::invalid_argument("configure_phenotypes requires at least one phenotype");
    }
    const std::size_t parameter_count = checked_product(
        phenotype_count, substrate_count, "phenotype parameter count");
    if (uptake_vmax.size() != parameter_count || uptake_km.size() != parameter_count ||
        inhibition_alpha.size() != parameter_count || inhibition_beta.size() != parameter_count ||
        max_inhibition.size() != parameter_count) {
        throw std::invalid_argument("phenotype parameter arrays must have shape (P, 3)");
    }
    for (std::int32_t model_index : model_indices) {
        if (model_index < 0 || static_cast<std::size_t>(model_index) >= dense_models_.size()) {
            throw std::invalid_argument("phenotype model indices must select a configured model");
        }
    }
    std::vector<double> multipliers(parameter_count);
    for (std::size_t index = 0; index < parameter_count; ++index) {
        require_finite_nonnegative(uptake_vmax[index], "uptake_vmax entries");
        require_finite_positive(uptake_km[index], "uptake_km entries");
        require_finite_positive(inhibition_alpha[index], "inhibition_alpha entries");
        require_finite_positive(inhibition_beta[index], "inhibition_beta entries");
        if (!is_finite_float64(max_inhibition[index]) || max_inhibition[index] < 0.0 ||
            max_inhibition[index] > 1.0) {
            throw std::invalid_argument("max_inhibition entries must be finite and in [0, 1]");
        }
        const double posterior_mass = inhibition_alpha[index] + inhibition_beta[index];
        require_finite_positive(posterior_mass, "alpha + beta entries");
        // The posterior mean active fraction applies deterministic inhibition without
        // sampling or hidden biological constants; see DOCUMENTATION.md §10.3.
        multipliers[index] = 1.0 -
            inhibition_alpha[index] / posterior_mass * max_inhibition[index];
    }
    phenotypes_ = {
        std::move(model_indices),
        std::move(uptake_vmax),
        std::move(uptake_km),
        std::move(multipliers),
        model_revision,
        prior_revision,
    };
    clear_cache();
}

void MetabolicKernel::clear_cache() noexcept {
    cache_.clear();
    cache_stats_ = {};
    sequence_ = 0;
}

CacheStats MetabolicKernel::cache_stats() const noexcept {
    CacheStats result = cache_stats_;
    result.entries = cache_.size();
    return result;
}

CacheSnapshot MetabolicKernel::snapshot_cache() const {
    return {cache_, cache_stats_, sequence_};
}

void MetabolicKernel::restore_cache(const CacheSnapshot& snapshot) {
    if (snapshot.entries.size() > model_.cache_capacity) {
        throw std::invalid_argument("cache snapshot exceeds configured capacity");
    }
    cache_ = snapshot.entries;
    cache_stats_ = snapshot.stats;
    sequence_ = snapshot.sequence;
}

PreparedSolve MetabolicKernel::prepare_solve(
    std::int32_t phenotype_index,
    const std::array<double, substrate_count>& concentrations_m) const {
    if (phenotype_index < 0 ||
        static_cast<std::size_t>(phenotype_index) >= phenotypes_.model_indices.size()) {
        throw std::invalid_argument("phenotype_index must select a configured phenotype");
    }
    PreparedSolve request;
    request.phenotype_index = phenotype_index;
    request.model_index = static_cast<std::size_t>(
        phenotypes_.model_indices[static_cast<std::size_t>(phenotype_index)]);
    request.concentrations = concentrations_m;
    const DenseModel& dense = dense_models_[request.model_index];
    const std::size_t parameter_start =
        static_cast<std::size_t>(phenotype_index) * substrate_count;
    for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
        const double concentration = concentrations_m[substrate];
        require_finite_nonnegative(concentration, "concentrations_m entries");
        const double scaled_bin = concentration / model_.cache_tolerance_molar;
        if (!is_finite_float64(scaled_bin) ||
            scaled_bin > static_cast<double>(std::numeric_limits<std::int64_t>::max() - 1)) {
            throw std::invalid_argument("concentrations_m exceed the int64 cache-bin range");
        }
        request.concentration_bins[substrate] =
            static_cast<std::int64_t>(std::floor(scaled_bin));
        const std::size_t parameter = parameter_start + substrate;
        // Saturating phenotype uptake is capped by the finite model exchange bound;
        // zero concentration therefore yields zero uptake; see DOCUMENTATION.md §10.3.
        const double dynamic_bound =
            phenotypes_.uptake_vmax[parameter] * concentration /
            (phenotypes_.uptake_km[parameter] + concentration) *
            phenotypes_.inhibition_multiplier[parameter];
        request.effective_exchange_bounds[substrate] = std::min(
            dense.upper_bounds[dense.exchange_local_indices[substrate]], dynamic_bound);
    }
    return request;
}

bool MetabolicKernel::lookup_cache(
    const PreparedSolve& request,
    double* flux_output) {
    CacheEntry* selected = nullptr;
    for (std::int64_t first = -1; first <= 1; ++first) {
        for (std::int64_t second = -1; second <= 1; ++second) {
            for (std::int64_t third = -1; third <= 1; ++third) {
                const CacheKey key{
                    phenotypes_.model_revision,
                    phenotypes_.prior_revision,
                    request.phenotype_index,
                    {
                        request.concentration_bins[0] + first,
                        request.concentration_bins[1] + second,
                        request.concentration_bins[2] + third,
                    },
                };
                auto iterator = cache_.find(key);
                if (iterator == cache_.end()) {
                    continue;
                }
                CacheEntry& entry = iterator->second;
                bool matches = true;
                for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
                    if (std::abs(entry.concentrations[substrate] -
                                 request.concentrations[substrate]) > model_.cache_tolerance_molar) {
                        matches = false;
                        break;
                    }
                }
                if (matches &&
                    (selected == nullptr ||
                     entry.insertion_sequence < selected->insertion_sequence)) {
                    selected = &entry;
                }
            }
        }
    }
    if (selected == nullptr) {
        return false;
    }
    ++sequence_;
    selected->last_used_sequence = sequence_;
    std::copy(selected->flux_values.begin(), selected->flux_values.end(), flux_output);
    ++cache_stats_.hits;
    return true;
}

void MetabolicKernel::insert_cache(
    const PreparedSolve& request,
    const double* flux_values) {
    const CacheKey key{
        phenotypes_.model_revision,
        phenotypes_.prior_revision,
        request.phenotype_index,
        request.concentration_bins,
    };
    auto existing = cache_.find(key);
    if (existing == cache_.end() && cache_.size() >= model_.cache_capacity) {
        auto eviction = cache_.end();
        for (auto iterator = cache_.begin(); iterator != cache_.end(); ++iterator) {
            if (eviction == cache_.end() ||
                iterator->second.last_used_sequence < eviction->second.last_used_sequence ||
                (iterator->second.last_used_sequence == eviction->second.last_used_sequence &&
                 iterator->second.insertion_sequence < eviction->second.insertion_sequence)) {
                eviction = iterator;
            }
        }
        cache_.erase(eviction);
        ++cache_stats_.evictions;
    }
    ++sequence_;
    CacheEntry entry;
    entry.concentrations = request.concentrations;
    entry.effective_exchange_bounds = request.effective_exchange_bounds;
    entry.flux_values.assign(flux_values, flux_values + model_.lower_bounds.size());
    entry.insertion_sequence = sequence_;
    entry.last_used_sequence = sequence_;
    cache_[key] = std::move(entry);
}

LpStatus MetabolicKernel::solve_prepared(
    const PreparedSolve& request,
    SimplexWorkspace& workspace,
    double* flux_output) const noexcept {
    std::fill(flux_output, flux_output + model_.lower_bounds.size(), 0.0);
    const DenseModel& dense = dense_models_[request.model_index];
    const std::size_t equality_rows = dense.metabolite_count;
    const std::size_t columns = dense.reaction_count;
    const std::size_t constraint_rows = 2 * equality_rows + columns;
    std::fill(workspace.constraint_matrix.begin(), workspace.constraint_matrix.end(), 0.0);
    for (std::size_t row = 0; row < equality_rows; ++row) {
        double equality_bound = 0.0;
        for (std::size_t column = 0; column < columns; ++column) {
            const double coefficient = dense.stoichiometry[row * columns + column];
            workspace.constraint_matrix[row * columns + column] = coefficient;
            workspace.constraint_matrix[(row + equality_rows) * columns + column] = -coefficient;
            equality_bound -= coefficient * dense.lower_bounds[column];
        }
        workspace.constraint_bounds[row] = equality_bound;
        workspace.constraint_bounds[row + equality_rows] = -equality_bound;
    }
    for (std::size_t column = 0; column < columns; ++column) {
        const std::size_t row = 2 * equality_rows + column;
        workspace.constraint_matrix[row * columns + column] = 1.0;
        double upper = dense.upper_bounds[column];
        for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
            if (dense.exchange_local_indices[substrate] == column) {
                upper = request.effective_exchange_bounds[substrate];
            }
        }
        workspace.constraint_bounds[row] = upper - dense.lower_bounds[column];
    }
    const LpStatus status = solve_nonnegative_lp(
        workspace, constraint_rows, columns, dense.objective.data());
    if (status != LpStatus::optimal) {
        return status;
    }
    for (std::size_t column = 0; column < columns; ++column) {
        double upper = dense.upper_bounds[column];
        for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
            if (dense.exchange_local_indices[substrate] == column) {
                upper = request.effective_exchange_bounds[substrate];
            }
        }
        double value = dense.lower_bounds[column] + workspace.shifted_solution[column];
        if (std::abs(value) <= simplex_tolerance) {
            value = 0.0;
        }
        value = std::min(upper, std::max(dense.lower_bounds[column], value));
        if (!is_finite_float64(value) || value < dense.lower_bounds[column] - 10.0 * simplex_tolerance ||
            value > upper + 10.0 * simplex_tolerance) {
            return LpStatus::numerical_error;
        }
        flux_output[dense.reaction_start + column] = value;
    }
    for (std::size_t row = 0; row < equality_rows; ++row) {
        double residual = 0.0;
        for (std::size_t column = 0; column < columns; ++column) {
            residual += dense.stoichiometry[row * columns + column] *
                flux_output[dense.reaction_start + column];
        }
        if (!is_finite_float64(residual) || std::abs(residual) > 10.0 * simplex_tolerance) {
            return LpStatus::numerical_error;
        }
    }
    return LpStatus::optimal;
}

bool MetabolicKernel::requests_match(
    const PreparedSolve& left,
    const PreparedSolve& right) const noexcept {
    return left.phenotype_index == right.phenotype_index &&
        left.concentrations == right.concentrations &&
        left.effective_exchange_bounds == right.effective_exchange_bounds;
}

FluxResult MetabolicKernel::solve_flux(
    std::int32_t phenotype_index,
    const std::array<double, substrate_count>& concentrations_m,
    bool cache_enabled) {
    const PreparedSolve request = prepare_solve(phenotype_index, concentrations_m);
    FluxResult result;
    result.flux_values.resize(model_.lower_bounds.size());
    if (cache_enabled && lookup_cache(request, result.flux_values.data())) {
        result.cache_hit = true;
        return result;
    }
    if (cache_enabled) {
        ++cache_stats_.misses;
    }
    SimplexWorkspace workspace;
    const DenseModel& dense = dense_models_[request.model_index];
    workspace.prepare(dense.metabolite_count, dense.reaction_count);
    const LpStatus status = solve_prepared(request, workspace, result.flux_values.data());
    ++cache_stats_.lp_solves;
    if (status != LpStatus::optimal) {
        throw std::runtime_error(std::string("metabolic LP status: ") + status_message(status));
    }
    if (cache_enabled) {
        insert_cache(request, result.flux_values.data());
    }
    return result;
}

BatchSummary MetabolicKernel::advance_batch(
    const BatchBuffers& buffers,
    double dt_seconds,
    double low_substrate_threshold_m,
    double depletion_safety,
    std::size_t max_substeps,
    bool cache_enabled) {
    if (phenotypes_.model_indices.empty()) {
        throw std::runtime_error("configure_phenotypes must be called before advancing");
    }
    require_finite_nonnegative(dt_seconds, "dt_seconds");
    require_finite_nonnegative(low_substrate_threshold_m, "low_substrate_threshold_m");
    if (!is_finite_float64(depletion_safety) || depletion_safety <= 0.0 || depletion_safety > 1.0) {
        throw std::invalid_argument("depletion_safety must be finite and in (0, 1]");
    }
    if (max_substeps == 0) {
        throw std::invalid_argument("max_substeps must be positive");
    }
    if (buffers.agent_count > static_cast<std::uintmax_t>(std::numeric_limits<std::int64_t>::max())) {
        throw std::overflow_error("agent count exceeds the OpenMP loop range");
    }
    const std::size_t concentration_count = checked_product(
        buffers.environment_count, substrate_count, "concentration count");
    const std::size_t substrate_output_count = checked_product(
        buffers.agent_count, substrate_count, "substrate output count");
    for (std::size_t environment = 0; environment < buffers.environment_count; ++environment) {
        require_finite_positive(
            buffers.environment_volume_l[environment], "environment_volume_l entries");
        for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
            require_finite_nonnegative(
                buffers.concentrations_m[environment * substrate_count + substrate],
                "concentrations_m entries");
        }
    }
    for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
        if (buffers.phenotype_id[agent] < 0 ||
            static_cast<std::size_t>(buffers.phenotype_id[agent]) >= phenotypes_.model_indices.size()) {
            throw std::invalid_argument("phenotype_id entries must select configured phenotypes");
        }
        if (buffers.environment_id[agent] < 0 ||
            static_cast<std::size_t>(buffers.environment_id[agent]) >= buffers.environment_count) {
            throw std::invalid_argument("environment_id entries must select environments");
        }
        require_finite_nonnegative(buffers.biomass_gdw[agent], "biomass_gdw entries");
        require_finite(buffers.growth_rate_s[agent], "growth_rate_s entries");
        require_finite(buffers.biomass_delta_gdw[agent], "biomass_delta_gdw entries");
        require_finite(buffers.acetate_secreted_mol[agent], "acetate_secreted_mol entries");
        for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
            require_finite(
                buffers.substrate_consumed_mol[agent * substrate_count + substrate],
                "substrate_consumed_mol entries");
        }
    }

    std::vector<double> original_concentrations(concentration_count);
    std::vector<double> original_biomass(buffers.agent_count);
    std::vector<double> original_growth(buffers.agent_count);
    std::vector<double> original_biomass_delta(buffers.agent_count);
    std::vector<double> original_consumed(substrate_output_count);
    std::vector<double> original_acetate(buffers.agent_count);
    if (concentration_count != 0) {
        std::copy_n(buffers.concentrations_m, concentration_count, original_concentrations.begin());
    }
    if (buffers.agent_count != 0) {
        std::copy_n(buffers.biomass_gdw, buffers.agent_count, original_biomass.begin());
        std::copy_n(buffers.growth_rate_s, buffers.agent_count, original_growth.begin());
        std::copy_n(buffers.biomass_delta_gdw, buffers.agent_count, original_biomass_delta.begin());
        std::copy_n(buffers.acetate_secreted_mol, buffers.agent_count, original_acetate.begin());
    }
    if (substrate_output_count != 0) {
        std::copy_n(buffers.substrate_consumed_mol, substrate_output_count, original_consumed.begin());
    }

    const std::size_t reaction_count = model_.lower_bounds.size();
    std::vector<double> fluxes(checked_product(buffers.agent_count, reaction_count, "batch fluxes"));
    std::vector<double> uptake_fluxes(substrate_output_count);
    std::vector<double> growth_fluxes(buffers.agent_count);
    std::vector<double> acetate_fluxes(buffers.agent_count);
    std::vector<double> concentration_rates(concentration_count);
    std::vector<double> step_consumption(substrate_output_count);
    std::vector<double> step_biomass_delta(buffers.agent_count);
    std::vector<double> step_acetate(buffers.agent_count);
    std::vector<double> growth_integral(buffers.agent_count, 0.0);
    std::vector<double> next_concentrations(concentration_count);
    std::vector<PreparedSolve> prepared(buffers.agent_count);
    std::vector<std::size_t> request_for_agent(buffers.agent_count);
    std::vector<std::size_t> unique_agents;
    unique_agents.reserve(buffers.agent_count);
    std::vector<SimplexWorkspace> workspaces(buffers.agent_count);
    std::vector<LpStatus> statuses(buffers.agent_count, LpStatus::optimal);
    auto original_cache = cache_;
    const CacheStats starting_stats = cache_stats();
    const std::uint64_t starting_sequence = sequence_;

    const auto rollback = [&]() noexcept {
        if (concentration_count != 0) {
            std::copy(original_concentrations.begin(), original_concentrations.end(), buffers.concentrations_m);
        }
        if (buffers.agent_count != 0) {
            std::copy(original_biomass.begin(), original_biomass.end(), buffers.biomass_gdw);
            std::copy(original_growth.begin(), original_growth.end(), buffers.growth_rate_s);
            std::copy(original_biomass_delta.begin(), original_biomass_delta.end(), buffers.biomass_delta_gdw);
            std::copy(original_acetate.begin(), original_acetate.end(), buffers.acetate_secreted_mol);
        }
        if (substrate_output_count != 0) {
            std::copy(original_consumed.begin(), original_consumed.end(), buffers.substrate_consumed_mol);
        }
        cache_.swap(original_cache);
        cache_stats_ = starting_stats;
        sequence_ = starting_sequence;
    };

    try {
        std::fill_n(buffers.growth_rate_s, buffers.agent_count, 0.0);
        std::fill_n(buffers.biomass_delta_gdw, buffers.agent_count, 0.0);
        std::fill_n(buffers.substrate_consumed_mol, substrate_output_count, 0.0);
        std::fill_n(buffers.acetate_secreted_mol, buffers.agent_count, 0.0);
        if (dt_seconds == 0.0) {
            return {};
        }

        BatchSummary summary;
        double elapsed = 0.0;
        double minimum_substep = std::numeric_limits<double>::infinity();
        while (elapsed < dt_seconds) {
            if (summary.substeps >= max_substeps) {
                throw std::runtime_error("max_substeps exhausted before dt_seconds completed");
            }
            unique_agents.clear();
            std::fill(request_for_agent.begin(), request_for_agent.end(),
                      std::numeric_limits<std::size_t>::max());
            for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
                const std::size_t environment =
                    static_cast<std::size_t>(buffers.environment_id[agent]);
                const std::array<double, substrate_count> concentrations{
                    buffers.concentrations_m[environment * substrate_count],
                    buffers.concentrations_m[environment * substrate_count + 1],
                    buffers.concentrations_m[environment * substrate_count + 2],
                };
                prepared[agent] = prepare_solve(buffers.phenotype_id[agent], concentrations);
                double* agent_flux = fluxes.data() + agent * reaction_count;
                if (cache_enabled && lookup_cache(prepared[agent], agent_flux)) {
                    continue;
                }
                if (cache_enabled) {
                    ++cache_stats_.misses;
                    for (std::size_t request = 0; request < unique_agents.size(); ++request) {
                        if (requests_match(prepared[agent], prepared[unique_agents[request]])) {
                            request_for_agent[agent] = request;
                            break;
                        }
                    }
                    if (request_for_agent[agent] != std::numeric_limits<std::size_t>::max()) {
                        continue;
                    }
                }
                request_for_agent[agent] = unique_agents.size();
                unique_agents.push_back(agent);
            }
            for (std::size_t request = 0; request < unique_agents.size(); ++request) {
                const DenseModel& dense = dense_models_[prepared[unique_agents[request]].model_index];
                workspaces[request].prepare(dense.metabolite_count, dense.reaction_count);
            }

#if BACSANDBOX_OPENMP_ENABLED
#pragma omp parallel for schedule(static)
#endif
            for (std::int64_t request = 0;
                 request < static_cast<std::int64_t>(unique_agents.size());
                 ++request) {
                const std::size_t request_index = static_cast<std::size_t>(request);
                const std::size_t agent = unique_agents[request_index];
                statuses[request_index] = solve_prepared(
                    prepared[agent],
                    workspaces[request_index],
                    fluxes.data() + agent * reaction_count);
            }
            cache_stats_.lp_solves += static_cast<std::uint64_t>(unique_agents.size());
            for (std::size_t request = 0; request < unique_agents.size(); ++request) {
                if (statuses[request] != LpStatus::optimal) {
                    throw std::runtime_error(
                        std::string("metabolic LP status: ") + status_message(statuses[request]));
                }
            }
            if (cache_enabled) {
                for (std::size_t request = 0; request < unique_agents.size(); ++request) {
                    const std::size_t agent = unique_agents[request];
                    insert_cache(prepared[agent], fluxes.data() + agent * reaction_count);
                }
            }
            for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
                const std::size_t request = request_for_agent[agent];
                if (request == std::numeric_limits<std::size_t>::max()) {
                    continue;
                }
                const std::size_t source_agent = unique_agents[request];
                if (source_agent != agent) {
                    std::copy_n(
                        fluxes.data() + source_agent * reaction_count,
                        reaction_count,
                        fluxes.data() + agent * reaction_count);
                }
            }

            std::fill(concentration_rates.begin(), concentration_rates.end(), 0.0);
            for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
                const std::size_t phenotype = static_cast<std::size_t>(buffers.phenotype_id[agent]);
                const DenseModel& dense = dense_models_[
                    static_cast<std::size_t>(phenotypes_.model_indices[phenotype])];
                const std::size_t environment =
                    static_cast<std::size_t>(buffers.environment_id[agent]);
                for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
                    const double uptake = fluxes[
                        agent * reaction_count +
                        static_cast<std::size_t>(dense.exchange_global_indices[substrate])];
                    uptake_fluxes[agent * substrate_count + substrate] = uptake;
                    // Specific uptake times dry biomass divided by liters gives M/s;
                    // fixed-order accumulation is deterministic; see DOCUMENTATION.md §10.3.
                    double& rate = concentration_rates[environment * substrate_count + substrate];
                    rate += uptake * buffers.biomass_gdw[agent] /
                        buffers.environment_volume_l[environment];
                    if (!is_finite_float64(rate) || rate < 0.0) {
                        throw std::runtime_error("environmental consumption rate overflowed");
                    }
                }
                growth_fluxes[agent] = fluxes[
                    agent * reaction_count + static_cast<std::size_t>(dense.biomass_global_index)];
                acetate_fluxes[agent] = fluxes[
                    agent * reaction_count + static_cast<std::size_t>(dense.acetate_global_index)];
            }

            const double remaining = dt_seconds - elapsed;
            double substep = remaining;
            for (std::size_t environment = 0; environment < buffers.environment_count; ++environment) {
                for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
                    const std::size_t index = environment * substrate_count + substrate;
                    const double rate = concentration_rates[index];
                    if (rate <= 0.0) {
                        continue;
                    }
                    // C/rate prevents depletion overshoot; below the reviewed glycan/MDCF
                    // threshold the safety factor tightens, never loosens, the limit;
                    // see DOCUMENTATION.md §10.3.
                    double depletion_limit = buffers.concentrations_m[index] / rate;
                    if (substrate >= 1 &&
                        buffers.concentrations_m[index] < low_substrate_threshold_m) {
                        depletion_limit *= depletion_safety;
                    }
                    substep = std::min(substep, depletion_limit);
                }
            }
            if (!is_finite_float64(substep) || substep <= 0.0 || elapsed + substep == elapsed) {
                throw std::runtime_error("adaptive metabolic substep made no numerical progress");
            }
            minimum_substep = std::min(minimum_substep, substep);

            for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
                for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
                    const std::size_t index = agent * substrate_count + substrate;
                    step_consumption[index] =
                        uptake_fluxes[index] * buffers.biomass_gdw[agent] * substep;
                    if (!is_finite_float64(step_consumption[index]) || step_consumption[index] < 0.0) {
                        throw std::runtime_error("substrate consumption overflowed");
                    }
                }
            }
            for (std::size_t environment = 0; environment < buffers.environment_count; ++environment) {
                for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
                    double requested = 0.0;
                    for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
                        if (static_cast<std::size_t>(buffers.environment_id[agent]) == environment) {
                            requested += step_consumption[agent * substrate_count + substrate];
                        }
                    }
                    const double available =
                        buffers.concentrations_m[environment * substrate_count + substrate] *
                        buffers.environment_volume_l[environment];
                    if (!is_finite_float64(requested) || !is_finite_float64(available)) {
                        throw std::runtime_error("substrate mass accounting overflowed");
                    }
                    if (requested > available && requested > 0.0) {
                        // Proportional scaling is a deterministic roundoff backstop that
                        // conserves mass across agents sharing an environment; see future
                        // DOCUMENTATION.md §10.3.
                        const double scale = available / requested;
                        for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
                            if (static_cast<std::size_t>(buffers.environment_id[agent]) == environment) {
                                step_consumption[agent * substrate_count + substrate] *= scale;
                            }
                        }
                    }
                }
            }
            for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
                step_biomass_delta[agent] =
                    growth_fluxes[agent] * buffers.biomass_gdw[agent] * substep;
                step_acetate[agent] =
                    acetate_fluxes[agent] * buffers.biomass_gdw[agent] * substep;
                if (!is_finite_float64(step_biomass_delta[agent]) ||
                    !is_finite_float64(step_acetate[agent])) {
                    throw std::runtime_error("metabolic output overflowed");
                }
            }
            for (std::size_t environment = 0; environment < buffers.environment_count; ++environment) {
                for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
                    double consumed = 0.0;
                    for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
                        if (static_cast<std::size_t>(buffers.environment_id[agent]) == environment) {
                            consumed += step_consumption[agent * substrate_count + substrate];
                        }
                    }
                    const std::size_t index = environment * substrate_count + substrate;
                    const double updated = buffers.concentrations_m[index] -
                        consumed / buffers.environment_volume_l[environment];
                    if (!is_finite_float64(updated)) {
                        throw std::runtime_error("concentration update overflowed");
                    }
                    next_concentrations[index] = std::max(0.0, updated);
                }
            }

            std::copy(next_concentrations.begin(), next_concentrations.end(), buffers.concentrations_m);
            for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
                for (std::size_t substrate = 0; substrate < substrate_count; ++substrate) {
                    const std::size_t index = agent * substrate_count + substrate;
                    buffers.substrate_consumed_mol[index] += step_consumption[index];
                }
                buffers.biomass_delta_gdw[agent] += step_biomass_delta[agent];
                buffers.acetate_secreted_mol[agent] += step_acetate[agent];
                growth_integral[agent] += growth_fluxes[agent] * substep;
                buffers.biomass_gdw[agent] += step_biomass_delta[agent];
                if (!is_finite_float64(buffers.biomass_gdw[agent]) ||
                    !is_finite_float64(buffers.biomass_delta_gdw[agent]) ||
                    !is_finite_float64(buffers.acetate_secreted_mol[agent]) ||
                    !is_finite_float64(growth_integral[agent])) {
                    throw std::runtime_error("cumulative metabolic output overflowed");
                }
            }
            elapsed += substep;
            ++summary.substeps;
        }
        for (std::size_t agent = 0; agent < buffers.agent_count; ++agent) {
            // The output is the time-weighted mean specific growth rate across internal
            // substeps, independent of render tick rate; see DOCUMENTATION.md §10.3.
            buffers.growth_rate_s[agent] = growth_integral[agent] / dt_seconds;
        }
        const CacheStats ending_stats = cache_stats();
        summary.lp_solves = ending_stats.lp_solves - starting_stats.lp_solves;
        summary.cache_hits = ending_stats.hits - starting_stats.hits;
        summary.cache_misses = ending_stats.misses - starting_stats.misses;
        summary.cache_evictions = ending_stats.evictions - starting_stats.evictions;
        summary.minimum_substep_seconds = minimum_substep;
        return summary;
    } catch (...) {
        rollback();
        throw;
    }
}

}  // namespace bacsandbox::metabolism