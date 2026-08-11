"""Versioned Phase C metabolic assets and solver-independent architecture."""

from __future__ import annotations

import copy
from dataclasses import dataclass
from enum import Enum
import importlib
import itertools
import json
import math
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping, Optional, Sequence, Tuple, Union

import numpy as np


SUBSTRATE_ORDER = (
    'generic_calorie',
    'complex_glycan',
    'mdcf_prebiotic',
)
MODEL_SCHEMA_VERSION = 1
MODEL_ID = 'reduced_metabolism_v1'
MODEL_REVISION = 1
PRIOR_SCHEMA_VERSION = 1
PRIOR_ID = 'metabolic_priors_v1'
PRIOR_REVISION = 1
PROVENANCE = 'Own (BacSandbox model)'
CACHE_TOLERANCE_MOLAR = 1.0e-5
DEFAULT_CACHE_CAPACITY = 256
MAX_CACHE_CAPACITY = 1_000_000
SIMPLEX_TOLERANCE = 1.0e-10
INT64_MAX = int(np.iinfo(np.int64).max)
MODEL_DIRECTORY = Path(__file__).resolve().parent / 'models'
DEFAULT_MODEL_PATH = MODEL_DIRECTORY / 'reduced_metabolism_v1.json'
DEFAULT_PRIOR_PATH = MODEL_DIRECTORY / 'metabolic_priors_v1.json'

# Own-model base capacities keep reduced whole-cell exchange in the single-digit
# micromole/gDW/s range and make generic monomer uptake less preprocessing-limited
# than glycan-derived pools; see DOCUMENTATION.md section 10.4.
SPECIES_BASE_UPTAKE_VMAX_MOL_GDW_S = (2.8e-6, 1.4e-6, 1.8e-6)
# Positive 0.10-0.20 mM half-saturation concentrations keep the Phase 1 MDCF bolus
# concentration-responsive while separating Km from capacity; see DOCUMENTATION.md section 10.4.
SPECIES_UPTAKE_KM_MOLAR = (1.0e-4, 2.0e-4, 1.5e-4)
SPECIES_AFFINITY_FIELDS = ('aff_calorie', 'aff_glycan', 'aff_mdcf')
# Missing legacy fields remain calorie-only, preserving the pre-diet fallback;
# see DOCUMENTATION.md section 10.4.
SPECIES_DEFAULT_AFFINITY_MULTIPLIERS = (1.0, 0.0, 0.0)
# Shipped CSV profiles use Own-model relative tiers: glycan 2, MAG-like MDCF 3,
# cultured BgF5_2 MDCF 4, and explicit MDCF non-utilization 0; see DOCUMENTATION.md section 10.4.

_MODEL_KEYS = frozenset((
    'schema_version',
    'model_id',
    'model_revision',
    'provenance',
    'scope',
    'design_rationale',
    'units',
    'substrate_order',
    'metabolite_order',
    'reaction_order',
    'stoichiometry',
    'lower_bounds',
    'upper_bounds',
    'objective_coefficients',
    'exchange_reaction_indices',
    'biomass_reaction_index',
    'acetate_reaction_index',
    'cache_tolerance_molar',
    'requires_explicit_environment_volume_liters',
    'finite_bound_semantics',
    'analytical_reference',
))
_SCOPE_KEYS = frozenset(('kind', 'genome_scale_model', 'statement'))
_DESIGN_RATIONALE_KEYS = frozenset((
    'routing',
    'acetate_coupling',
    'finite_bounds',
    'objective',
    'cache_tolerance',
))
_UNIT_KEYS = frozenset((
    'substrate_concentration',
    'environment_volume',
    'reaction_flux',
    'stoichiometric_coefficient',
))
_STOICHIOMETRY_KEYS = frozenset(('metabolite_index', 'reaction_index', 'value'))
_ANALYTICAL_REFERENCE_KEYS = frozenset((
    'condition',
    'biomass_flux',
    'acetate_flux',
    'statement',
))
_PRIOR_KEYS = frozenset((
    'schema_version',
    'prior_id',
    'prior_revision',
    'provenance',
    'scope',
    'parameter_rationale',
    'substrate_order',
    'rows',
))
_PARAMETER_RATIONALE_KEYS = frozenset(('alpha_beta', 'max_inhibition'))
_PRIOR_ROW_KEYS = frozenset((
    'phenotype_id',
    'substrate',
    'alpha',
    'beta',
    'max_inhibition',
    'evidence_scope',
))
_EXPECTED_UNITS = {
    'substrate_concentration': 'M glucose-equivalent monomer',
    'environment_volume': 'L',
    'reaction_flux': (
        'Exchange and route fluxes are mol glucose-equivalent monomer / gDW / s; '
        'biomass pseudo-flux is a specific growth rate in 1/s; acetate export is mol '
        'glucose-equivalent / gDW / s. These are reduced bookkeeping units, not a '
        'claim of dimensional GEM rigor.'
    ),
    'stoichiometric_coefficient': 'dimensionless glucose-equivalent bookkeeping coefficient',
}


class MetabolicAssetError(ValueError):
    """Raised when a versioned metabolic asset violates its schema."""


class MetabolicSolverUnavailable(NotImplementedError):
    """Raised when a Phase C solver operation is requested before Phase 2."""


class MetabolicSolverError(RuntimeError):
    """Raised when the embedded LP solver cannot return a feasible optimum."""


class MetabolicMode(str, Enum):
    LEGACY = 'legacy'
    SHADOW = 'shadow'
    FBA = 'fba'


class MetabolicBackend(str, Enum):
    AUTO = 'auto'
    PYTHON = 'python'
    CPP = 'cpp'


@dataclass(frozen=True)
class MetabolicModelData:
    schema_version: int
    model_id: str
    model_revision: int
    provenance: str
    scope_statement: str
    design_rationale: Mapping[str, str]
    units: Mapping[str, str]
    substrate_order: Tuple[str, str, str]
    metabolite_order: Tuple[str, ...]
    reaction_order: Tuple[str, ...]
    cache_tolerance_molar: float
    requires_explicit_environment_volume_liters: bool
    finite_bound_semantics: str
    analytical_full_uptake_biomass_flux: float
    analytical_full_uptake_acetate_flux: float
    model_reaction_offsets: np.ndarray
    model_metabolite_offsets: np.ndarray
    stoich_row_offsets: np.ndarray
    stoich_column_indices: np.ndarray
    stoich_values: np.ndarray
    lower_bounds: np.ndarray
    upper_bounds: np.ndarray
    objective_coefficients: np.ndarray
    exchange_reaction_indices: np.ndarray
    biomass_reaction_indices: np.ndarray
    acetate_reaction_indices: np.ndarray

    @property
    def model_count(self) -> int:
        return int(self.model_reaction_offsets.size - 1)


@dataclass(frozen=True)
class MetabolicPriorRow:
    phenotype_id: str
    substrate: str
    alpha: float
    beta: float
    max_inhibition: float
    evidence_scope: str


@dataclass(frozen=True)
class MetabolicPriorData:
    schema_version: int
    prior_id: str
    prior_revision: int
    provenance: str
    scope: str
    parameter_rationale: Mapping[str, str]
    substrate_order: Tuple[str, str, str]
    rows: Tuple[MetabolicPriorRow, ...]


def _reject_json_constant(value: str) -> None:
    raise MetabolicAssetError(f'JSON constants must be finite numbers, got {value}')


def _reject_duplicate_keys(pairs: Sequence[Tuple[str, Any]]) -> Mapping[str, Any]:
    result = {}
    for key, value in pairs:
        if key in result:
            raise MetabolicAssetError(f'duplicate JSON key: {key}')
        result[key] = value
    return result


def _read_json_object(path: Union[str, Path]) -> Mapping[str, Any]:
    asset_path = Path(path)
    try:
        source = asset_path.read_text(encoding='utf-8')
        value = json.loads(
            source,
            object_pairs_hook=_reject_duplicate_keys,
            parse_constant=_reject_json_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise MetabolicAssetError(f'could not load metabolic asset {asset_path}: {exc}') from exc
    if not isinstance(value, Mapping):
        raise MetabolicAssetError('metabolic asset root must be a JSON object')
    return value


def _require_exact_keys(value: Any, expected: frozenset, context: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise MetabolicAssetError(f'{context} must be an object')
    actual = frozenset(value)
    if actual != expected:
        missing = sorted(expected - actual)
        unexpected = sorted(actual - expected)
        raise MetabolicAssetError(
            f'{context} keys differ: missing={missing}, unexpected={unexpected}'
        )
    return value


def _require_int(value: Any, context: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise MetabolicAssetError(f'{context} must be an integer')
    return value


def _require_number(value: Any, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise MetabolicAssetError(f'{context} must be a number')
    result = float(value)
    if not math.isfinite(result):
        raise MetabolicAssetError(f'{context} must be finite')
    return result


def _require_string(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise MetabolicAssetError(f'{context} must be a nonempty string')
    return value


def _require_string_order(value: Any, context: str) -> Tuple[str, ...]:
    if not isinstance(value, list) or not value:
        raise MetabolicAssetError(f'{context} must be a nonempty array')
    result = tuple(_require_string(item, f'{context}[{index}]') for index, item in enumerate(value))
    if len(set(result)) != len(result):
        raise MetabolicAssetError(f'{context} entries must be unique')
    return result


def _require_index(value: Any, size: int, context: str) -> int:
    index = _require_int(value, context)
    if index < 0 or index >= size:
        raise MetabolicAssetError(f'{context} must be in [0, {size})')
    return index


def _require_numeric_vector(value: Any, size: int, context: str) -> Tuple[float, ...]:
    if not isinstance(value, list) or len(value) != size:
        raise MetabolicAssetError(f'{context} must have length {size}')
    return tuple(_require_number(item, f'{context}[{index}]') for index, item in enumerate(value))


def _readonly_array(value: Any, dtype: Any) -> np.ndarray:
    array = np.ascontiguousarray(value, dtype=dtype)
    array.setflags(write=False)
    return array


def load_metabolic_model(
    path: Union[str, Path] = DEFAULT_MODEL_PATH,
) -> MetabolicModelData:
    root = _require_exact_keys(_read_json_object(path), _MODEL_KEYS, 'model')
    schema_version = _require_int(root['schema_version'], 'schema_version')
    model_id = _require_string(root['model_id'], 'model_id')
    model_revision = _require_int(root['model_revision'], 'model_revision')
    if (schema_version, model_id, model_revision) != (
        MODEL_SCHEMA_VERSION,
        MODEL_ID,
        MODEL_REVISION,
    ):
        raise MetabolicAssetError(
            'unsupported model schema/id/revision: '
            f'{schema_version}/{model_id}/{model_revision}'
        )
    provenance = _require_string(root['provenance'], 'provenance')
    if provenance != PROVENANCE:
        raise MetabolicAssetError(f'provenance must be {PROVENANCE!r}')

    scope = _require_exact_keys(root['scope'], _SCOPE_KEYS, 'scope')
    if scope['kind'] != 'reduced_three_substrate':
        raise MetabolicAssetError('scope.kind must be reduced_three_substrate')
    if scope['genome_scale_model'] is not False:
        raise MetabolicAssetError('scope.genome_scale_model must be false')
    scope_statement = _require_string(scope['statement'], 'scope.statement')
    if 'not a genome-scale model' not in scope_statement.lower():
        raise MetabolicAssetError('scope.statement must explicitly disclaim GEM scope')
    if (
        'glucose-equivalent monomer concentrations' not in scope_statement.lower()
        or 'never literal polymer molarity' not in scope_statement.lower()
    ):
        raise MetabolicAssetError(
            'scope.statement must define glucose-equivalent, non-polymer-molarity substrates'
        )

    design_rationale = _require_exact_keys(
        root['design_rationale'],
        _DESIGN_RATIONALE_KEYS,
        'design_rationale',
    )
    for rationale_name, rationale in design_rationale.items():
        _require_string(rationale, f'design_rationale.{rationale_name}')
    if 'not empirical' not in design_rationale['acetate_coupling'].lower():
        raise MetabolicAssetError('acetate coupling rationale must disclaim empirical status')
    if 'not a biological threshold' not in design_rationale['cache_tolerance'].lower():
        raise MetabolicAssetError('cache tolerance rationale must disclaim biological-threshold status')

    units = _require_exact_keys(root['units'], _UNIT_KEYS, 'units')
    for unit_name, expected_unit in _EXPECTED_UNITS.items():
        if units[unit_name] != expected_unit:
            raise MetabolicAssetError(f'units.{unit_name} must be {expected_unit!r}')

    substrate_order = _require_string_order(root['substrate_order'], 'substrate_order')
    if substrate_order != SUBSTRATE_ORDER:
        raise MetabolicAssetError(f'substrate_order must be exactly {SUBSTRATE_ORDER!r}')
    metabolite_order = _require_string_order(root['metabolite_order'], 'metabolite_order')
    reaction_order = _require_string_order(root['reaction_order'], 'reaction_order')
    metabolite_count = len(metabolite_order)
    reaction_count = len(reaction_order)

    lower_bounds = _require_numeric_vector(root['lower_bounds'], reaction_count, 'lower_bounds')
    upper_bounds = _require_numeric_vector(root['upper_bounds'], reaction_count, 'upper_bounds')
    objective_coefficients = _require_numeric_vector(
        root['objective_coefficients'],
        reaction_count,
        'objective_coefficients',
    )
    for index, (lower_bound, upper_bound) in enumerate(zip(lower_bounds, upper_bounds)):
        if lower_bound > upper_bound:
            raise MetabolicAssetError(f'reaction bound ordering is invalid at index {index}')

    exchange_values = root['exchange_reaction_indices']
    if not isinstance(exchange_values, list) or len(exchange_values) != len(SUBSTRATE_ORDER):
        raise MetabolicAssetError('exchange_reaction_indices must have length 3')
    exchange_indices = tuple(
        _require_index(value, reaction_count, f'exchange_reaction_indices[{index}]')
        for index, value in enumerate(exchange_values)
    )
    if len(set(exchange_indices)) != len(exchange_indices):
        raise MetabolicAssetError('exchange_reaction_indices entries must be unique')
    biomass_index = _require_index(
        root['biomass_reaction_index'],
        reaction_count,
        'biomass_reaction_index',
    )
    acetate_index = _require_index(
        root['acetate_reaction_index'],
        reaction_count,
        'acetate_reaction_index',
    )
    if biomass_index == acetate_index or biomass_index in exchange_indices or acetate_index in exchange_indices:
        raise MetabolicAssetError('biomass, acetate, and exchange reaction indices must be distinct')
    if objective_coefficients[biomass_index] <= 0.0:
        raise MetabolicAssetError('biomass objective coefficient must be positive')
    if any(
        coefficient != 0.0
        for index, coefficient in enumerate(objective_coefficients)
        if index != biomass_index
    ):
        raise MetabolicAssetError('only the biomass reaction may have a nonzero objective coefficient')

    entries = root['stoichiometry']
    if not isinstance(entries, list) or not entries:
        raise MetabolicAssetError('stoichiometry must be a nonempty array')
    coordinate_values = {}
    covered_reactions = set()
    for entry_index, raw_entry in enumerate(entries):
        entry = _require_exact_keys(raw_entry, _STOICHIOMETRY_KEYS, f'stoichiometry[{entry_index}]')
        metabolite_index = _require_index(
            entry['metabolite_index'],
            metabolite_count,
            f'stoichiometry[{entry_index}].metabolite_index',
        )
        reaction_index = _require_index(
            entry['reaction_index'],
            reaction_count,
            f'stoichiometry[{entry_index}].reaction_index',
        )
        coefficient = _require_number(entry['value'], f'stoichiometry[{entry_index}].value')
        if coefficient == 0.0:
            raise MetabolicAssetError('stoichiometry entries must be nonzero')
        coordinate = (metabolite_index, reaction_index)
        if coordinate in coordinate_values:
            raise MetabolicAssetError(f'duplicate stoichiometry coordinate: {coordinate}')
        coordinate_values[coordinate] = coefficient
        covered_reactions.add(reaction_index)
    if covered_reactions != set(range(reaction_count)):
        raise MetabolicAssetError('every reaction must have at least one stoichiometry entry')

    stoich_row_offsets = [0]
    stoich_column_indices = []
    stoich_values = []
    for metabolite_index in range(metabolite_count):
        row_entries = sorted(
            (reaction_index, coefficient)
            for (row_index, reaction_index), coefficient in coordinate_values.items()
            if row_index == metabolite_index
        )
        if not row_entries:
            raise MetabolicAssetError(f'metabolite index {metabolite_index} has no stoichiometry')
        for reaction_index, coefficient in row_entries:
            stoich_column_indices.append(reaction_index)
            stoich_values.append(coefficient)
        stoich_row_offsets.append(len(stoich_values))

    cache_tolerance = _require_number(root['cache_tolerance_molar'], 'cache_tolerance_molar')
    if cache_tolerance != CACHE_TOLERANCE_MOLAR:
        raise MetabolicAssetError(
            f'cache_tolerance_molar must be {CACHE_TOLERANCE_MOLAR!r} for revision 1'
        )
    if root['requires_explicit_environment_volume_liters'] is not True:
        raise MetabolicAssetError('explicit environment volume in liters must be required')
    finite_bound_semantics = _require_string(
        root['finite_bound_semantics'],
        'finite_bound_semantics',
    )
    if 'finite' not in finite_bound_semantics.lower():
        raise MetabolicAssetError('finite_bound_semantics must explain finite bounds')

    analytical = _require_exact_keys(
        root['analytical_reference'],
        _ANALYTICAL_REFERENCE_KEYS,
        'analytical_reference',
    )
    _require_string(analytical['condition'], 'analytical_reference.condition')
    analytical_biomass = _require_number(
        analytical['biomass_flux'],
        'analytical_reference.biomass_flux',
    )
    analytical_acetate = _require_number(
        analytical['acetate_flux'],
        'analytical_reference.acetate_flux',
    )
    analytical_statement = _require_string(
        analytical['statement'],
        'analytical_reference.statement',
    )
    if 'not an empirical prediction' not in analytical_statement.lower():
        raise MetabolicAssetError('analytical_reference must disclaim empirical prediction')

    return MetabolicModelData(
        schema_version=schema_version,
        model_id=model_id,
        model_revision=model_revision,
        provenance=provenance,
        scope_statement=scope_statement,
        design_rationale=MappingProxyType(dict(design_rationale)),
        units=MappingProxyType(dict(units)),
        substrate_order=substrate_order,
        metabolite_order=metabolite_order,
        reaction_order=reaction_order,
        cache_tolerance_molar=cache_tolerance,
        requires_explicit_environment_volume_liters=True,
        finite_bound_semantics=finite_bound_semantics,
        analytical_full_uptake_biomass_flux=analytical_biomass,
        analytical_full_uptake_acetate_flux=analytical_acetate,
        model_reaction_offsets=_readonly_array((0, reaction_count), np.int64),
        model_metabolite_offsets=_readonly_array((0, metabolite_count), np.int64),
        stoich_row_offsets=_readonly_array(stoich_row_offsets, np.int64),
        stoich_column_indices=_readonly_array(stoich_column_indices, np.int32),
        stoich_values=_readonly_array(stoich_values, np.float64),
        lower_bounds=_readonly_array(lower_bounds, np.float64),
        upper_bounds=_readonly_array(upper_bounds, np.float64),
        objective_coefficients=_readonly_array(objective_coefficients, np.float64),
        exchange_reaction_indices=_readonly_array((exchange_indices,), np.int32),
        biomass_reaction_indices=_readonly_array((biomass_index,), np.int32),
        acetate_reaction_indices=_readonly_array((acetate_index,), np.int32),
    )


def load_metabolic_priors(
    path: Union[str, Path] = DEFAULT_PRIOR_PATH,
) -> MetabolicPriorData:
    root = _require_exact_keys(_read_json_object(path), _PRIOR_KEYS, 'priors')
    schema_version = _require_int(root['schema_version'], 'schema_version')
    prior_id = _require_string(root['prior_id'], 'prior_id')
    prior_revision = _require_int(root['prior_revision'], 'prior_revision')
    if (schema_version, prior_id, prior_revision) != (
        PRIOR_SCHEMA_VERSION,
        PRIOR_ID,
        PRIOR_REVISION,
    ):
        raise MetabolicAssetError(
            'unsupported prior schema/id/revision: '
            f'{schema_version}/{prior_id}/{prior_revision}'
        )
    provenance = _require_string(root['provenance'], 'provenance')
    if provenance != PROVENANCE:
        raise MetabolicAssetError(f'provenance must be {PROVENANCE!r}')
    scope = _require_string(root['scope'], 'scope')
    lowered_scope = scope.lower()
    if 'neutral' not in lowered_scope or 'no empirical interaction' not in lowered_scope:
        raise MetabolicAssetError('prior scope must state neutral/no-empirical-interaction status')
    parameter_rationale = _require_exact_keys(
        root['parameter_rationale'],
        _PARAMETER_RATIONALE_KEYS,
        'parameter_rationale',
    )
    for rationale_name, rationale in parameter_rationale.items():
        _require_string(rationale, f'parameter_rationale.{rationale_name}')
    if 'no empirical interaction' not in parameter_rationale['alpha_beta'].lower():
        raise MetabolicAssetError('alpha/beta rationale must disclaim empirical interaction evidence')
    if 'exactly neutral' not in parameter_rationale['max_inhibition'].lower():
        raise MetabolicAssetError('max-inhibition rationale must state neutral production behavior')
    substrate_order = _require_string_order(root['substrate_order'], 'substrate_order')
    if substrate_order != SUBSTRATE_ORDER:
        raise MetabolicAssetError(f'substrate_order must be exactly {SUBSTRATE_ORDER!r}')

    raw_rows = root['rows']
    if not isinstance(raw_rows, list) or not raw_rows:
        raise MetabolicAssetError('rows must be a nonempty array')
    rows = []
    row_keys = set()
    for row_index, raw_row in enumerate(raw_rows):
        row = _require_exact_keys(raw_row, _PRIOR_ROW_KEYS, f'rows[{row_index}]')
        phenotype_id = _require_string(row['phenotype_id'], f'rows[{row_index}].phenotype_id')
        substrate = _require_string(row['substrate'], f'rows[{row_index}].substrate')
        if substrate not in SUBSTRATE_ORDER:
            raise MetabolicAssetError(f'rows[{row_index}].substrate is unknown')
        key = (phenotype_id, substrate)
        if key in row_keys:
            raise MetabolicAssetError(f'duplicate prior row: {key}')
        row_keys.add(key)
        alpha = _require_number(row['alpha'], f'rows[{row_index}].alpha')
        beta = _require_number(row['beta'], f'rows[{row_index}].beta')
        max_inhibition = _require_number(
            row['max_inhibition'],
            f'rows[{row_index}].max_inhibition',
        )
        if alpha <= 0.0 or beta <= 0.0:
            raise MetabolicAssetError('prior alpha and beta must be positive')
        if max_inhibition != 0.0:
            raise MetabolicAssetError('production prior max_inhibition must be zero')
        evidence_scope = _require_string(
            row['evidence_scope'],
            f'rows[{row_index}].evidence_scope',
        )
        if 'no empirical interaction' not in evidence_scope.lower():
            raise MetabolicAssetError('prior rows must disclaim empirical interaction evidence')
        rows.append(MetabolicPriorRow(
            phenotype_id=phenotype_id,
            substrate=substrate,
            alpha=alpha,
            beta=beta,
            max_inhibition=max_inhibition,
            evidence_scope=evidence_scope,
        ))
    if tuple(row.substrate for row in rows) != SUBSTRATE_ORDER:
        raise MetabolicAssetError('revision 1 prior rows must follow exact substrate order')

    return MetabolicPriorData(
        schema_version=schema_version,
        prior_id=prior_id,
        prior_revision=prior_revision,
        provenance=provenance,
        scope=scope,
        parameter_rationale=MappingProxyType(dict(parameter_rationale)),
        substrate_order=substrate_order,
        rows=tuple(rows),
    )


def _finite_float(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise TypeError(f'{name} must be a real number')
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f'{name} must be finite')
    return result


def _three_concentrations(concentrations: Sequence[float]) -> Tuple[float, float, float]:
    if isinstance(concentrations, (str, bytes)):
        raise TypeError('concentrations must be a three-value sequence')
    try:
        values = tuple(concentrations)
    except TypeError as exc:
        raise TypeError('concentrations must be a three-value sequence') from exc
    if len(values) != len(SUBSTRATE_ORDER):
        raise ValueError('concentrations must contain exactly three substrate values')
    result = tuple(_finite_float(value, f'concentrations[{index}]') for index, value in enumerate(values))
    if any(value < 0.0 for value in result):
        raise ValueError('concentrations must be nonnegative')
    return result  # type: ignore[return-value]


def quantize_cache_key(
    concentrations: Sequence[float],
    tolerance_molar: float = CACHE_TOLERANCE_MOLAR,
) -> Tuple[int, int, int]:
    values = _three_concentrations(concentrations)
    tolerance = _finite_float(tolerance_molar, 'tolerance_molar')
    if tolerance <= 0.0:
        raise ValueError('tolerance_molar must be positive')
    # Flooring each concentration bin creates deterministic one-sided cache cells at exact
    # tolerance boundaries; see DOCUMENTATION.md §10.3.
    scaled_values = tuple(value / tolerance for value in values)
    if any(not math.isfinite(value) or value > INT64_MAX for value in scaled_values):
        raise ValueError('quantized concentrations must fit signed int64 cache bins')
    return tuple(math.floor(value) for value in scaled_values)  # type: ignore[return-value]


def enumerate_neighbor_cache_keys(cache_key: Sequence[int]) -> Tuple[Tuple[int, int, int], ...]:
    if isinstance(cache_key, (str, bytes)):
        raise TypeError('cache_key must be a three-integer sequence')
    try:
        values = tuple(cache_key)
    except TypeError as exc:
        raise TypeError('cache_key must be a three-integer sequence') from exc
    if len(values) != len(SUBSTRATE_ORDER):
        raise ValueError('cache_key must contain exactly three values')
    for index, value in enumerate(values):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise TypeError(f'cache_key[{index}] must be an integer')
        if int(value) < 0:
            raise ValueError('cache_key entries must be nonnegative')
        if int(value) >= INT64_MAX:
            raise ValueError('cache_key entries must leave room for the positive neighbor')
    normalized = tuple(int(value) for value in values)
    # The Cartesian {-1, 0, 1} offsets inspect all 3^3 adjacent bins in stable substrate
    # order, including the center bin; see DOCUMENTATION.md §10.3.
    return tuple(
        tuple(normalized[index] + offsets[index] for index in range(3))
        for offsets in itertools.product((-1, 0, 1), repeat=3)
    )  # type: ignore[return-value]


def posterior_mean_inhibition_multiplier(
    alpha: float,
    beta: float,
    max_inhibition: float,
) -> float:
    alpha_value = _finite_float(alpha, 'alpha')
    beta_value = _finite_float(beta, 'beta')
    inhibition_value = _finite_float(max_inhibition, 'max_inhibition')
    if alpha_value <= 0.0 or beta_value <= 0.0:
        raise ValueError('alpha and beta must be positive')
    if inhibition_value < 0.0 or inhibition_value > 1.0:
        raise ValueError('max_inhibition must be in [0, 1]')
    posterior_mass = alpha_value + beta_value
    if not math.isfinite(posterior_mass):
        raise ValueError('alpha + beta must remain finite')
    # The beta posterior mean alpha/(alpha+beta) scales the configured maximum inhibition,
    # leaving a multiplicative active fraction; see DOCUMENTATION.md §10.3.
    return 1.0 - alpha_value / posterior_mass * inhibition_value


def concentration_dependent_uptake_bound(
    vmax: float,
    concentration: float,
    km: float,
    multiplier: float = 1.0,
) -> float:
    vmax_value = _finite_float(vmax, 'vmax')
    concentration_value = _finite_float(concentration, 'concentration')
    km_value = _finite_float(km, 'km')
    multiplier_value = _finite_float(multiplier, 'multiplier')
    if vmax_value < 0.0 or concentration_value < 0.0:
        raise ValueError('vmax and concentration must be nonnegative')
    if km_value <= 0.0:
        raise ValueError('km must be positive')
    if multiplier_value < 0.0 or multiplier_value > 1.0:
        raise ValueError('multiplier must be in [0, 1]')
    denominator = km_value + concentration_value
    if not math.isfinite(denominator):
        raise ValueError('km + concentration must remain finite')
    # Saturating uptake allocates vmax*C/(km+C), then applies the independently bounded
    # inhibition multiplier; see DOCUMENTATION.md §10.3.
    return concentration_value / denominator * vmax_value * multiplier_value


def _coerce_enum(value: Any, enum_type: Any, name: str) -> Any:
    try:
        return enum_type(value)
    except (TypeError, ValueError) as exc:
        choices = ', '.join(member.value for member in enum_type)
        raise ValueError(f'{name} must be one of {choices}') from exc


class _LPStatus(str, Enum):
    OPTIMAL = 'optimal'
    INFEASIBLE = 'infeasible'
    UNBOUNDED = 'unbounded'
    NUMERICAL_ERROR = 'numerical_error'


@dataclass(frozen=True)
class _DenseLPModel:
    reaction_start: int
    reaction_stop: int
    stoichiometry: np.ndarray
    lower_bounds: np.ndarray
    upper_bounds: np.ndarray
    objective: np.ndarray
    exchange_indices: Tuple[int, int, int]
    biomass_index: int
    acetate_index: int


@dataclass
class _CacheEntry:
    concentrations: np.ndarray
    effective_exchange_bounds: np.ndarray
    flux_values: np.ndarray
    insertion_sequence: int
    last_used_sequence: int


def _require_exact_array(
    value: Any,
    dtype: Any,
    ndim: int,
    name: str,
    shape: Optional[Tuple[Optional[int], ...]] = None,
    writable: bool = False,
) -> np.ndarray:
    if not isinstance(value, np.ndarray):
        raise TypeError(f'{name} must be a NumPy array')
    if value.dtype != np.dtype(dtype):
        raise TypeError(f'{name} has an invalid dtype')
    if value.ndim != ndim:
        raise ValueError(f'{name} has an invalid shape')
    if shape is not None and any(
        expected is not None and value.shape[index] != expected
        for index, expected in enumerate(shape)
    ):
        raise ValueError(f'{name} has an invalid shape')
    if not value.flags.c_contiguous:
        raise ValueError(f'{name} must be C-contiguous')
    if writable and not value.flags.writeable:
        raise ValueError(f'{name} must be writable')
    return value


def _simplex_pivot(
    tableau: np.ndarray,
    basis: np.ndarray,
    nonbasis: np.ndarray,
    pivot_row: int,
    pivot_column: int,
) -> None:
    row_count = basis.size
    column_count = nonbasis.size - 1
    inverse = 1.0 / tableau[pivot_row, pivot_column]
    # Gauss-Jordan elimination preserves the current basic feasible dictionary; the
    # fixed row/column order makes every degenerate pivot reproducible; see future
    # DOCUMENTATION.md §10.3 (Embedded Metabolic LP Solver).
    for row in range(row_count + 2):
        if row == pivot_row:
            continue
        for column in range(column_count + 2):
            if column != pivot_column:
                tableau[row, column] -= (
                    tableau[pivot_row, column]
                    * tableau[row, pivot_column]
                    * inverse
                )
    for column in range(column_count + 2):
        if column != pivot_column:
            tableau[pivot_row, column] *= inverse
    for row in range(row_count + 2):
        if row != pivot_row:
            tableau[row, pivot_column] *= -inverse
    tableau[pivot_row, pivot_column] = inverse
    basis[pivot_row], nonbasis[pivot_column] = (
        nonbasis[pivot_column],
        basis[pivot_row],
    )


def _simplex_phase(
    tableau: np.ndarray,
    basis: np.ndarray,
    nonbasis: np.ndarray,
    phase: int,
    tolerance: float,
) -> bool:
    row_count = basis.size
    column_count = nonbasis.size - 1
    objective_row = row_count + 1 if phase == 1 else row_count
    while True:
        entering = -1
        for column in range(column_count + 1):
            if phase == 2 and nonbasis[column] == -1:
                continue
            if entering < 0 or (
                tableau[objective_row, column]
                < tableau[objective_row, entering] - tolerance
            ) or (
                abs(
                    tableau[objective_row, column]
                    - tableau[objective_row, entering]
                ) <= tolerance
                and nonbasis[column] < nonbasis[entering]
            ):
                entering = column
        if entering < 0 or tableau[objective_row, entering] >= -tolerance:
            return True

        leaving = -1
        for row in range(row_count):
            if tableau[row, entering] <= tolerance:
                continue
            if leaving < 0:
                leaving = row
                continue
            candidate_ratio = tableau[row, column_count + 1] / tableau[row, entering]
            current_ratio = (
                tableau[leaving, column_count + 1]
                / tableau[leaving, entering]
            )
            if candidate_ratio < current_ratio - tolerance or (
                abs(candidate_ratio - current_ratio) <= tolerance
                and basis[row] < basis[leaving]
            ):
                leaving = row
        if leaving < 0:
            return False
        _simplex_pivot(tableau, basis, nonbasis, leaving, entering)


def _solve_nonnegative_lp(
    constraint_matrix: np.ndarray,
    constraint_bounds: np.ndarray,
    objective: np.ndarray,
    tolerance: float = SIMPLEX_TOLERANCE,
) -> Tuple[_LPStatus, np.ndarray]:
    row_count, column_count = constraint_matrix.shape
    tableau = np.zeros((row_count + 2, column_count + 2), dtype=np.float64)
    basis = np.empty(row_count, dtype=np.int64)
    nonbasis = np.empty(column_count + 1, dtype=np.int64)
    tableau[:row_count, :column_count] = constraint_matrix
    for row in range(row_count):
        basis[row] = column_count + row
        tableau[row, column_count] = -1.0
        tableau[row, column_count + 1] = constraint_bounds[row]
    for column in range(column_count):
        nonbasis[column] = column
        tableau[row_count, column] = -objective[column]
    nonbasis[column_count] = -1
    tableau[row_count + 1, column_count] = 1.0

    if row_count:
        pivot_row = min(
            range(row_count),
            key=lambda row: (tableau[row, column_count + 1], basis[row]),
        )
        if tableau[pivot_row, column_count + 1] < -tolerance:
            _simplex_pivot(tableau, basis, nonbasis, pivot_row, column_count)
            if not _simplex_phase(tableau, basis, nonbasis, 1, tolerance):
                return _LPStatus.INFEASIBLE, np.zeros(column_count, dtype=np.float64)
            if abs(tableau[row_count + 1, column_count + 1]) > tolerance:
                return _LPStatus.INFEASIBLE, np.zeros(column_count, dtype=np.float64)
            artificial_rows = np.flatnonzero(basis == -1)
            for row_value in artificial_rows:
                row = int(row_value)
                pivot_column = min(
                    range(column_count + 1),
                    key=lambda column: (
                        abs(tableau[row, column]) <= tolerance,
                        nonbasis[column],
                    ),
                )
                if abs(tableau[row, pivot_column]) > tolerance:
                    _simplex_pivot(tableau, basis, nonbasis, row, pivot_column)

    if not _simplex_phase(tableau, basis, nonbasis, 2, tolerance):
        return _LPStatus.UNBOUNDED, np.zeros(column_count, dtype=np.float64)
    solution = np.zeros(column_count, dtype=np.float64)
    for row in range(row_count):
        if 0 <= basis[row] < column_count:
            solution[basis[row]] = tableau[row, column_count + 1]
    return _LPStatus.OPTIMAL, solution


def _solve_bounded_equality_lp(
    stoichiometry: np.ndarray,
    lower_bounds: np.ndarray,
    upper_bounds: np.ndarray,
    objective: np.ndarray,
    tolerance: float = SIMPLEX_TOLERANCE,
) -> Tuple[_LPStatus, np.ndarray]:
    reaction_count = lower_bounds.size
    widths = upper_bounds - lower_bounds
    # Shifting x=l+y converts finite reaction bounds to y>=0 and y<=u-l; paired
    # inequalities retain exact steady-state equalities Sx=0; see future
    # DOCUMENTATION.md §10.3 (Embedded Metabolic LP Solver).
    equality_bounds = -(stoichiometry @ lower_bounds)
    constraint_matrix = np.vstack((
        stoichiometry,
        -stoichiometry,
        np.eye(reaction_count, dtype=np.float64),
    ))
    constraint_bounds = np.concatenate((equality_bounds, -equality_bounds, widths))
    status, shifted = _solve_nonnegative_lp(
        constraint_matrix,
        constraint_bounds,
        objective,
        tolerance,
    )
    if status is not _LPStatus.OPTIMAL:
        return status, np.zeros(reaction_count, dtype=np.float64)
    solution = lower_bounds + shifted
    solution[np.abs(solution) <= tolerance] = 0.0
    solution = np.minimum(upper_bounds, np.maximum(lower_bounds, solution))
    residual = stoichiometry @ solution
    if (
        np.max(np.abs(residual), initial=0.0) > 10.0 * tolerance
        or np.any(solution < lower_bounds - 10.0 * tolerance)
        or np.any(solution > upper_bounds + 10.0 * tolerance)
    ):
        return _LPStatus.NUMERICAL_ERROR, np.zeros(reaction_count, dtype=np.float64)
    return status, solution


class PythonMetabolicKernel:
    """Deterministic Python reference implementation of the final metabolic ABI."""

    def __init__(
        self,
        model_reaction_offsets: np.ndarray,
        model_metabolite_offsets: np.ndarray,
        stoich_row_offsets: np.ndarray,
        stoich_column_indices: np.ndarray,
        stoich_values: np.ndarray,
        lower_bounds: np.ndarray,
        upper_bounds: np.ndarray,
        objective_coefficients: np.ndarray,
        exchange_reaction_indices: np.ndarray,
        biomass_reaction_indices: np.ndarray,
        acetate_reaction_indices: np.ndarray,
        cache_tolerance_molar: float,
        cache_capacity: int,
    ) -> None:
        reaction_offsets = _require_exact_array(
            model_reaction_offsets, np.int64, 1, 'model_reaction_offsets'
        )
        metabolite_offsets = _require_exact_array(
            model_metabolite_offsets, np.int64, 1, 'model_metabolite_offsets'
        )
        row_offsets = _require_exact_array(
            stoich_row_offsets, np.int64, 1, 'stoich_row_offsets'
        )
        column_indices = _require_exact_array(
            stoich_column_indices, np.int32, 1, 'stoich_column_indices'
        )
        values = _require_exact_array(stoich_values, np.float64, 1, 'stoich_values')
        lower = _require_exact_array(lower_bounds, np.float64, 1, 'lower_bounds')
        upper = _require_exact_array(upper_bounds, np.float64, 1, 'upper_bounds')
        objective = _require_exact_array(
            objective_coefficients, np.float64, 1, 'objective_coefficients'
        )
        if reaction_offsets.size < 2 or reaction_offsets[0] != 0:
            raise ValueError('model_reaction_offsets must start at zero and contain a model')
        model_count = reaction_offsets.size - 1
        if metabolite_offsets.shape != (model_count + 1,) or metabolite_offsets[0] != 0:
            raise ValueError('model_metabolite_offsets must have shape (M+1,) and start at zero')
        exchange = _require_exact_array(
            exchange_reaction_indices,
            np.int32,
            2,
            'exchange_reaction_indices',
            (model_count, 3),
        )
        biomass = _require_exact_array(
            biomass_reaction_indices,
            np.int32,
            1,
            'biomass_reaction_indices',
            (model_count,),
        )
        acetate = _require_exact_array(
            acetate_reaction_indices,
            np.int32,
            1,
            'acetate_reaction_indices',
            (model_count,),
        )
        reaction_count = int(reaction_offsets[-1])
        metabolite_count = int(metabolite_offsets[-1])
        if reaction_count <= 0 or metabolite_count <= 0:
            raise ValueError('model offsets must describe nonempty models')
        if (
            np.any(reaction_offsets[1:] < reaction_offsets[:-1])
            or np.any(metabolite_offsets[1:] < metabolite_offsets[:-1])
            or lower.shape != (reaction_count,)
            or upper.shape != (reaction_count,)
            or objective.shape != (reaction_count,)
            or row_offsets.shape != (metabolite_count + 1,)
            or row_offsets[0] != 0
            or row_offsets[-1] != values.size
            or column_indices.shape != values.shape
            or np.any(row_offsets[1:] < row_offsets[:-1])
            or np.any(column_indices < 0)
            or np.any(column_indices >= reaction_count)
        ):
            raise ValueError('model arrays have inconsistent extents or indices')
        if (
            not np.all(np.isfinite(values))
            or not np.all(np.isfinite(lower))
            or not np.all(np.isfinite(upper))
            or not np.all(np.isfinite(objective))
            or np.any(lower > upper)
        ):
            raise ValueError('model coefficients and finite bounds must be valid')
        tolerance = _finite_float(cache_tolerance_molar, 'cache_tolerance_molar')
        if tolerance <= 0.0:
            raise ValueError('cache_tolerance_molar must be positive')
        if isinstance(cache_capacity, bool) or not isinstance(cache_capacity, (int, np.integer)):
            raise TypeError('cache_capacity must be an integer')
        capacity = int(cache_capacity)
        if capacity <= 0 or capacity > MAX_CACHE_CAPACITY:
            raise ValueError(f'cache_capacity must be in [1, {MAX_CACHE_CAPACITY}]')

        self._reaction_count = reaction_count
        self._cache_tolerance = tolerance
        self._cache_capacity = capacity
        self._dense_models = []
        for model_index in range(model_count):
            reaction_start = int(reaction_offsets[model_index])
            reaction_stop = int(reaction_offsets[model_index + 1])
            metabolite_start = int(metabolite_offsets[model_index])
            metabolite_stop = int(metabolite_offsets[model_index + 1])
            dense = np.zeros(
                (metabolite_stop - metabolite_start, reaction_stop - reaction_start),
                dtype=np.float64,
            )
            for global_row in range(metabolite_start, metabolite_stop):
                for entry in range(int(row_offsets[global_row]), int(row_offsets[global_row + 1])):
                    global_column = int(column_indices[entry])
                    if global_column < reaction_start or global_column >= reaction_stop:
                        raise ValueError('stoichiometry columns must remain within each model')
                    dense[global_row - metabolite_start, global_column - reaction_start] = values[entry]
            exchange_values = tuple(int(index) for index in exchange[model_index])
            marker_values = (int(biomass[model_index]), int(acetate[model_index]))
            if any(
                index < reaction_start or index >= reaction_stop
                for index in exchange_values + marker_values
            ):
                raise ValueError('reaction marker indices must remain within each model')
            self._dense_models.append(_DenseLPModel(
                reaction_start=reaction_start,
                reaction_stop=reaction_stop,
                stoichiometry=dense,
                lower_bounds=np.array(lower[reaction_start:reaction_stop], copy=True),
                upper_bounds=np.array(upper[reaction_start:reaction_stop], copy=True),
                objective=np.array(objective[reaction_start:reaction_stop], copy=True),
                exchange_indices=exchange_values,
                biomass_index=marker_values[0],
                acetate_index=marker_values[1],
            ))
        self._model_indices: Optional[np.ndarray] = None
        self._uptake_vmax: Optional[np.ndarray] = None
        self._uptake_km: Optional[np.ndarray] = None
        self._inhibition_multiplier: Optional[np.ndarray] = None
        self._model_revision = 0
        self._prior_revision = 0
        self._cache: dict = {}
        self.clear_cache()

    def configure_phenotypes(
        self,
        model_indices: np.ndarray,
        uptake_vmax: np.ndarray,
        uptake_km: np.ndarray,
        inhibition_alpha: np.ndarray,
        inhibition_beta: np.ndarray,
        max_inhibition: np.ndarray,
        model_revision: int,
        prior_revision: int,
    ) -> None:
        models = _require_exact_array(model_indices, np.int32, 1, 'model_indices')
        phenotype_count = models.size
        if phenotype_count == 0 or np.any(models < 0) or np.any(models >= len(self._dense_models)):
            raise ValueError('model_indices must select at least one configured model')
        parameters = (
            _require_exact_array(uptake_vmax, np.float64, 2, 'uptake_vmax', (phenotype_count, 3)),
            _require_exact_array(uptake_km, np.float64, 2, 'uptake_km', (phenotype_count, 3)),
            _require_exact_array(inhibition_alpha, np.float64, 2, 'inhibition_alpha', (phenotype_count, 3)),
            _require_exact_array(inhibition_beta, np.float64, 2, 'inhibition_beta', (phenotype_count, 3)),
            _require_exact_array(max_inhibition, np.float64, 2, 'max_inhibition', (phenotype_count, 3)),
        )
        vmax, km, alpha, beta, inhibition = parameters
        if (
            not all(np.all(np.isfinite(parameter)) for parameter in parameters)
            or np.any(vmax < 0.0)
            or np.any(km <= 0.0)
            or np.any(alpha <= 0.0)
            or np.any(beta <= 0.0)
            or np.any(inhibition < 0.0)
            or np.any(inhibition > 1.0)
            or not np.all(np.isfinite(alpha + beta))
        ):
            raise ValueError('phenotype parameters violate their finite ranges')
        for revision, name in ((model_revision, 'model_revision'), (prior_revision, 'prior_revision')):
            if isinstance(revision, bool) or not isinstance(revision, (int, np.integer)):
                raise TypeError(f'{name} must be an integer')
            if int(revision) < 0 or int(revision) > np.iinfo(np.uint64).max:
                raise ValueError(f'{name} must fit uint64')
        self._model_indices = np.array(models, copy=True)
        self._uptake_vmax = np.array(vmax, copy=True)
        self._uptake_km = np.array(km, copy=True)
        # The posterior mean active fraction applies deterministic inhibition without
        # sampling or adding biological constants; see DOCUMENTATION.md §10.3
        # (Posterior-Mean Exchange Bounds).
        self._inhibition_multiplier = 1.0 - alpha / (alpha + beta) * inhibition
        self._model_revision = int(model_revision)
        self._prior_revision = int(prior_revision)
        self.clear_cache()

    def clear_cache(self) -> None:
        self._cache.clear()
        self._hits = 0
        self._misses = 0
        self._evictions = 0
        self._lp_solves = 0
        self._sequence = 0

    def cache_stats(self) -> Tuple[int, int, int, int, int]:
        return len(self._cache), self._hits, self._misses, self._evictions, self._lp_solves

    def snapshot_cache(self) -> Tuple[dict, int, int, int, int, int]:
        """Return an opaque deep copy of cache payload, counters, and LRU sequence."""
        return (
            copy.deepcopy(self._cache),
            self._hits,
            self._misses,
            self._evictions,
            self._lp_solves,
            self._sequence,
        )

    def restore_cache(self, snapshot: Tuple[dict, int, int, int, int, int]) -> None:
        """Restore a snapshot produced by this kernel implementation."""
        if not isinstance(snapshot, tuple) or len(snapshot) != 6:
            raise TypeError('cache snapshot is invalid')
        cache, hits, misses, evictions, lp_solves, sequence = snapshot
        if not isinstance(cache, dict):
            raise TypeError('cache snapshot payload is invalid')
        counters = (hits, misses, evictions, lp_solves, sequence)
        if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in counters):
            raise ValueError('cache snapshot counters must be nonnegative integers')
        self._cache = copy.deepcopy(cache)
        self._hits = hits
        self._misses = misses
        self._evictions = evictions
        self._lp_solves = lp_solves
        self._sequence = sequence

    def _prepare_solve(
        self,
        phenotype_index: int,
        concentrations_m: np.ndarray,
    ) -> Tuple[int, _DenseLPModel, np.ndarray, np.ndarray, Tuple[int, int, int]]:
        if self._model_indices is None:
            raise RuntimeError('configure_phenotypes must be called before solving')
        if isinstance(phenotype_index, bool) or not isinstance(phenotype_index, (int, np.integer)):
            raise TypeError('phenotype_index must be an integer')
        phenotype = int(phenotype_index)
        if phenotype < 0 or phenotype >= self._model_indices.size:
            raise ValueError('phenotype_index must select a configured phenotype')
        concentrations = _require_exact_array(
            concentrations_m, np.float64, 1, 'concentrations_m', (3,)
        )
        if not np.all(np.isfinite(concentrations)) or np.any(concentrations < 0.0):
            raise ValueError('concentrations_m must be finite and nonnegative')
        model = self._dense_models[int(self._model_indices[phenotype])]
        vmax = self._uptake_vmax[phenotype]
        km = self._uptake_km[phenotype]
        multiplier = self._inhibition_multiplier[phenotype]
        # Saturating phenotype capacity is capped by each finite model exchange bound;
        # zero concentration therefore has exactly zero uptake; see future
        # DOCUMENTATION.md §10.3 (Posterior-Mean Exchange Bounds).
        dynamic = vmax * concentrations / (km + concentrations) * multiplier
        model_exchange_upper = np.asarray(
            [model.upper_bounds[index - model.reaction_start] for index in model.exchange_indices],
            dtype=np.float64,
        )
        effective_bounds = np.minimum(model_exchange_upper, dynamic)
        bins = quantize_cache_key(concentrations, self._cache_tolerance)
        return phenotype, model, np.array(concentrations, copy=True), effective_bounds, bins

    def _cache_key(self, phenotype: int, bins: Tuple[int, int, int]) -> Tuple[Any, ...]:
        return self._model_revision, self._prior_revision, phenotype, *bins

    def _lookup_cache(
        self,
        phenotype: int,
        concentrations: np.ndarray,
        effective_bounds: np.ndarray,
        bins: Tuple[int, int, int],
    ) -> Optional[np.ndarray]:
        matches = []
        for neighbor in enumerate_neighbor_cache_keys(bins):
            entry = self._cache.get(self._cache_key(phenotype, neighbor))
            if entry is None:
                continue
            if (
                np.max(np.abs(entry.concentrations - concentrations), initial=0.0)
                <= self._cache_tolerance
            ):
                matches.append(entry)
        if not matches:
            return None
        entry = min(matches, key=lambda candidate: candidate.insertion_sequence)
        self._sequence += 1
        entry.last_used_sequence = self._sequence
        return np.array(entry.flux_values, copy=True)

    def _insert_cache(
        self,
        phenotype: int,
        concentrations: np.ndarray,
        effective_bounds: np.ndarray,
        bins: Tuple[int, int, int],
        flux_values: np.ndarray,
    ) -> None:
        key = self._cache_key(phenotype, bins)
        if key not in self._cache and len(self._cache) >= self._cache_capacity:
            eviction_key = min(
                self._cache,
                key=lambda candidate_key: (
                    self._cache[candidate_key].last_used_sequence,
                    self._cache[candidate_key].insertion_sequence,
                ),
            )
            del self._cache[eviction_key]
            self._evictions += 1
        self._sequence += 1
        self._cache[key] = _CacheEntry(
            concentrations=np.array(concentrations, copy=True),
            effective_exchange_bounds=np.array(effective_bounds, copy=True),
            flux_values=np.array(flux_values, copy=True),
            insertion_sequence=self._sequence,
            last_used_sequence=self._sequence,
        )

    def _solve_prepared(
        self,
        model: _DenseLPModel,
        effective_bounds: np.ndarray,
    ) -> np.ndarray:
        upper = np.array(model.upper_bounds, copy=True)
        for substrate, global_index in enumerate(model.exchange_indices):
            upper[global_index - model.reaction_start] = effective_bounds[substrate]
        status, local_flux = _solve_bounded_equality_lp(
            model.stoichiometry,
            model.lower_bounds,
            upper,
            model.objective,
        )
        if status is not _LPStatus.OPTIMAL:
            raise MetabolicSolverError(f'metabolic LP status: {status.value}')
        flux_values = np.zeros(self._reaction_count, dtype=np.float64)
        flux_values[model.reaction_start:model.reaction_stop] = local_flux
        return flux_values

    def solve_flux(
        self,
        phenotype_index: int,
        concentrations_m: np.ndarray,
        cache_enabled: bool = True,
    ) -> Tuple[np.ndarray, bool]:
        if not isinstance(cache_enabled, (bool, np.bool_)):
            raise TypeError('cache_enabled must be a bool')
        phenotype, model, concentrations, effective_bounds, bins = self._prepare_solve(
            phenotype_index,
            concentrations_m,
        )
        if cache_enabled:
            cached = self._lookup_cache(
                phenotype, concentrations, effective_bounds, bins
            )
            if cached is not None:
                self._hits += 1
                return cached, True
            self._misses += 1
        flux_values = self._solve_prepared(model, effective_bounds)
        self._lp_solves += 1
        if cache_enabled:
            self._insert_cache(
                phenotype,
                concentrations,
                effective_bounds,
                bins,
                flux_values,
            )
        return flux_values, False

    def _resolve_batch_fluxes(
        self,
        phenotype_id: np.ndarray,
        environment_id: np.ndarray,
        concentrations_m: np.ndarray,
        cache_enabled: bool,
    ) -> np.ndarray:
        agent_count = phenotype_id.size
        fluxes = np.zeros((agent_count, self._reaction_count), dtype=np.float64)
        prepared = [
            self._prepare_solve(
                int(phenotype_id[agent]),
                concentrations_m[int(environment_id[agent])],
            )
            for agent in range(agent_count)
        ]
        unique_requests = []
        request_for_agent = np.full(agent_count, -1, dtype=np.int64)
        for agent, request in enumerate(prepared):
            phenotype, _, concentrations, effective_bounds, bins = request
            if cache_enabled:
                cached = self._lookup_cache(
                    phenotype,
                    concentrations,
                    effective_bounds,
                    bins,
                )
                if cached is not None:
                    self._hits += 1
                    fluxes[agent] = cached
                    continue
                self._misses += 1
                duplicate = next((
                    index
                    for index, candidate in enumerate(unique_requests)
                    if candidate[0] == phenotype
                    and np.array_equal(candidate[2], concentrations)
                    and np.array_equal(candidate[3], effective_bounds)
                ), None)
                if duplicate is not None:
                    request_for_agent[agent] = duplicate
                    continue
            request_for_agent[agent] = len(unique_requests)
            unique_requests.append(request)

        solved = []
        for request in unique_requests:
            phenotype, model, concentrations, effective_bounds, bins = request
            flux_values = self._solve_prepared(model, effective_bounds)
            self._lp_solves += 1
            if cache_enabled:
                self._insert_cache(
                    phenotype,
                    concentrations,
                    effective_bounds,
                    bins,
                    flux_values,
                )
            solved.append(flux_values)
        for agent, request_index in enumerate(request_for_agent):
            if request_index >= 0:
                fluxes[agent] = solved[int(request_index)]
        return fluxes

    def advance_batch(
        self,
        phenotype_id: np.ndarray,
        environment_id: np.ndarray,
        concentrations_m: np.ndarray,
        environment_volume_l: np.ndarray,
        biomass_gdw: np.ndarray,
        growth_rate_s: np.ndarray,
        biomass_delta_gdw: np.ndarray,
        substrate_consumed_mol: np.ndarray,
        acetate_secreted_mol: np.ndarray,
        dt_seconds: float,
        *,
        low_substrate_threshold_m: float = 1.0e-4,
        depletion_safety: float,
        max_substeps: int,
        cache_enabled: bool,
    ) -> Tuple[int, int, int, int, int, float]:
        if self._model_indices is None:
            raise RuntimeError('configure_phenotypes must be called before advancing')
        phenotypes = _require_exact_array(phenotype_id, np.int32, 1, 'phenotype_id')
        environments = _require_exact_array(environment_id, np.int32, 1, 'environment_id')
        agent_count = phenotypes.size
        if environments.shape != (agent_count,):
            raise ValueError('environment_id must have shape (N,)')
        concentrations = _require_exact_array(
            concentrations_m, np.float64, 2, 'concentrations_m', writable=True
        )
        if concentrations.shape[1:] != (3,):
            raise ValueError('concentrations_m must have shape (E, 3)')
        environment_count = concentrations.shape[0]
        volumes = _require_exact_array(
            environment_volume_l,
            np.float64,
            1,
            'environment_volume_l',
            (environment_count,),
        )
        biomass = _require_exact_array(
            biomass_gdw, np.float64, 1, 'biomass_gdw', (agent_count,), writable=True
        )
        growth = _require_exact_array(
            growth_rate_s, np.float64, 1, 'growth_rate_s', (agent_count,), writable=True
        )
        biomass_delta = _require_exact_array(
            biomass_delta_gdw,
            np.float64,
            1,
            'biomass_delta_gdw',
            (agent_count,),
            writable=True,
        )
        substrate_consumed = _require_exact_array(
            substrate_consumed_mol,
            np.float64,
            2,
            'substrate_consumed_mol',
            (agent_count, 3),
            writable=True,
        )
        acetate_secreted = _require_exact_array(
            acetate_secreted_mol,
            np.float64,
            1,
            'acetate_secreted_mol',
            (agent_count,),
            writable=True,
        )
        all_arrays = (
            ('phenotype_id', phenotypes, False),
            ('environment_id', environments, False),
            ('concentrations_m', concentrations, True),
            ('environment_volume_l', volumes, False),
            ('biomass_gdw', biomass, True),
            ('growth_rate_s', growth, True),
            ('biomass_delta_gdw', biomass_delta, True),
            ('substrate_consumed_mol', substrate_consumed, True),
            ('acetate_secreted_mol', acetate_secreted, True),
        )
        for left_index, (left_name, left, left_mutable) in enumerate(all_arrays):
            for right_name, right, right_mutable in all_arrays[left_index + 1:]:
                if (left_mutable or right_mutable) and np.shares_memory(left, right):
                    raise ValueError(f'{left_name} and {right_name} must not overlap')
        if (
            np.any(phenotypes < 0)
            or np.any(phenotypes >= self._model_indices.size)
            or np.any(environments < 0)
            or np.any(environments >= environment_count)
        ):
            raise ValueError('phenotype_id and environment_id entries must be in range')
        if (
            not np.all(np.isfinite(concentrations))
            or np.any(concentrations < 0.0)
            or not np.all(np.isfinite(volumes))
            or np.any(volumes <= 0.0)
            or not np.all(np.isfinite(biomass))
            or np.any(biomass < 0.0)
        ):
            raise ValueError('concentrations, volumes, and biomass violate finite ranges')
        if not all(np.all(np.isfinite(array)) for array in (
            growth,
            biomass_delta,
            substrate_consumed,
            acetate_secreted,
        )):
            raise ValueError('mutable metabolic output arrays must contain finite values')
        elapsed_total = _finite_float(dt_seconds, 'dt_seconds')
        threshold = _finite_float(low_substrate_threshold_m, 'low_substrate_threshold_m')
        safety = _finite_float(depletion_safety, 'depletion_safety')
        if elapsed_total < 0.0 or threshold < 0.0:
            raise ValueError('dt_seconds and low_substrate_threshold_m must be nonnegative')
        if safety <= 0.0 or safety > 1.0:
            raise ValueError('depletion_safety must be in (0, 1]')
        if isinstance(max_substeps, bool) or not isinstance(max_substeps, (int, np.integer)):
            raise TypeError('max_substeps must be an integer')
        if int(max_substeps) <= 0:
            raise ValueError('max_substeps must be positive')
        if not isinstance(cache_enabled, (bool, np.bool_)):
            raise TypeError('cache_enabled must be a bool')

        mutable_arrays = (
            concentrations,
            biomass,
            growth,
            biomass_delta,
            substrate_consumed,
            acetate_secreted,
        )
        snapshots = tuple(np.array(array, copy=True) for array in mutable_arrays)
        cache_snapshot = copy.deepcopy(self._cache)
        starting_stats = self.cache_stats()
        starting_sequence = self._sequence
        try:
            growth.fill(0.0)
            biomass_delta.fill(0.0)
            substrate_consumed.fill(0.0)
            acetate_secreted.fill(0.0)
            if elapsed_total == 0.0:
                return 0, 0, 0, 0, 0, 0.0

            growth_integral = np.zeros(agent_count, dtype=np.float64)
            elapsed = 0.0
            substeps = 0
            minimum_substep = math.inf
            while elapsed < elapsed_total:
                if substeps >= int(max_substeps):
                    raise RuntimeError('max_substeps exhausted before dt_seconds completed')
                fluxes = self._resolve_batch_fluxes(
                    phenotypes,
                    environments,
                    concentrations,
                    bool(cache_enabled),
                )
                concentration_rates = np.zeros((environment_count, 3), dtype=np.float64)
                growth_fluxes = np.zeros(agent_count, dtype=np.float64)
                acetate_fluxes = np.zeros(agent_count, dtype=np.float64)
                uptake_fluxes = np.zeros((agent_count, 3), dtype=np.float64)
                for agent in range(agent_count):
                    phenotype = int(phenotypes[agent])
                    model = self._dense_models[int(self._model_indices[phenotype])]
                    environment = int(environments[agent])
                    for substrate, reaction_index in enumerate(model.exchange_indices):
                        uptake = fluxes[agent, reaction_index]
                        uptake_fluxes[agent, substrate] = uptake
                        # Specific uptake times dry biomass divided by liters is the
                        # deterministic environmental depletion rate in M/s; see future
                        # DOCUMENTATION.md §10.3 (Adaptive Metabolic Integration).
                        concentration_rates[environment, substrate] += (
                            uptake * biomass[agent] / volumes[environment]
                        )
                    growth_fluxes[agent] = fluxes[agent, model.biomass_index]
                    acetate_fluxes[agent] = fluxes[agent, model.acetate_index]

                remaining = elapsed_total - elapsed
                substep = remaining
                for environment in range(environment_count):
                    for substrate in range(3):
                        rate = concentration_rates[environment, substrate]
                        if rate <= 0.0:
                            continue
                        # The concentration/rate depletion time prevents overshoot; the
                        # safety factor only tightens glycan/MDCF steps below the reviewed
                        # low-substrate threshold; see DOCUMENTATION.md §10.3
                        # (Adaptive Metabolic Integration).
                        depletion_limit = concentrations[environment, substrate] / rate
                        if substrate in (1, 2) and concentrations[environment, substrate] < threshold:
                            depletion_limit *= safety
                        substep = min(substep, depletion_limit)
                if not math.isfinite(substep) or substep <= 0.0 or elapsed + substep == elapsed:
                    raise RuntimeError('adaptive metabolic substep made no numerical progress')
                minimum_substep = min(minimum_substep, substep)

                step_consumption = uptake_fluxes * biomass[:, None] * substep
                for environment in range(environment_count):
                    members = np.flatnonzero(environments == environment)
                    for substrate in range(3):
                        requested = float(np.sum(step_consumption[members, substrate]))
                        available = concentrations[environment, substrate] * volumes[environment]
                        if requested > available and requested > 0.0:
                            # Proportional scaling is a deterministic roundoff backstop that
                            # conserves available substrate across all agents sharing an
                            # environment; see DOCUMENTATION.md §10.3 (Adaptive
                            # Metabolic Integration).
                            step_consumption[members, substrate] *= available / requested

                step_biomass_delta = growth_fluxes * biomass * substep
                step_acetate = acetate_fluxes * biomass * substep
                for environment in range(environment_count):
                    members = np.flatnonzero(environments == environment)
                    totals = np.sum(step_consumption[members], axis=0)
                    concentrations[environment] = np.maximum(
                        0.0,
                        concentrations[environment] - totals / volumes[environment],
                    )
                substrate_consumed += step_consumption
                biomass_delta += step_biomass_delta
                acetate_secreted += step_acetate
                growth_integral += growth_fluxes * substep
                biomass += step_biomass_delta
                elapsed += substep
                substeps += 1
            growth[:] = growth_integral / elapsed_total
            ending_stats = self.cache_stats()
            return (
                substeps,
                ending_stats[4] - starting_stats[4],
                ending_stats[1] - starting_stats[1],
                ending_stats[2] - starting_stats[2],
                ending_stats[3] - starting_stats[3],
                minimum_substep,
            )
        except BaseException:
            for array, snapshot in zip(mutable_arrays, snapshots):
                array[...] = snapshot
            self._cache = cache_snapshot
            _, self._hits, self._misses, self._evictions, self._lp_solves = starting_stats
            self._sequence = starting_sequence
            raise


_PHENOTYPE_CONFIG_KEYS = frozenset((
    'model_indices',
    'uptake_vmax',
    'uptake_km',
    'inhibition_alpha',
    'inhibition_beta',
    'max_inhibition',
    'model_revision',
    'prior_revision',
    'cache_capacity',
))


def build_neutral_phenotype_config(
    uptake_vmax: np.ndarray,
    uptake_km: np.ndarray,
    model: Optional[MetabolicModelData] = None,
    priors: Optional[MetabolicPriorData] = None,
    cache_capacity: int = DEFAULT_CACHE_CAPACITY,
) -> Mapping[str, Any]:
    model_data = load_metabolic_model() if model is None else model
    prior_data = load_metabolic_priors() if priors is None else priors
    vmax = _require_exact_array(uptake_vmax, np.float64, 2, 'uptake_vmax', (1, 3))
    km = _require_exact_array(uptake_km, np.float64, 2, 'uptake_km', (1, 3))
    if tuple(row.substrate for row in prior_data.rows) != SUBSTRATE_ORDER:
        raise ValueError('neutral prior rows must follow the fixed substrate order')
    return MappingProxyType({
        'model_indices': np.asarray((0,), dtype=np.int32),
        'uptake_vmax': np.array(vmax, copy=True, order='C'),
        'uptake_km': np.array(km, copy=True, order='C'),
        'inhibition_alpha': np.asarray(
            (tuple(row.alpha for row in prior_data.rows),), dtype=np.float64
        ),
        'inhibition_beta': np.asarray(
            (tuple(row.beta for row in prior_data.rows),), dtype=np.float64
        ),
        'max_inhibition': np.asarray(
            (tuple(row.max_inhibition for row in prior_data.rows),), dtype=np.float64
        ),
        'model_revision': model_data.model_revision,
        'prior_revision': prior_data.prior_revision,
        'cache_capacity': cache_capacity,
    })


def build_species_phenotype_config(
    species_profiles: Sequence[Any],
    model: Optional[MetabolicModelData] = None,
    priors: Optional[MetabolicPriorData] = None,
    cache_capacity: int = DEFAULT_CACHE_CAPACITY,
) -> Mapping[str, Any]:
    """Build the reduced-model phenotype matrix from legacy species affinities.

    ``aff_*`` values are dimensionless multipliers of substrate-specific base
    capacities in mol glucose-equivalent / gDW / s. ``Km`` remains an independent
    positive molar concentration. The numerical calibration is Own (BacSandbox
    model), not an empirical whole-cell kinetic measurement.
    """
    if isinstance(species_profiles, (str, bytes)) or not isinstance(species_profiles, Sequence):
        raise TypeError('species_profiles must be a sequence of mappings or species objects')
    profiles = tuple(species_profiles)
    if not profiles:
        raise ValueError('species_profiles must contain at least one species')

    model_data = load_metabolic_model() if model is None else model
    prior_data = load_metabolic_priors() if priors is None else priors
    if tuple(row.substrate for row in prior_data.rows) != SUBSTRATE_ORDER:
        raise ValueError('species prior rows must follow the fixed substrate order')

    affinities = np.empty((len(profiles), len(SUBSTRATE_ORDER)), dtype=np.float64)
    for profile_index, profile in enumerate(profiles):
        for substrate_index, (field_name, default_value) in enumerate(zip(
            SPECIES_AFFINITY_FIELDS,
            SPECIES_DEFAULT_AFFINITY_MULTIPLIERS,
        )):
            if isinstance(profile, Mapping):
                raw_value = profile.get(field_name, default_value)
            else:
                raw_value = getattr(profile, field_name, default_value)
            if raw_value is None or (isinstance(raw_value, str) and not raw_value.strip()):
                raw_value = default_value
            if isinstance(raw_value, bool):
                raise TypeError(f'species_profiles[{profile_index}].{field_name} must be numeric')
            try:
                affinity = float(raw_value)
            except (TypeError, ValueError) as exc:
                raise TypeError(
                    f'species_profiles[{profile_index}].{field_name} must be numeric'
                ) from exc
            if not math.isfinite(affinity) or affinity < 0.0:
                raise ValueError(
                    f'species_profiles[{profile_index}].{field_name} must be finite and nonnegative'
                )
            affinities[profile_index, substrate_index] = affinity

    base_vmax = np.asarray(SPECIES_BASE_UPTAKE_VMAX_MOL_GDW_S, dtype=np.float64)
    # A legacy affinity scales capacity, so zero produces an exact zero exchange
    # upper bound without changing Michaelis-Menten Km; see DOCUMENTATION.md section 10.4.
    uptake_vmax = affinities * base_vmax
    if not np.all(np.isfinite(uptake_vmax)):
        raise ValueError('species affinity scaling produced nonfinite uptake_vmax')
    # Every phenotype shares the positive Own-model substrate Km calibration while
    # retaining an independent writable row; see DOCUMENTATION.md section 10.4.
    uptake_km = np.tile(
        np.asarray(SPECIES_UPTAKE_KM_MOLAR, dtype=np.float64),
        (len(profiles), 1),
    )
    prior_alpha = np.asarray(tuple(row.alpha for row in prior_data.rows), dtype=np.float64)
    prior_beta = np.asarray(tuple(row.beta for row in prior_data.rows), dtype=np.float64)
    prior_inhibition = np.asarray(
        tuple(row.max_inhibition for row in prior_data.rows), dtype=np.float64
    )
    # Every species selects the shared reduced bookkeeping model at index zero, and
    # neutral prior rows repeat per phenotype; see DOCUMENTATION.md section 10.4.
    return MappingProxyType({
        'model_indices': np.zeros(len(profiles), dtype=np.int32),
        'uptake_vmax': np.array(uptake_vmax, copy=True, order='C'),
        'uptake_km': np.array(uptake_km, copy=True, order='C'),
        'inhibition_alpha': np.tile(prior_alpha, (len(profiles), 1)),
        'inhibition_beta': np.tile(prior_beta, (len(profiles), 1)),
        'max_inhibition': np.tile(prior_inhibition, (len(profiles), 1)),
        'model_revision': model_data.model_revision,
        'prior_revision': prior_data.prior_revision,
        'cache_capacity': cache_capacity,
    })


def _validate_phenotype_config(
    phenotype_config: Mapping[str, Any],
    model: MetabolicModelData,
    priors: MetabolicPriorData,
) -> Mapping[str, Any]:
    config = _require_exact_keys(phenotype_config, _PHENOTYPE_CONFIG_KEYS, 'phenotype_config')
    models = _require_exact_array(config['model_indices'], np.int32, 1, 'model_indices')
    phenotype_count = models.size
    if phenotype_count == 0:
        raise ValueError('phenotype_config must contain at least one phenotype')
    for name in (
        'uptake_vmax',
        'uptake_km',
        'inhibition_alpha',
        'inhibition_beta',
        'max_inhibition',
    ):
        _require_exact_array(config[name], np.float64, 2, name, (phenotype_count, 3))
    if config['model_revision'] != model.model_revision:
        raise ValueError('phenotype_config model_revision does not match the model')
    if config['prior_revision'] != priors.prior_revision:
        raise ValueError('phenotype_config prior_revision does not match the priors')
    capacity = config['cache_capacity']
    if isinstance(capacity, bool) or not isinstance(capacity, (int, np.integer)):
        raise TypeError('cache_capacity must be an integer')
    if int(capacity) <= 0 or int(capacity) > MAX_CACHE_CAPACITY:
        raise ValueError(f'cache_capacity must be in [1, {MAX_CACHE_CAPACITY}]')
    return MappingProxyType(dict(config))


class MetabolicOrchestrator:
    """Lazy legacy/shadow/FBA boundary for the Phase C rollout."""

    def __init__(
        self,
        mode: Union[str, MetabolicMode] = MetabolicMode.LEGACY,
        backend: Union[str, MetabolicBackend] = MetabolicBackend.AUTO,
        model: Optional[MetabolicModelData] = None,
        priors: Optional[MetabolicPriorData] = None,
        phenotype_config: Optional[Mapping[str, Any]] = None,
    ) -> None:
        self._mode = _coerce_enum(mode, MetabolicMode, 'mode')
        self._backend = _coerce_enum(backend, MetabolicBackend, 'backend')
        if model is not None and not isinstance(model, MetabolicModelData):
            raise TypeError('model must be MetabolicModelData')
        if priors is not None and not isinstance(priors, MetabolicPriorData):
            raise TypeError('priors must be MetabolicPriorData')
        if phenotype_config is not None and not isinstance(phenotype_config, Mapping):
            raise TypeError('phenotype_config must be a mapping')
        if self._mode is not MetabolicMode.LEGACY:
            missing = [
                name
                for name, value in (
                    ('model', model),
                    ('priors', priors),
                    ('phenotype_config', phenotype_config),
                )
                if value is None
            ]
            if missing:
                raise ValueError(
                    f'{self._mode.value} mode requires explicit ' + ', '.join(missing)
                )
        self._model = model
        self._priors = priors
        self._phenotype_config = None
        if phenotype_config is not None:
            if model is None or priors is None:
                raise ValueError('phenotype_config requires explicit model and priors')
            self._phenotype_config = _validate_phenotype_config(
                phenotype_config,
                model,
                priors,
            )
        self._native_available: Optional[bool] = None
        self._backend_used: Optional[str] = None
        self._kernel: Optional[Any] = None

    @property
    def mode(self) -> str:
        return self._mode.value

    @property
    def backend(self) -> str:
        return self._backend.value

    @property
    def backend_used(self) -> Optional[str]:
        return self._backend_used

    @property
    def model_revision(self) -> Optional[int]:
        return None if self._model is None else self._model.model_revision

    @property
    def prior_revision(self) -> Optional[int]:
        return None if self._priors is None else self._priors.prior_revision

    @property
    def native_available(self) -> Optional[bool]:
        if self._mode is MetabolicMode.LEGACY:
            return None
        if self._backend is MetabolicBackend.PYTHON:
            return False
        return self._native_available

    def _constructor_arguments(self) -> Mapping[str, Any]:
        return {
            'model_reaction_offsets': self._model.model_reaction_offsets,
            'model_metabolite_offsets': self._model.model_metabolite_offsets,
            'stoich_row_offsets': self._model.stoich_row_offsets,
            'stoich_column_indices': self._model.stoich_column_indices,
            'stoich_values': self._model.stoich_values,
            'lower_bounds': self._model.lower_bounds,
            'upper_bounds': self._model.upper_bounds,
            'objective_coefficients': self._model.objective_coefficients,
            'exchange_reaction_indices': self._model.exchange_reaction_indices,
            'biomass_reaction_indices': self._model.biomass_reaction_indices,
            'acetate_reaction_indices': self._model.acetate_reaction_indices,
            'cache_tolerance_molar': self._model.cache_tolerance_molar,
            'cache_capacity': int(self._phenotype_config['cache_capacity']),
        }

    def _configure_kernel(self, kernel: Any) -> None:
        config = self._phenotype_config
        kernel.configure_phenotypes(
            config['model_indices'],
            config['uptake_vmax'],
            config['uptake_km'],
            config['inhibition_alpha'],
            config['inhibition_beta'],
            config['max_inhibition'],
            config['model_revision'],
            config['prior_revision'],
        )

    def _get_kernel(self) -> Any:
        if self._kernel is not None:
            return self._kernel
        kernel_type = None
        if self._backend is not MetabolicBackend.PYTHON:
            try:
                native_module = importlib.import_module('bacteria_c_kernel')
            except ModuleNotFoundError as exc:
                if exc.name != 'bacteria_c_kernel':
                    raise
                self._native_available = False
                if self._backend is MetabolicBackend.CPP:
                    raise MetabolicSolverUnavailable('native MetabolicKernel is unavailable')
            else:
                kernel_type = getattr(native_module, 'MetabolicKernel', None)
                self._native_available = kernel_type is not None
                if kernel_type is None and self._backend is MetabolicBackend.CPP:
                    raise MetabolicSolverUnavailable('native MetabolicKernel is unavailable')
        if kernel_type is None:
            kernel_type = PythonMetabolicKernel
            self._backend_used = MetabolicBackend.PYTHON.value
        else:
            self._backend_used = MetabolicBackend.CPP.value
        kernel = kernel_type(**self._constructor_arguments())
        self._configure_kernel(kernel)
        self._kernel = kernel
        return kernel

    def solve_flux(
        self,
        phenotype_index: int = 0,
        concentrations_m: Optional[np.ndarray] = None,
        cache_enabled: bool = True,
    ) -> Optional[Tuple[np.ndarray, bool]]:
        if self._mode is MetabolicMode.LEGACY:
            return None
        if concentrations_m is None:
            raise ValueError('concentrations_m is required outside legacy mode')
        return self._get_kernel().solve_flux(
            phenotype_index,
            concentrations_m,
            cache_enabled,
        )

    def clear_cache(self) -> None:
        """Clear an already-created solver cache without forcing backend initialization."""
        if self._kernel is not None:
            self._kernel.clear_cache()

    def snapshot_cache(self) -> Optional[Any]:
        """Snapshot an initialized backend cache without forcing lazy construction."""
        if self._kernel is None:
            return None
        return self._kernel.snapshot_cache()

    def restore_cache(self, snapshot: Optional[Any]) -> None:
        """Restore the initialized backend cache captured before a live transaction."""
        if snapshot is None:
            if self._kernel is not None:
                self._kernel.clear_cache()
            return
        if self._kernel is None:
            raise RuntimeError('cannot restore a cache before backend initialization')
        self._kernel.restore_cache(snapshot)

    def advance_batch(
        self,
        phenotype_id: Optional[np.ndarray] = None,
        environment_id: Optional[np.ndarray] = None,
        concentrations_m: Optional[np.ndarray] = None,
        environment_volume_l: Optional[np.ndarray] = None,
        biomass_gdw: Optional[np.ndarray] = None,
        growth_rate_s: Optional[np.ndarray] = None,
        biomass_delta_gdw: Optional[np.ndarray] = None,
        substrate_consumed_mol: Optional[np.ndarray] = None,
        acetate_secreted_mol: Optional[np.ndarray] = None,
        dt_seconds: Optional[float] = None,
        low_substrate_threshold_m: float = 1.0e-4,
        depletion_safety: float = 0.9,
        max_substeps: int = 256,
        cache_enabled: bool = True,
    ) -> Optional[Tuple[int, int, int, int, int, float]]:
        if self._mode is MetabolicMode.LEGACY:
            return None
        arguments = (
            phenotype_id,
            environment_id,
            concentrations_m,
            environment_volume_l,
            biomass_gdw,
            growth_rate_s,
            biomass_delta_gdw,
            substrate_consumed_mol,
            acetate_secreted_mol,
            dt_seconds,
        )
        if any(argument is None for argument in arguments):
            raise ValueError('all metabolic batch arguments are required outside legacy mode')
        original_concentrations = np.array(concentrations_m, copy=True)
        original_biomass = np.array(biomass_gdw, copy=True)
        summary = self._get_kernel().advance_batch(
            phenotype_id,
            environment_id,
            concentrations_m,
            environment_volume_l,
            biomass_gdw,
            growth_rate_s,
            biomass_delta_gdw,
            substrate_consumed_mol,
            acetate_secreted_mol,
            dt_seconds,
            low_substrate_threshold_m=low_substrate_threshold_m,
            depletion_safety=depletion_safety,
            max_substeps=max_substeps,
            cache_enabled=cache_enabled,
        )
        if self._mode is MetabolicMode.SHADOW:
            concentrations_m[...] = original_concentrations
            biomass_gdw[...] = original_biomass
        return summary