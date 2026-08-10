"""Schema-driven translation of dry diet mass into Phase C substrate molarity."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from numbers import Real
from pathlib import Path
from types import MappingProxyType
from typing import Any, Dict, Mapping, Sequence, Tuple, Union

import numpy as np

from .metabolism import SUBSTRATE_ORDER


DIET_SCHEMA_VERSION = 1
DIET_PROFILE_ID = 'diet_profiles_v1'
DIET_PROFILE_REVISION = 1
DIET_PROVENANCE = 'Own (BacSandbox model)'
DEFAULT_DIET_PROFILES_PATH = Path(__file__).resolve().parent / 'models' / 'diet_profiles.json'
_REQUIRED_DIET_ORDER = ('malnutrition', 'rusf', 'mdcf_2')
# This tolerance only absorbs binary representation error in source decimal sums; see
# DOCUMENTATION.md section 10.4.1.
_NUMERIC_EQUALITY_TOLERANCE = 1.0e-12

_TOP_LEVEL_KEYS = frozenset((
    'schema_version',
    'profile_id',
    'profile_revision',
    'provenance',
    'documentation_anchor',
    'scope',
    'basis',
    'constants',
    'substrate_order',
    'component_order',
    'components',
    'ingredient_order',
    'ingredients',
    'diet_order',
    'diets',
    'source_hooks',
))
_BASIS_KEYS = frozenset((
    'feed_mass',
    'ingredient_amount',
    'component_fraction',
    'volume',
    'output',
    'feed_semantics',
))
_CONSTANT_KEYS = frozenset((
    'glucose_equivalent_molar_mass_g_per_mol',
    'recipe_rounding_tolerance_g_per_100g',
))
_COMPONENT_KEYS = frozenset(('display_name', 'substrate', 'source_hook'))
_INGREDIENT_KEYS = frozenset(('display_name', 'source_hook', 'component_mass_fractions'))
_DIET_KEYS = frozenset((
    'description',
    'source_hook',
    'reported_total_g_per_100g',
    'ingredient_mass_g_per_100g',
))


class DietProfileValidationError(ValueError):
    """Raised when a diet profile asset violates its versioned schema."""


@dataclass(frozen=True)
class DietProfileData:
    """Immutable, validated arrays used by :class:`DietTranslationEngine`."""

    substrate_order: Tuple[str, ...]
    component_order: Tuple[str, ...]
    ingredient_order: Tuple[str, ...]
    diet_order: Tuple[str, ...]
    diet_index: Mapping[str, int]
    glucose_equivalent_molar_mass_g_per_mol: float
    component_to_substrate: np.ndarray
    ingredient_component_mass_fractions: np.ndarray
    diet_ingredient_mass_fractions: np.ndarray
    diet_component_mass_fractions: np.ndarray


def _reject_duplicate_keys(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    """Build one JSON object while rejecting duplicate keys at every depth."""
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DietProfileValidationError(f'duplicate JSON key {key!r}')
        result[key] = value
    return result


def _require_object(value: Any, path: str) -> Dict[str, Any]:
    """Return a decoded JSON object or raise a path-specific schema error."""
    if not isinstance(value, dict):
        raise DietProfileValidationError(f'{path} must be an object')
    return value


def _require_exact_keys(value: Mapping[str, Any], expected: frozenset, path: str) -> None:
    """Reject missing and unknown object fields."""
    keys = frozenset(value)
    missing = sorted(expected - keys)
    unknown = sorted(keys - expected)
    if missing or unknown:
        raise DietProfileValidationError(
            f'{path} has missing keys {missing} and unknown keys {unknown}'
        )


def _require_string(value: Any, path: str) -> str:
    """Return a non-empty string or raise a path-specific schema error."""
    if not isinstance(value, str) or not value.strip():
        raise DietProfileValidationError(f'{path} must be a non-empty string')
    return value


def _require_string_order(value: Any, path: str) -> Tuple[str, ...]:
    """Validate an ordered JSON string array with unique entries."""
    if not isinstance(value, list):
        raise DietProfileValidationError(f'{path} must be an array')
    result = tuple(_require_string(item, f'{path}[{index}]') for index, item in enumerate(value))
    if len(set(result)) != len(result):
        raise DietProfileValidationError(f'{path} must not contain duplicate entries')
    return result


def _require_finite_number(
    value: Any,
    path: str,
    *,
    minimum: float = 0.0,
    maximum: Union[float, None] = None,
) -> float:
    """Validate a finite JSON number within an inclusive interval."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise DietProfileValidationError(f'{path} must be a number')
    result = float(value)
    if not math.isfinite(result):
        raise DietProfileValidationError(f'{path} must be finite')
    if result < minimum or (maximum is not None and result > maximum):
        upper = 'infinity' if maximum is None else str(maximum)
        raise DietProfileValidationError(f'{path} must be in [{minimum}, {upper}]')
    return result


class DietTranslationEngine:
    """Translate versioned dry-food profiles into Phase C concentration increments."""

    def __init__(
        self,
        profiles_path: Union[str, Path, None] = None,
    ) -> None:
        """Eagerly load and validate one immutable diet profile asset."""
        path = DEFAULT_DIET_PROFILES_PATH if profiles_path is None else Path(profiles_path)
        self._data = self._load_and_validate(path)

    @property
    def substrate_order(self) -> Tuple[str, ...]:
        """Return the fixed Phase C substrate ABI order."""
        return self._data.substrate_order

    @property
    def component_order(self) -> Tuple[str, ...]:
        """Return the biochemical component order used for diagnostic output."""
        return self._data.component_order

    @property
    def diet_names(self) -> Tuple[str, ...]:
        """Return available diet identifiers in deterministic asset order."""
        return self._data.diet_order

    def get_molar_substrate_bounds(
        self,
        diet_name: str,
        feed_mass_g: float,
        volume_L: float,
        dt_s: float,
    ) -> np.ndarray:
        """Return one environment's additive three-substrate concentration increment in M."""
        diet_index, feed_mass, volume = self._validate_request(
            diet_name,
            feed_mass_g,
            volume_L,
            dt_s,
        )
        components = self._component_molarities(diet_index, feed_mass, volume)
        # Summing components through the one-hot ABI map conserves glucose-equivalent
        # monomer amount across substrate classes; see DOCUMENTATION.md section 10.4.2.
        substrate_molarities = components @ self._data.component_to_substrate
        result = np.array(substrate_molarities, dtype=np.float64, order='C', copy=True)
        if result.shape != (3,) or not np.all(np.isfinite(result)):
            raise OverflowError('translated substrate concentrations are not finite')
        return result

    def get_molar_component_concentrations(
        self,
        diet_name: str,
        feed_mass_g: float,
        volume_L: float,
        dt_s: float,
    ) -> np.ndarray:
        """Return component-level molarity in :attr:`component_order` for diagnostics."""
        diet_index, feed_mass, volume = self._validate_request(
            diet_name,
            feed_mass_g,
            volume_L,
            dt_s,
        )
        return self._component_molarities(diet_index, feed_mass, volume)

    def _component_molarities(
        self,
        diet_index: int,
        feed_mass_g: float,
        volume_L: float,
    ) -> np.ndarray:
        """Apply the documented monomer-equivalent concentration equation."""
        try:
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                # Feed mass times each class fraction gives component grams; division by
                # monomer molar mass and liters gives molarity. See DOCUMENTATION.md
                # section 10.4.2.
                scale = feed_mass_g / (
                    self._data.glucose_equivalent_molar_mass_g_per_mol * volume_L
                )
                result = np.array(
                    self._data.diet_component_mass_fractions[diet_index] * scale,
                    dtype=np.float64,
                    order='C',
                    copy=True,
                )
        except FloatingPointError as error:
            raise OverflowError('translated component concentrations are not representable') from error
        if not np.all(np.isfinite(result)):
            raise OverflowError('translated component concentrations are not finite')
        return result

    def _validate_request(
        self,
        diet_name: Any,
        feed_mass_g: Any,
        volume_L: Any,
        dt_s: Any,
    ) -> Tuple[int, float, float]:
        """Validate public request values and resolve the diet row."""
        if not isinstance(diet_name, str):
            raise TypeError('diet_name must be a string')
        if diet_name not in self._data.diet_index:
            raise ValueError(
                f'unknown diet {diet_name!r}; expected one of {self._data.diet_order}'
            )
        feed_mass = self._require_request_number(feed_mass_g, 'feed_mass_g', minimum=0.0)
        volume = self._require_request_number(volume_L, 'volume_L', minimum=0.0, strict=True)
        self._require_request_number(dt_s, 'dt_s', minimum=0.0, strict=True)
        # dt_s validates real-time scheduling, but a bolus already states total delivered
        # mass and must remain timestep invariant; see DOCUMENTATION.md section 10.4.2.
        return self._data.diet_index[diet_name], feed_mass, volume

    @staticmethod
    def _require_request_number(
        value: Any,
        name: str,
        *,
        minimum: float,
        strict: bool = False,
    ) -> float:
        """Validate a finite real-valued request scalar."""
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError(f'{name} must be a real number')
        result = float(value)
        if not math.isfinite(result):
            raise ValueError(f'{name} must be finite')
        if result < minimum or (strict and result == minimum):
            relation = 'greater than' if strict else 'greater than or equal to'
            raise ValueError(f'{name} must be {relation} {minimum}')
        return result

    @classmethod
    def _load_and_validate(cls, path: Path) -> DietProfileData:
        """Decode a profile asset and materialize strictly validated immutable arrays."""
        try:
            with path.open('r', encoding='utf-8') as profile_file:
                raw = json.load(profile_file, object_pairs_hook=_reject_duplicate_keys)
        except DietProfileValidationError:
            raise
        except (OSError, json.JSONDecodeError) as error:
            raise DietProfileValidationError(f'cannot load diet profile asset {path}: {error}') from error

        root = _require_object(raw, 'root')
        _require_exact_keys(root, _TOP_LEVEL_KEYS, 'root')
        cls._validate_identity(root)

        source_hooks = _require_object(root['source_hooks'], 'source_hooks')
        if not source_hooks:
            raise DietProfileValidationError('source_hooks must not be empty')
        for hook, target in source_hooks.items():
            _require_string(hook, 'source_hooks key')
            _require_string(target, f'source_hooks.{hook}')

        basis = _require_object(root['basis'], 'basis')
        _require_exact_keys(basis, _BASIS_KEYS, 'basis')
        for key, value in basis.items():
            _require_string(value, f'basis.{key}')
        if basis['feed_semantics'] != 'interval_total_bolus':
            raise DietProfileValidationError(
                'basis.feed_semantics must be interval_total_bolus for this API'
            )

        constants = _require_object(root['constants'], 'constants')
        _require_exact_keys(constants, _CONSTANT_KEYS, 'constants')
        molar_mass = _require_finite_number(
            constants['glucose_equivalent_molar_mass_g_per_mol'],
            'constants.glucose_equivalent_molar_mass_g_per_mol',
        )
        if molar_mass == 0.0:
            raise DietProfileValidationError(
                'constants.glucose_equivalent_molar_mass_g_per_mol must be positive'
            )
        recipe_tolerance = _require_finite_number(
            constants['recipe_rounding_tolerance_g_per_100g'],
            'constants.recipe_rounding_tolerance_g_per_100g',
        )

        substrate_order = _require_string_order(root['substrate_order'], 'substrate_order')
        if substrate_order != SUBSTRATE_ORDER:
            raise DietProfileValidationError(
                f'substrate_order must exactly match Phase C order {SUBSTRATE_ORDER}'
            )
        component_order = _require_string_order(root['component_order'], 'component_order')
        ingredient_order = _require_string_order(root['ingredient_order'], 'ingredient_order')
        diet_order = _require_string_order(root['diet_order'], 'diet_order')
        if diet_order != _REQUIRED_DIET_ORDER:
            raise DietProfileValidationError(f'diet_order must equal {_REQUIRED_DIET_ORDER}')

        component_map = cls._validate_components(
            root['components'],
            component_order,
            substrate_order,
            source_hooks,
        )
        ingredient_components = cls._validate_ingredients(
            root['ingredients'],
            ingredient_order,
            component_order,
            source_hooks,
        )
        diet_ingredients = cls._validate_diets(
            root['diets'],
            diet_order,
            ingredient_order,
            source_hooks,
            recipe_tolerance,
        )

        component_to_substrate = np.zeros(
            (len(component_order), len(substrate_order)),
            dtype=np.float64,
        )
        for component_index, substrate_index in enumerate(component_map):
            component_to_substrate[component_index, substrate_index] = 1.0
        # Linear mixing is the mass-conserving ingredient-to-diet aggregation in the
        # documented Own model; see DOCUMENTATION.md section 10.4.2.
        diet_components = diet_ingredients @ ingredient_components

        arrays = (
            component_to_substrate,
            ingredient_components,
            diet_ingredients,
            diet_components,
        )
        for array in arrays:
            if not array.flags.c_contiguous or not np.all(np.isfinite(array)):
                raise DietProfileValidationError('validated diet arrays must be finite and C-contiguous')
            array.setflags(write=False)

        return DietProfileData(
            substrate_order=substrate_order,
            component_order=component_order,
            ingredient_order=ingredient_order,
            diet_order=diet_order,
            diet_index=MappingProxyType({name: index for index, name in enumerate(diet_order)}),
            glucose_equivalent_molar_mass_g_per_mol=molar_mass,
            component_to_substrate=component_to_substrate,
            ingredient_component_mass_fractions=ingredient_components,
            diet_ingredient_mass_fractions=diet_ingredients,
            diet_component_mass_fractions=diet_components,
        )

    @staticmethod
    def _validate_identity(root: Mapping[str, Any]) -> None:
        """Validate the immutable identity and descriptive fields of the asset."""
        integer_fields = (
            ('schema_version', DIET_SCHEMA_VERSION),
            ('profile_revision', DIET_PROFILE_REVISION),
        )
        for field, expected in integer_fields:
            value = root[field]
            if isinstance(value, bool) or not isinstance(value, int) or value != expected:
                raise DietProfileValidationError(f'{field} must equal {expected}')
        if root['profile_id'] != DIET_PROFILE_ID:
            raise DietProfileValidationError(f'profile_id must equal {DIET_PROFILE_ID!r}')
        if root['provenance'] != DIET_PROVENANCE:
            raise DietProfileValidationError(f'provenance must equal {DIET_PROVENANCE!r}')
        _require_string(root['documentation_anchor'], 'documentation_anchor')
        _require_string(root['scope'], 'scope')

    @staticmethod
    def _validate_components(
        value: Any,
        component_order: Tuple[str, ...],
        substrate_order: Tuple[str, ...],
        source_hooks: Mapping[str, str],
    ) -> Tuple[int, ...]:
        """Validate component records and return their substrate indices."""
        components = _require_object(value, 'components')
        if tuple(components) != component_order:
            raise DietProfileValidationError('components keys must exactly match component_order')
        result = []
        for component in component_order:
            record = _require_object(components[component], f'components.{component}')
            _require_exact_keys(record, _COMPONENT_KEYS, f'components.{component}')
            _require_string(record['display_name'], f'components.{component}.display_name')
            substrate = _require_string(record['substrate'], f'components.{component}.substrate')
            if substrate not in substrate_order:
                raise DietProfileValidationError(
                    f'components.{component}.substrate references unknown substrate {substrate!r}'
                )
            hook = _require_string(record['source_hook'], f'components.{component}.source_hook')
            if hook not in source_hooks:
                raise DietProfileValidationError(
                    f'components.{component}.source_hook references unknown hook {hook!r}'
                )
            result.append(substrate_order.index(substrate))
        return tuple(result)

    @staticmethod
    def _validate_ingredients(
        value: Any,
        ingredient_order: Tuple[str, ...],
        component_order: Tuple[str, ...],
        source_hooks: Mapping[str, str],
    ) -> np.ndarray:
        """Validate ingredient component fractions and return an ``(I, C)`` matrix."""
        ingredients = _require_object(value, 'ingredients')
        if tuple(ingredients) != ingredient_order:
            raise DietProfileValidationError('ingredients keys must exactly match ingredient_order')
        result = np.empty((len(ingredient_order), len(component_order)), dtype=np.float64)
        for ingredient_index, ingredient in enumerate(ingredient_order):
            record = _require_object(ingredients[ingredient], f'ingredients.{ingredient}')
            _require_exact_keys(record, _INGREDIENT_KEYS, f'ingredients.{ingredient}')
            _require_string(record['display_name'], f'ingredients.{ingredient}.display_name')
            hook = _require_string(record['source_hook'], f'ingredients.{ingredient}.source_hook')
            if hook not in source_hooks:
                raise DietProfileValidationError(
                    f'ingredients.{ingredient}.source_hook references unknown hook {hook!r}'
                )
            fractions = _require_object(
                record['component_mass_fractions'],
                f'ingredients.{ingredient}.component_mass_fractions',
            )
            if tuple(fractions) != component_order:
                raise DietProfileValidationError(
                    f'ingredients.{ingredient}.component_mass_fractions keys must match component_order'
                )
            for component_index, component in enumerate(component_order):
                result[ingredient_index, component_index] = _require_finite_number(
                    fractions[component],
                    f'ingredients.{ingredient}.component_mass_fractions.{component}',
                    maximum=1.0,
                )
            if math.fsum(result[ingredient_index]) > 1.0:
                raise DietProfileValidationError(
                    f'ingredients.{ingredient}.component_mass_fractions sum above 1.0'
                )
        return result

    @staticmethod
    def _validate_diets(
        value: Any,
        diet_order: Tuple[str, ...],
        ingredient_order: Tuple[str, ...],
        source_hooks: Mapping[str, str],
        recipe_tolerance: float,
    ) -> np.ndarray:
        """Validate recipe parts and return normalized ``(D, I)`` dry-mass fractions."""
        diets = _require_object(value, 'diets')
        if tuple(diets) != diet_order:
            raise DietProfileValidationError('diets keys must exactly match diet_order')
        result = np.empty((len(diet_order), len(ingredient_order)), dtype=np.float64)
        for diet_index, diet in enumerate(diet_order):
            record = _require_object(diets[diet], f'diets.{diet}')
            _require_exact_keys(record, _DIET_KEYS, f'diets.{diet}')
            _require_string(record['description'], f'diets.{diet}.description')
            hook = _require_string(record['source_hook'], f'diets.{diet}.source_hook')
            if hook not in source_hooks:
                raise DietProfileValidationError(
                    f'diets.{diet}.source_hook references unknown hook {hook!r}'
                )
            reported_total = _require_finite_number(
                record['reported_total_g_per_100g'],
                f'diets.{diet}.reported_total_g_per_100g',
            )
            if reported_total == 0.0:
                raise DietProfileValidationError(f'diets.{diet} must contain positive dry mass')
            parts = _require_object(
                record['ingredient_mass_g_per_100g'],
                f'diets.{diet}.ingredient_mass_g_per_100g',
            )
            if tuple(parts) != ingredient_order:
                raise DietProfileValidationError(
                    f'diets.{diet}.ingredient_mass_g_per_100g keys must match ingredient_order'
                )
            raw_parts = np.empty(len(ingredient_order), dtype=np.float64)
            for ingredient_index, ingredient in enumerate(ingredient_order):
                raw_parts[ingredient_index] = _require_finite_number(
                    parts[ingredient],
                    f'diets.{diet}.ingredient_mass_g_per_100g.{ingredient}',
                )
            calculated_total = math.fsum(raw_parts)
            if not math.isclose(
                calculated_total,
                reported_total,
                rel_tol=0.0,
                abs_tol=_NUMERIC_EQUALITY_TOLERANCE,
            ):
                raise DietProfileValidationError(
                    f'diets.{diet} ingredient sum {calculated_total} does not match '
                    f'reported_total_g_per_100g {reported_total}'
                )
            if abs(reported_total - 100.0) > recipe_tolerance:
                raise DietProfileValidationError(
                    f'diets.{diet} reported total is outside source rounding tolerance'
                )
            # Published recipes round individual g/100 g entries; explicit division by
            # their reported total restores exact dry-mass conservation. See
            # DOCUMENTATION.md section 10.4.3.
            result[diet_index] = raw_parts / reported_total
        return result