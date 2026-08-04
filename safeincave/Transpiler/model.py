# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Parse a YAML case file into a validated intermediate representation.

All validation happens here, before a single line of Python is generated:
type names are resolved through :mod:`.registry`, keyword blocks are checked
against the real constructor signatures through :mod:`.signatures`, and all
numeric leaves are coerced to plain numbers (everything is SI; there is no
expression language).

Reading the file, including any ``!include`` references, is
:mod:`.include`'s job; by the time anything here runs the case is a single
plain dict.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from . import include, registry, signatures
from .errors import TranspileError
from .suggest import closest


# --------------------------------------------------------------------------- IR

@dataclass
class ObjectSpec:
    """One YAML block that instantiates a class of the safeincave API."""

    section: str
    type_name: str
    cls: type
    kwargs: dict = field(default_factory=dict)         # plain values, emitted verbatim
    tensor_kwargs: dict = field(default_factory=dict)  # scalar or {region: value} map


@dataclass
class OutputSpec:
    equation: str          # "momentum" | "heat"
    output_format: str
    folder: str
    fields: dict           # {field_name: label_name}


@dataclass
class MaterialSpec:
    name: str | None = None
    properties: dict = field(default_factory=dict)  # {setter_suffix: scalar|region map}
    elastic: list = field(default_factory=list)
    non_elastic: list = field(default_factory=list)
    thermoelastic: list = field(default_factory=list)
    sets: list = field(default_factory=lambda: ["all_elements"])  # element regions this spec covers


@dataclass
class StageSpec:
    name: str
    time: ObjectSpec
    bcs: dict              # {"momentum": [ObjectSpec], "heat": [ObjectSpec]}
    caverns: list
    logging: ObjectSpec | None
    outputs: list
    simulator: ObjectSpec


@dataclass
class CaseModel:
    source: Path
    sources: list          # every file the case was built from, root first
    grid: ObjectSpec
    equations: dict        # {"momentum": ObjectSpec, "heat": ObjectSpec}
    material: MaterialSpec
    body_force: list
    initial_temperature: float | None
    stages: list


# ------------------------------------------------------------------- primitives

def _coerce_scalar(value, context: str):
    """Coerce a YAML leaf to a number where it looks like one.

    PyYAML (YAML 1.1) parses e.g. ``1e-20`` as a string, so numeric-looking
    strings are converted to floats. Non-numeric strings are returned as-is.
    """
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return value
    if isinstance(value, list):
        return [_coerce_scalar(v, context) for v in value]
    raise TranspileError(f"{context}: unsupported value {value!r}.")


def _require_number(value, context: str) -> float:
    value = _coerce_scalar(value, context)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TranspileError(f"{context}: expected a number, got {value!r}.")
    return float(value)


def _require_mapping(node, context: str) -> dict:
    if not isinstance(node, dict):
        raise TranspileError(f"{context}: expected a mapping, got {type(node).__name__}.")
    return node


def _require_list(node, context: str) -> list:
    if not isinstance(node, list):
        raise TranspileError(f"{context}: expected a list, got {type(node).__name__}.")
    return node


def _check_keys(node: dict, allowed, context: str) -> None:
    unknown = [key for key in node if key not in allowed]
    if unknown:
        raise TranspileError(
            f"{context}: unknown key(s) {', '.join(map(repr, unknown))}. "
            f"Allowed keys: {', '.join(sorted(allowed))}."
        )


def _region_map(node: dict, context: str) -> dict:
    values = {}
    for region, value in node.items():
        if not isinstance(region, str):
            raise TranspileError(f"{context}: region name {region!r} must be a string.")
        values[region] = _require_number(value, f"{context}.{region}")
    return values


# ------------------------------------------------------------------ object specs

def _finalize_object(section_key: str, type_name: str, cls: type, node: dict, context: str,
                      tensor_aware: bool = False) -> ObjectSpec:
    """Validate a resolved ``(cls, kwargs)`` pair and build its :class:`ObjectSpec`."""
    node = signatures.alias_kwargs(node, cls)
    signatures.validate_kwargs(cls, node, context)

    tensor_names = signatures.constitutive_tensor_params(cls) if tensor_aware else set()
    kwargs, tensor_kwargs = {}, {}
    for key, value in node.items():
        key_context = f"{context}.{key}"
        if key in tensor_names:
            if isinstance(value, dict):
                tensor_kwargs[key] = _region_map(value, key_context)
            else:
                tensor_kwargs[key] = _require_number(value, key_context)
        else:
            if isinstance(value, dict):
                raise TranspileError(
                    f"{key_context}: mappings are only valid for per-element "
                    f"(tensor) parameters; '{key}' of {type_name} is not one."
                )
            kwargs[key] = _coerce_scalar(value, key_context)

    return ObjectSpec(section=section_key, type_name=type_name, cls=cls,
                      kwargs=kwargs, tensor_kwargs=tensor_kwargs)


def build_object(section_key: str, node, context: str, tensor_aware: bool = False) -> ObjectSpec:
    """Resolve and validate one ``{type: ..., <kwargs>}`` YAML block."""
    node = dict(_require_mapping(node, context))
    type_name = node.pop("type", None)
    if not isinstance(type_name, str):
        raise TranspileError(f"{context}: a 'type' key naming a class is required.")

    cls = registry.resolve(section_key, type_name, context)
    return _finalize_object(section_key, type_name, cls, node, context, tensor_aware=tensor_aware)


# --------------------------------------------------------------------- sections

def _parse_equations(node) -> dict:
    node = _require_mapping(node, "equations")
    _check_keys(node, ("momentum", "heat"), "equations")
    if not node:
        raise TranspileError("equations: define at least one of 'momentum', 'heat'.")
    return {
        name: build_object(f"equations.{name}", spec, f"equations.{name}")
        for name, spec in node.items()
    }


def _material_property_names() -> list:
    """Material scalar properties = the real ``set_*`` methods of Material."""
    from safeincave import Material
    return sorted(
        name[len("set_"):]
        for name in vars(Material)
        if name.startswith("set_") and callable(getattr(Material, name))
    )


_ELEMENT_GROUPS = ("elastic", "non_elastic", "thermoelastic")


def _merge_mapping_list(node, context: str) -> dict:
    """Merge a list of one-or-more-key mappings into a single flat dict.

    Also accepts a single mapping directly (already "merged"), so both
    ``elastic: [{Youngs_modulus: ...}, {Poissons_ratio: ...}]`` and
    ``elastic: {Youngs_modulus: ..., Poissons_ratio: ...}`` are legal.
    """
    if isinstance(node, dict):
        return dict(node)
    merged = {}
    for i, item in enumerate(_require_list(node, context)):
        item = _require_mapping(item, f"{context}[{i}]")
        for key, value in item.items():
            if key in merged:
                raise TranspileError(f"{context}[{i}]: duplicate key {key!r}.")
            merged[key] = value
    return merged


def _resolve_material_element(key: str, context: str):
    """Resolve a material key (e.g. ``plasticDPR`` or ``elastic``) to its class.

    Direct, case-insensitive type names (``plasticDPR`` -> ``PlasticDPR``) are
    tried first; a bare category name (``elastic``) is accepted as shorthand
    only when that category has exactly one legal type (currently true for
    'elastic' -> Spring and 'thermoelastic' -> Thermoelastic). ``key`` is
    first translated through ``signatures.TYPE_ALIASES`` (e.g.
    ``plastic_drucker_prager`` -> ``PlasticDPR``), if it matches one.
    """
    key = signatures.resolve_type_alias(key)
    matches = []
    for category in _ELEMENT_GROUPS:
        for legal in registry.legal_names(f"material.{category}"):
            if legal.lower() == key.lower():
                matches.append((category, legal))
    if len(matches) == 1:
        category, type_name = matches[0]
        return category, type_name, registry.resolve(f"material.{category}", type_name, context)
    if len(matches) > 1:
        options = ", ".join(f"{cat}.{name}" for cat, name in matches)
        raise TranspileError(f"{context}: '{key}' is ambiguous ({options}).")

    if key in _ELEMENT_GROUPS:
        legal = registry.legal_names(f"material.{key}")
        if len(legal) == 1:
            type_name = legal[0]
            return key, type_name, registry.resolve(f"material.{key}", type_name, context)
        raise TranspileError(
            f"{context}: '{key}' has {len(legal)} legal types ({', '.join(legal)}); "
            f"use the specific type name as the key instead, e.g. '{legal[0]}'."
        )

    all_legal = sorted(
        {name for cat in _ELEMENT_GROUPS for name in registry.legal_names(f"material.{cat}")}
        | set(_ELEMENT_GROUPS)
    )
    lines = [f"{context}: unknown material key '{key}'."]
    suggestions = closest(key, all_legal, n=3)
    if suggestions:
        lines.append(f"Did you mean: {', '.join(suggestions)}?")
    lines.append(f"Legal keys: {', '.join(all_legal)}.")
    raise TranspileError("\n".join(lines))


def _parse_material(node, equations: dict) -> MaterialSpec:
    node = _require_mapping(node, "material")
    properties_allowed = _material_property_names()
    spec = MaterialSpec()

    if "name" in node:
        name = node["name"]
        if not isinstance(name, str):
            raise TranspileError("material.name: expected a string.")
        spec.name = name

    for key, value in node.items():
        if key == "name":
            continue

        context = f"material.{key}"
        if key in properties_allowed:
            if isinstance(value, dict):
                spec.properties[key] = _region_map(value, context)
            else:
                spec.properties[key] = _require_number(value, context)
            continue

        category, type_name, cls = _resolve_material_element(key, context)
        raw_kwargs = _merge_mapping_list(value, context)
        elem = _finalize_object(
            f"material.{category}", type_name, cls, raw_kwargs, context, tensor_aware=True,
        )
        getattr(spec, category).append(elem)

    spec.properties.setdefault("density", 0.0)

    if "momentum" in equations and not spec.elastic:
        raise TranspileError(
            "material: at least one 'elastic' element is required for a momentum equation."
        )
    if "heat" in equations:
        for prop in ("specific_heat_capacity", "thermal_conductivity"):
            if prop not in spec.properties:
                raise TranspileError(
                    f"material: '{prop}' is required for a heat equation."
                )
    return spec


# ------------------------------------------------------- named-block schema
#
# An alternate top-level schema (selected by the presence of a 'steps' key):
# 'materials'/'boundaries'/'loads'/'outputs'/'logging' are named collections
# defined once, and 'steps' (an ordered mapping, replacing 'stages') composes
# them by name. A step that omits a field inherits it verbatim from the
# previous step, so only genuine deltas need to be restated.
#
# Kept in sync with the identical schema in the
# extensions/FullNewtonRaphsonSolver/Transpiler/model.py overlay, which adds
# 'initial_stress'-kind loads (a Newton-path-only feature core does not have
# a CaseModel field for) on top of this.

_DIRECTION_COMPONENTS = {"x": 0, "y": 1, "z": 2}


def _by_name(items, context: str) -> dict:
    """A list of ``{name: ..., ...}`` blocks, keyed by their 'name'."""
    parsed = {}
    for i, item in enumerate(_require_list(items, context)):
        item = _require_mapping(item, f"{context}[{i}]")
        name = item.get("name")
        if not isinstance(name, str):
            raise TranspileError(f"{context}[{i}]: a string 'name' is required.")
        if name in parsed:
            raise TranspileError(f"{context}: duplicate name {name!r}.")
        parsed[name] = item
    return parsed


def _select_names(value, context: str) -> list:
    """A step reference field: a bare name, or a list of names."""
    if isinstance(value, str):
        return [value]
    return _require_list(value, context)


def _lookup_all(names: list, collection: dict, context: str) -> list:
    missing = [name for name in names if name not in collection]
    if missing:
        raise TranspileError(
            f"{context}: unknown name(s) {', '.join(map(repr, missing))}. "
            f"Defined: {', '.join(sorted(collection)) or '(none)'}."
        )
    return [collection[name] for name in names]


def _boundary_names(value, context: str) -> list:
    """'sets:' on a boundary/load block: a single boundary name or a list.

    Unlike 'materials.sets', 'ALL' has no meaning here (there is no "every
    boundary" concept) so it is rejected with a clear message instead of
    silently doing nothing.
    """
    if value == "ALL":
        raise TranspileError(f"{context}: 'ALL' is not valid here; list specific boundary names.")
    if isinstance(value, str):
        return [value]
    return _require_list(value, context)


def _type_and_kwargs(value, context: str) -> tuple:
    """A bare type name, or a dict of {type: ..., **kwargs}."""
    if isinstance(value, str):
        return value, {}
    node = dict(_require_mapping(value, context))
    type_name = node.pop("type", None)
    if not isinstance(type_name, str):
        raise TranspileError(f"{context}: a 'type' key naming a class is required.")
    return type_name, node


def _parse_materials(node, equations: dict, context: str = "materials") -> dict:
    parsed = {}
    for name, entry in _by_name(node, context).items():
        entry_context = f"{context}.{name}"
        entry = dict(entry)
        entry.pop("name")
        sets_raw = entry.pop("sets", ["all_elements"])
        sets = [str(s) for s in _require_list(sets_raw, f"{entry_context}.sets")]
        models = _require_list(entry.pop("models", []), f"{entry_context}.models")
        _check_keys(entry, (), entry_context)

        spec = MaterialSpec(name=name, sets=sets)
        for i, model in enumerate(models):
            model_context = f"{entry_context}.models[{i}]"
            model = dict(_require_mapping(model, model_context))
            model_type = model.pop("type", None)
            if not isinstance(model_type, str):
                raise TranspileError(f"{model_context}: a 'type' key is required.")
            category, type_name, cls = _resolve_material_element(model_type, model_context)
            elem = _finalize_object(
                f"material.{category}", type_name, cls, model, model_context, tensor_aware=True,
            )
            getattr(spec, category).append(elem)

        if "momentum" in equations and not spec.elastic:
            raise TranspileError(
                f"{entry_context}: at least one elastic model is required for a "
                f"momentum equation."
            )
        spec.properties.setdefault("density", 0.0)
        parsed[name] = spec
    return parsed


def _merge_material_specs(named_specs: list, context: str) -> "MaterialSpec":
    """Combine several region-scoped :class:`MaterialSpec`\\ s into one.

    Each input material's per-parameter values are turned into a
    ``{region: value}`` tensor map (or kept a plain scalar, when every
    combined material happens to agree on that parameter's value) -- the
    same region-map shape ``_finalize_object``/codegen already understand
    for a single material's parameters, just built automatically here from
    several named, region-scoped materials instead of by hand.

    ``named_specs`` is a list of ``(name, sets, MaterialSpec)`` tuples. All
    of them must define the same elastic/non_elastic/thermoelastic model
    types, in the same order, and the same set of properties -- only the
    *values* may differ by region -- and no two may claim the same region
    (nor use the ``ALL`` sentinel, which cannot be combined with real
    region names).
    """
    first_name, _, first_spec = named_specs[0]
    seen_regions: dict = {}
    for name, sets, spec in named_specs:
        for group in _ELEMENT_GROUPS:
            types = [e.type_name for e in getattr(spec, group)]
            first_types = [e.type_name for e in getattr(first_spec, group)]
            if types != first_types:
                raise TranspileError(
                    f"{context}: materials {first_name!r} and {name!r} must define the "
                    f"same {group} model(s) in the same order to be combined by region "
                    f"(got {first_types} vs {types})."
                )
        if set(spec.properties) != set(first_spec.properties):
            raise TranspileError(
                f"{context}: materials {first_name!r} and {name!r} must define the same "
                f"properties (got {sorted(first_spec.properties)} vs {sorted(spec.properties)})."
            )
        for region in sets:
            if region == "ALL":
                raise TranspileError(
                    f"{context}: material {name!r} uses sets: ['ALL'], which cannot be "
                    f"combined with other region-scoped materials; give it explicit "
                    f"region names instead."
                )
            if region in seen_regions:
                raise TranspileError(
                    f"{context}: region {region!r} is claimed by both material "
                    f"{seen_regions[region]!r} and {name!r}."
                )
            seen_regions[region] = name

    def _region_value(values_by_region: dict):
        distinct = set(values_by_region.values())
        return next(iter(distinct)) if len(distinct) == 1 else dict(values_by_region)

    merged = MaterialSpec(name="+".join(name for name, _, _ in named_specs))
    for prop in first_spec.properties:
        values_by_region = {
            region: spec.properties[prop]
            for name, sets, spec in named_specs for region in sets
        }
        merged.properties[prop] = _region_value(values_by_region)

    for group in _ELEMENT_GROUPS:
        for idx in range(len(getattr(first_spec, group))):
            elems = [getattr(spec, group)[idx] for _, _, spec in named_specs]
            type_name = elems[0].type_name
            all_keys = {key for elem in elems for key in (*elem.kwargs, *elem.tensor_kwargs)}
            merged_kwargs, merged_tensor_kwargs = {}, {}
            for key in all_keys:
                # A key is a tensor (constitutive, per-element) parameter for
                # this class if _finalize_object put it in tensor_kwargs for
                # *any* elem -- that's a property of the class, not of the
                # individual call, so it is consistent across all of them.
                # Its value there may be a plain scalar (the common case,
                # mergeable below) or already a {region: value} map (not
                # mergeable further -- rejected below).
                is_tensor_key = any(key in elem.tensor_kwargs for elem in elems)
                values_by_region = {}
                for (name, sets, _), elem in zip(named_specs, elems):
                    bucket = elem.tensor_kwargs if is_tensor_key else elem.kwargs
                    if key not in bucket:
                        raise TranspileError(
                            f"{context}: {type_name}.{key} is given in some but not all "
                            f"combined materials (missing from {name!r}); give it "
                            f"explicitly everywhere being combined, even if to the same "
                            f"value."
                        )
                    value = bucket[key]
                    if isinstance(value, dict):
                        raise TranspileError(
                            f"{context}: {type_name}.{key} in material {name!r} is "
                            f"already a per-region map; combining region-scoped "
                            f"materials that also use per-parameter region maps is not "
                            f"supported."
                        )
                    for region in sets:
                        values_by_region[region] = value
                resolved = _region_value(values_by_region)
                if is_tensor_key:
                    merged_tensor_kwargs[key] = resolved
                elif isinstance(resolved, dict):
                    raise TranspileError(
                        f"{context}: {type_name}.{key} differs by region but is not a "
                        f"per-element (constitutive tensor) parameter; only those can "
                        f"vary by region."
                    )
                else:
                    merged_kwargs[key] = resolved
            merged_elem = ObjectSpec(
                section=elems[0].section, type_name=type_name, cls=elems[0].cls,
                kwargs=merged_kwargs, tensor_kwargs=merged_tensor_kwargs,
            )
            getattr(merged, group).append(merged_elem)

    return merged


def _parse_boundaries(node, context: str = "boundaries") -> dict:
    parsed = {}
    for name, entry in _by_name(node, context).items():
        entry_context = f"{context}.{name}"
        entry = dict(entry)
        entry.pop("name")
        _check_keys(entry, ("type", "value", "sets"), entry_context)
        if entry.get("type") != "displacement":
            raise TranspileError(
                f"{entry_context}.type: only 'displacement' is currently supported, "
                f"got {entry.get('type')!r}."
            )
        value_node = _require_mapping(entry.get("value"), f"{entry_context}.value")
        axes = [axis for axis in ("x", "y", "z") if axis in value_node]
        if len(axes) != 1:
            raise TranspileError(
                f"{entry_context}.value: expected exactly one of 'x'/'y'/'z', "
                f"got {list(value_node)}."
            )
        axis = axes[0]
        value = _require_number(value_node[axis], f"{entry_context}.value.{axis}")
        sets = _boundary_names(entry.get("sets"), f"{entry_context}.sets")
        parsed[name] = {"sets": sets, "component": _DIRECTION_COMPONENTS[axis], "value": value}
    return parsed


def _parse_loads(node, context: str = "loads") -> dict:
    parsed = {}
    for name, entry in _by_name(node, context).items():
        entry_context = f"{context}.{name}"
        entry = dict(entry)
        entry.pop("name")
        load_type = entry.pop("type", None)
        if load_type == "initial_stress":
            raise TranspileError(
                f"{entry_context}: 'initial_stress' loads require the Newton-path "
                f"extension (extensions/FullNewtonRaphsonSolver), which adds "
                f"'apply_initial_stress' support; core does not have it."
            )
        elif load_type == "pressure":
            _check_keys(entry, ("value", "time", "sets"), entry_context)
            sets = _boundary_names(entry.get("sets"), f"{entry_context}.sets")
            value = entry.get("value")
            if isinstance(value, list):
                pressure = [
                    _require_number(v, f"{entry_context}.value[{i}]")
                    for i, v in enumerate(value)
                ]
            else:
                scalar = _require_number(value, f"{entry_context}.value")
                pressure = [scalar, scalar]
            own_time = None
            if "time" in entry:
                time_raw = _require_list(entry["time"], f"{entry_context}.time")
                if len(time_raw) != 2:
                    raise TranspileError(f"{entry_context}.time: expected [start, end].")
                own_time = (
                    _require_number(time_raw[0], f"{entry_context}.time[0]"),
                    _require_number(time_raw[1], f"{entry_context}.time[1]"),
                )
            parsed[name] = {
                "kind": "pressure", "sets": sets, "pressure": pressure, "time": own_time,
            }
        else:
            raise TranspileError(f"{entry_context}.type: expected 'pressure', got {load_type!r}.")
    return parsed


def _parse_outputs_defs(node, context: str = "outputs") -> dict:
    parsed = {}
    for name, entry in _by_name(node, context).items():
        entry_context = f"{context}.{name}"
        entry = dict(entry)
        entry.pop("name")
        kind = entry.pop("type", None)
        if kind == "results":
            _check_keys(entry, ("fields", "merged_solutions", "smooth_output"), entry_context)
            fields_raw = entry.get("fields")
            if isinstance(fields_raw, dict):
                fields = {
                    signatures.resolve_field_alias(str(k)): str(v)
                    for k, v in fields_raw.items()
                }
            else:
                fields = {
                    signatures.resolve_field_alias(str(f)): str(f)
                    for f in _require_list(fields_raw, f"{entry_context}.fields")
                }
            if not fields:
                raise TranspileError(f"{entry_context}: 'fields' must not be empty.")
            parsed[name] = {
                "kind": "results",
                "fields": fields,
                "merged_solutions": bool(entry.get("merged_solutions", False)),
                "smooth_output": bool(entry.get("smooth_output", False)),
            }
        elif kind == "logging":
            _check_keys(entry, ("target_point", "variables_to_track"), entry_context)
            logging_node = dict(entry)
            logging_node["type"] = "SimulationLogging"
            parsed[name] = {
                "kind": "logging",
                "spec": build_object("logging", logging_node, entry_context),
            }
        else:
            raise TranspileError(
                f"{entry_context}.type: expected 'results' or 'logging', got {kind!r}."
            )
    return parsed


def _default_time_controller_kwargs(start: float, end: float) -> dict:
    span = end - start
    return {
        "initial_time": start,
        "final_time": end,
        "time_unit": "second",
        "initial_dt": span / 10.0,
        "dt_min": span / 1000.0,
        "dt_max": span / 2.0,
        "shrink_factor": 0.5,
        "growth_factor": 1.5,
        "easy_ratio_threshold": 0.25,
        "hard_ratio_threshold": 0.5,
        "max_bisections": 10,
    }


def _boundary_bcs(names: list, boundaries: dict, start: float, end: float, context: str) -> list:
    defs = _lookup_all(names, boundaries, context)
    cls = registry.resolve("bcs.momentum", "DirichletBC", context)
    specs = []
    for bdef in defs:
        for boundary_name in bdef["sets"]:
            kwargs = {
                "boundary_name": boundary_name,
                "component": bdef["component"],
                "values": [bdef["value"], bdef["value"]],
                "time_values": [start, end],
            }
            specs.append(_finalize_object("bcs.momentum", "DirichletBC", cls, kwargs, context))
    return specs


def _load_bcs(names: list, loads: dict, start: float, end: float, context: str) -> list:
    """Fan out 'pressure' loads to NeumannBCs.

    A pressure load's own 'time:' (if given) overrides the containing step's
    (start, end) for that load's NeumannBC time_values.
    """
    defs = _lookup_all(names, loads, context)
    cls = registry.resolve("bcs.momentum", "NeumannBC", context)
    specs = []
    for ldef in defs:
        load_start, load_end = ldef["time"] or (start, end)
        for boundary_name in ldef["sets"]:
            kwargs = {
                "boundary_name": boundary_name,
                "direction": 2,
                "values": ldef["pressure"],
                "time_values": [load_start, load_end],
            }
            specs.append(_finalize_object("bcs.momentum", "NeumannBC", cls, kwargs, context))
    return specs


def _parse_solvers(node, context: str = "solvers") -> dict:
    parsed = {}
    for name, entry in _by_name(node, context).items():
        entry_context = f"{context}.{name}"
        entry = dict(entry)
        entry.pop("name")

        solver_type = entry.pop("type", None)
        if not isinstance(solver_type, str):
            raise TranspileError(f"{entry_context}: a 'type' key is required.")

        if "time" not in entry:
            raise TranspileError(f"{entry_context}: 'time' is required.")
        time_raw = _require_list(entry.pop("time"), f"{entry_context}.time")
        if len(time_raw) != 2:
            raise TranspileError(f"{entry_context}.time: expected [start, end].")
        start = _require_number(time_raw[0], f"{entry_context}.time[0]")
        end = _require_number(time_raw[1], f"{entry_context}.time[1]")

        momentum_spec = None
        if "momentum" in entry:
            mom_context = f"{entry_context}.momentum"
            mom_type, mom_kwargs = _type_and_kwargs(entry.pop("momentum"), mom_context)
            mom_cls = registry.resolve("equations.momentum", mom_type, mom_context)
            momentum_spec = _finalize_object(
                "equations.momentum", mom_type, mom_cls, mom_kwargs, mom_context
            )

        time_type, time_overrides = _type_and_kwargs(
            entry.pop("time_step", "TimeControllerAdaptive"), f"{entry_context}.time_step"
        )
        time_kwargs = _default_time_controller_kwargs(start, end)
        time_kwargs.update(time_overrides)
        time_cls = registry.resolve("time", time_type, f"{entry_context}.time_step")
        time_spec = _finalize_object(
            "time", time_type, time_cls, time_kwargs, f"{entry_context}.time_step"
        )

        parsed[name] = {
            "type": solver_type, "kwargs": entry,
            "start": start, "end": end,
            "momentum": momentum_spec, "time_spec": time_spec,
        }
    return parsed


def _resolve_equations(solvers: dict, steps_node, context: str = "steps") -> dict:
    """Determine the momentum equation from the solvers actually referenced by steps.

    Equations are built once, before any step runs (the real API always
    constructs the momentum equation once and passes it by reference into
    every stage's simulator), so every solver referenced across all steps
    (inheriting a step's 'solver:' name forward when omitted) must agree.
    """
    steps_list = _require_list(steps_node, context)
    carried_solver = None
    referenced = []
    for i, step_node in enumerate(steps_list):
        step_node = _require_mapping(step_node, f"{context}[{i}]")
        if "solver" in step_node:
            carried_solver = step_node["solver"]
        if carried_solver is not None and carried_solver not in referenced:
            referenced.append(carried_solver)

    momentum_spec = None
    for solver_name in referenced:
        if solver_name not in solvers:
            raise TranspileError(
                f"{context}: unknown solver {solver_name!r}. "
                f"Defined: {', '.join(sorted(solvers)) or '(none)'}."
            )
        spec = solvers[solver_name]["momentum"]
        if spec is None:
            continue
        if momentum_spec is None:
            momentum_spec = spec
        elif (spec.type_name, spec.kwargs) != (momentum_spec.type_name, momentum_spec.kwargs):
            raise TranspileError(
                f"{context}: solvers referenced across steps define different 'momentum' "
                f"specs; the momentum equation is built once and must be identical."
            )
    return {"momentum": momentum_spec} if momentum_spec is not None else {}


def _parse_steps(node, equations: dict, materials: dict, boundaries: dict, loads: dict,
                 outputs_defs: dict, solvers: dict, context: str = "steps") -> tuple:
    steps_list = _require_list(node, context)
    if not steps_list:
        raise TranspileError(f"{context}: at least one step is required.")

    stages = []
    resolved_material_names = None
    material_spec = MaterialSpec()
    carried: dict = {}

    for i, step_node in enumerate(steps_list):
        step_node = _require_mapping(step_node, f"{context}[{i}]")
        step_name = step_node.get("name")
        if not isinstance(step_name, str):
            raise TranspileError(f"{context}[{i}]: a string 'name' is required.")
        step_context = f"{context}.{step_name}"
        _check_keys(
            step_node, ("name", "materials", "boundaries", "loads", "solver", "outputs"),
            step_context,
        )
        merged = {**carried, **{k: v for k, v in step_node.items() if k != "name"}}
        carried = merged

        if "solver" not in merged:
            raise TranspileError(f"{step_context}: 'solver' is required (directly or inherited).")
        solver_name = merged["solver"]
        if solver_name not in solvers:
            raise TranspileError(
                f"{step_context}.solver: unknown solver {solver_name!r}. "
                f"Defined: {', '.join(sorted(solvers)) or '(none)'}."
            )
        solver_def = solvers[solver_name]
        start, end, time_spec = solver_def["start"], solver_def["end"], solver_def["time_spec"]

        material_names = (
            _select_names(merged["materials"], f"{step_context}.materials")
            if "materials" in merged else []
        )
        if "momentum" in equations and not material_names:
            raise TranspileError(
                f"{step_context}: 'materials' is required (directly or inherited)."
            )
        if material_names:
            step_materials = _lookup_all(material_names, materials, f"{step_context}.materials")
            if len(step_materials) > 1:
                named_specs = list(zip(material_names, (m.sets for m in step_materials), step_materials))
                combined_material = _merge_material_specs(named_specs, f"{step_context}.materials")
            else:
                combined_material = step_materials[0]
            if resolved_material_names is None:
                resolved_material_names = material_names
                material_spec = combined_material
            elif resolved_material_names != material_names:
                raise TranspileError(
                    f"{step_context}.materials: changing the material between steps "
                    f"is not yet supported (material is built once, before any step runs)."
                )

        bcs_momentum = []
        if "momentum" in equations:
            boundary_names = (
                _select_names(merged["boundaries"], f"{step_context}.boundaries")
                if "boundaries" in merged else []
            )
            bcs_momentum += _boundary_bcs(
                boundary_names, boundaries, start, end, f"{step_context}.boundaries"
            )

            load_names = (
                _select_names(merged["loads"], f"{step_context}.loads")
                if "loads" in merged else []
            )
            bcs_momentum += _load_bcs(load_names, loads, start, end, f"{step_context}.loads")

        bcs = {"momentum": bcs_momentum} if "momentum" in equations else {}

        output_names = (
            _select_names(merged["outputs"], f"{step_context}.outputs")
            if "outputs" in merged else []
        )
        if not output_names:
            raise TranspileError(f"{step_context}: 'outputs' is required (directly or inherited).")
        output_defs = _lookup_all(output_names, outputs_defs, f"{step_context}.outputs")
        fields = {}
        merged_solutions_values = set()
        smooth_output_values = set()
        logging_specs = []
        for odef in output_defs:
            if odef["kind"] == "results":
                fields.update(odef["fields"])
                merged_solutions_values.add(odef["merged_solutions"])
                smooth_output_values.add(odef["smooth_output"])
            else:
                logging_specs.append(odef["spec"])
        if len(merged_solutions_values) > 1 or len(smooth_output_values) > 1:
            raise TranspileError(
                f"{step_context}.outputs: combined output profiles disagree on "
                f"'merged_solutions'/'smooth_output'."
            )
        if len(logging_specs) > 1:
            raise TranspileError(
                f"{step_context}.outputs: only one logging profile per step is supported."
            )
        if not fields:
            raise TranspileError(f"{step_context}.outputs: no 'results'-type output referenced.")
        logging = logging_specs[0] if logging_specs else None
        outputs = [OutputSpec(
            equation="momentum" if "momentum" in equations else next(iter(equations)),
            output_format="xdmf", folder=f"output/{step_name}", fields=fields,
        )]

        solver_kwargs = dict(solver_def["kwargs"])
        solver_kwargs["merged_solutions"] = merged_solutions_values.pop() if merged_solutions_values else False
        solver_kwargs["smooth_output"] = smooth_output_values.pop() if smooth_output_values else False
        solver_cls = registry.resolve("simulator", solver_def["type"], f"{step_context}.solver")
        simulator = _finalize_object(
            "simulator", solver_def["type"], solver_cls, solver_kwargs, f"{step_context}.solver"
        )
        _check_simulator_wiring(simulator, equations, [], logging, f"{step_context}.solver")

        stages.append(StageSpec(
            name=step_name, time=time_spec, bcs=bcs, caverns=[], logging=logging,
            outputs=outputs, simulator=simulator,
        ))

    return stages, material_spec


def _parse_output(node, equations: dict, context: str) -> OutputSpec:
    from safeincave.Output.SaveFields import WRITER_BACKENDS

    node = _require_mapping(node, context)
    _check_keys(node, ("equation", "output_format", "folder", "fields"), context)

    if "equation" in node:
        equation = node["equation"]
    elif len(equations) == 1:
        equation = next(iter(equations))
    else:
        raise TranspileError(
            f"{context}: 'equation' is required when several equations are defined."
        )
    if equation not in equations:
        raise TranspileError(
            f"{context}: equation '{equation}' is not defined under 'equations'. "
            f"Defined: {', '.join(equations)}."
        )

    output_format = node.get("output_format", "xdmf")
    if output_format not in WRITER_BACKENDS:
        raise TranspileError(
            f"{context}: unknown output_format '{output_format}'. "
            f"Available: {', '.join(sorted(WRITER_BACKENDS))}."
        )

    folder = node.get("folder")
    if not isinstance(folder, str):
        raise TranspileError(f"{context}: 'folder' (string) is required.")

    fields = _require_mapping(node.get("fields", {}), f"{context}.fields")
    if not fields:
        raise TranspileError(f"{context}: 'fields' must map field names to labels.")
    fields = {str(name): str(label) for name, label in fields.items()}

    return OutputSpec(equation=equation, output_format=output_format,
                      folder=folder, fields=fields)


def _parse_stage(node, index: int, equations: dict) -> StageSpec:
    context = f"stages[{index}]"
    node = _require_mapping(node, context)
    _check_keys(
        node,
        ("name", "time", "bcs", "caverns", "logging", "outputs", "simulator"),
        context,
    )
    name = str(node.get("name", f"stage_{index}"))

    if "time" not in node:
        raise TranspileError(f"{context}: a 'time' block is required.")
    time = build_object("time", node["time"], f"{context}.time")

    bcs_node = _require_mapping(node.get("bcs", {}), f"{context}.bcs")
    _check_keys(bcs_node, ("momentum", "heat"), f"{context}.bcs")
    bcs = {}
    for eq_name in bcs_node:
        if eq_name not in equations:
            raise TranspileError(
                f"{context}.bcs: equation '{eq_name}' is not defined under 'equations'."
            )
        bcs[eq_name] = [
            build_object(f"bcs.{eq_name}", bc_node, f"{context}.bcs.{eq_name}[{i}]")
            for i, bc_node in enumerate(
                _require_list(bcs_node[eq_name], f"{context}.bcs.{eq_name}")
            )
        ]
    for eq_name in equations:
        if eq_name not in bcs:
            raise TranspileError(
                f"{context}.bcs: boundary conditions for '{eq_name}' are required "
                f"(use an empty list to apply none)."
            )

    caverns = [
        build_object("caverns", cavern_node, f"{context}.caverns[{i}]")
        for i, cavern_node in enumerate(
            _require_list(node.get("caverns", []), f"{context}.caverns")
        )
    ]

    logging = None
    if "logging" in node:
        logging = build_object("logging", node["logging"], f"{context}.logging")

    outputs = [
        _parse_output(out_node, equations, f"{context}.outputs[{i}]")
        for i, out_node in enumerate(
            _require_list(node.get("outputs", []), f"{context}.outputs")
        )
    ]
    if not outputs:
        raise TranspileError(f"{context}: at least one entry under 'outputs' is required.")

    if "simulator" not in node:
        raise TranspileError(f"{context}: a 'simulator' block is required.")
    simulator = build_object("simulator", node["simulator"], f"{context}.simulator")
    _check_simulator_wiring(simulator, equations, caverns, logging, f"{context}.simulator")

    return StageSpec(name=name, time=time, bcs=bcs, caverns=caverns,
                     logging=logging, outputs=outputs, simulator=simulator)


def _check_simulator_wiring(simulator: ObjectSpec, equations: dict,
                            caverns: list, logging, context: str) -> None:
    """The simulator's auto-wired parameters define what the stage needs."""
    params = signatures.init_parameters(simulator.cls)
    for param_name, eq_name in (("eq_mom", "momentum"), ("eq_heat", "heat")):
        param = params.get(param_name)
        if param is not None and param.default is param.empty and eq_name not in equations:
            raise TranspileError(
                f"{context}: {simulator.type_name} requires a '{eq_name}' equation "
                f"(parameter '{param_name}'), but none is defined under 'equations'."
            )
    if caverns and "caverns" not in params:
        raise TranspileError(
            f"{context}: {simulator.type_name} does not accept caverns, "
            f"but the stage defines some."
        )
    if logging is not None and "simulation_logger" not in params:
        raise TranspileError(
            f"{context}: {simulator.type_name} does not accept a simulation logger, "
            f"but the stage defines one."
        )


# ------------------------------------------------------------------------ parse

_INITIAL_TEMPERATURE_HELP = (
    "initial_temperature: {problem}. It must be a uniform value in Kelvin; "
    "spatially varying initial fields cannot be expressed in YAML. Use "
    "'sic y2p' and set the field in the generated Python script instead."
)

_TOP_LEVEL_KEYS = (
    "grid",
    "mesh",
    "equations",
    "material",
    "body_force",
    "initial_temperature",
    "stages",
    "materials",
    "boundaries",
    "loads",
    "solvers",
    "outputs",
    "steps",
)


def _parse_grid(root: dict, context: str) -> ObjectSpec:
    """Build the grid ObjectSpec from either 'grid' (full form) or 'mesh' (shorthand).

    'mesh: path/to/name.msh' (or 'mesh: {file: path/to/name.msh}') is sugar
    for a GridHandlerGMSH block, splitting the path into grid_folder and
    geometry_name.
    """
    if "grid" in root and "mesh" in root:
        raise TranspileError(f"{context}: specify either 'grid' or 'mesh', not both.")
    if "mesh" in root:
        mesh_value = root["mesh"]
        if isinstance(mesh_value, dict):
            _check_keys(mesh_value, ("file",), f"{context}.mesh")
            mesh_path = mesh_value.get("file")
        else:
            mesh_path = mesh_value
        if not isinstance(mesh_path, str):
            raise TranspileError(
                f"{context}.mesh: expected a string path (or {{file: path}}), got {mesh_path!r}."
            )
        grid_folder, filename = _posixpath_split(mesh_path)
        geometry_name, ext = _posixpath_splitext(filename)
        if ext.lower() != ".msh":
            raise TranspileError(f"{context}.mesh: expected a '.msh' file, got {mesh_path!r}.")
        cls = registry.resolve("grid", "GridHandlerGMSH", f"{context}.mesh")
        return ObjectSpec(
            section="grid", type_name="GridHandlerGMSH", cls=cls,
            kwargs={"geometry_name": geometry_name, "grid_folder": grid_folder or "."},
        )
    return build_object("grid", root["grid"], "grid")


def _posixpath_split(path: str) -> tuple[str, str]:
    head, _, tail = path.rpartition("/")
    return head, tail


def _posixpath_splitext(name: str) -> tuple[str, str]:
    stem, dot, ext = name.rpartition(".")
    if not dot:
        return name, ""
    return stem, f".{ext}"


def parse(yaml_path) -> CaseModel:
    """Load and fully validate a YAML case file."""
    yaml_path = Path(yaml_path)
    root, sources = include.load(yaml_path)

    root = _require_mapping(root, str(yaml_path))
    _check_keys(root, _TOP_LEVEL_KEYS, str(yaml_path))
    if "grid" not in root and "mesh" not in root:
        raise TranspileError(
            f"{yaml_path}: top-level section 'grid' (or its 'mesh' shorthand) is required."
        )

    grid = _parse_grid(root, str(yaml_path))

    uses_new_schema = "steps" in root
    if uses_new_schema:
        if "equations" in root:
            raise TranspileError(
                f"{yaml_path}: top-level 'equations' is not used with the 'steps' schema; "
                f"define a 'solvers' entry's 'momentum:' field instead."
            )
        solvers = _parse_solvers(root.get("solvers", []))
        equations = _resolve_equations(solvers, root["steps"])
        if not equations:
            raise TranspileError(
                f"{yaml_path}: at least one step must define 'momentum' "
                f"(directly or inherited from an earlier step)."
            )
    else:
        if "equations" not in root:
            raise TranspileError(f"{yaml_path}: top-level section 'equations' is required.")
        equations = _parse_equations(root["equations"])

    body_force = root.get("body_force", [0.0, 0.0, 0.0])
    body_force = _require_list(body_force, "body_force")
    if len(body_force) != 3:
        raise TranspileError("body_force: expected exactly three components [gx, gy, gz].")
    body_force = [_require_number(v, f"body_force[{i}]") for i, v in enumerate(body_force)]

    try:
        initial_temperature = _require_number(
            root.get("initial_temperature", 293.0), "initial_temperature"
        )
    except TranspileError:
        raise TranspileError(
            _INITIAL_TEMPERATURE_HELP.format(
                problem=f"expected a number, got {root['initial_temperature']!r}"
            )
        ) from None

    if uses_new_schema:
        if "stages" in root or "material" in root:
            raise TranspileError(
                f"{yaml_path}: specify either the legacy 'material'/'stages' schema or the "
                f"'materials'/'boundaries'/'loads'/'solvers'/'outputs'/'steps' schema, not a "
                f"mix of both."
            )
        materials = _parse_materials(root.get("materials", []), equations)
        boundaries = _parse_boundaries(root.get("boundaries", []))
        loads = _parse_loads(root.get("loads", []))
        outputs_defs = _parse_outputs_defs(root.get("outputs", []))
        stages, material = _parse_steps(
            root["steps"], equations, materials, boundaries, loads, outputs_defs, solvers
        )
    else:
        if "material" not in root:
            raise TranspileError(f"{yaml_path}: top-level section 'material' is required.")
        if "stages" not in root:
            raise TranspileError(f"{yaml_path}: top-level section 'stages' is required.")

        material = _parse_material(root["material"], equations)

        stages_node = _require_list(root["stages"], "stages")
        if not stages_node:
            raise TranspileError("stages: at least one stage is required.")
        stages = [_parse_stage(node, i, equations) for i, node in enumerate(stages_node)]

    return CaseModel(
        source=yaml_path,
        sources=sources,
        grid=grid,
        equations=equations,
        material=material,
        body_force=body_force,
        initial_temperature=initial_temperature,
        stages=stages,
    )
