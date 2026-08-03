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
    properties: dict = field(default_factory=dict)  # {setter_suffix: scalar|region map}
    elastic: list = field(default_factory=list)
    non_elastic: list = field(default_factory=list)
    thermoelastic: list = field(default_factory=list)


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

def build_object(section_key: str, node, context: str, tensor_aware: bool = False) -> ObjectSpec:
    """Resolve and validate one ``{type: ..., <kwargs>}`` YAML block."""
    node = dict(_require_mapping(node, context))
    type_name = node.pop("type", None)
    if not isinstance(type_name, str):
        raise TranspileError(f"{context}: a 'type' key naming a class is required.")

    cls = registry.resolve(section_key, type_name, context)
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


def _parse_material(node, equations: dict) -> MaterialSpec:
    node = _require_mapping(node, "material")
    properties_allowed = _material_property_names()
    _check_keys(node, list(_ELEMENT_GROUPS) + properties_allowed, "material")

    spec = MaterialSpec()
    for prop in properties_allowed:
        if prop in node:
            value = node[prop]
            context = f"material.{prop}"
            if isinstance(value, dict):
                spec.properties[prop] = _region_map(value, context)
            else:
                spec.properties[prop] = _require_number(value, context)

    for group in _ELEMENT_GROUPS:
        for i, elem_node in enumerate(_require_list(node.get(group, []), f"material.{group}")):
            elem = build_object(
                f"material.{group}", elem_node, f"material.{group}[{i}]", tensor_aware=True
            )
            getattr(spec, group).append(elem)

    if "density" not in spec.properties:
        raise TranspileError("material: 'density' is required.")
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
    "equations",
    "material",
    "body_force",
    "initial_temperature",
    "stages",
)


def parse(yaml_path) -> CaseModel:
    """Load and fully validate a YAML case file."""
    yaml_path = Path(yaml_path)
    root, sources = include.load(yaml_path)

    root = _require_mapping(root, str(yaml_path))
    _check_keys(root, _TOP_LEVEL_KEYS, str(yaml_path))
    for key in ("grid", "equations", "material", "stages"):
        if key not in root:
            raise TranspileError(f"{yaml_path}: top-level section '{key}' is required.")

    grid = build_object("grid", root["grid"], "grid")
    equations = _parse_equations(root["equations"])
    material = _parse_material(root["material"], equations)

    body_force = root.get("body_force", [0.0, 0.0, 0.0])
    body_force = _require_list(body_force, "body_force")
    if len(body_force) != 3:
        raise TranspileError("body_force: expected exactly three components [gx, gy, gz].")
    body_force = [_require_number(v, f"body_force[{i}]") for i, v in enumerate(body_force)]

    if "initial_temperature" not in root:
        raise TranspileError(_INITIAL_TEMPERATURE_HELP.format(problem="is required"))
    try:
        initial_temperature = _require_number(
            root["initial_temperature"], "initial_temperature"
        )
    except TranspileError:
        raise TranspileError(
            _INITIAL_TEMPERATURE_HELP.format(
                problem=f"expected a number, got {root['initial_temperature']!r}"
            )
        ) from None

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
