# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Reflection-driven type-name resolution for the YAML transpiler.

Every ``type:`` value in a YAML case file is resolved against the *live*
safeincave package at transpile time — never against a hardcoded table — so
classes added by the extensions overlay (:mod:`safeincave.extensions`) or by
the auto-discovery in :mod:`safeincave.Materials.Constitutive` are legal YAML
keywords automatically, with their real constructor parameter names.

Each YAML section resolves names in a specific namespace and only accepts
classes of the matching category (checked via ``issubclass``), so e.g. a time
controller name cannot be used where a boundary condition is expected.
"""

from __future__ import annotations

import importlib
import inspect
from dataclasses import dataclass

from .errors import TranspileError
from .suggest import closest


@dataclass(frozen=True)
class Section:
    """Namespace + category rules for one kind of YAML ``type:`` entry."""

    key: str            # section identifier used in error messages
    module: str         # dotted module in which type names are resolved
    prefix: str         # spelling of that namespace in the generated script
    bases: tuple = ()   # dotted paths of accepted base classes (any match)
    name_prefixes: tuple = ()  # accepted class-name prefixes when no base fits
    exclude: tuple = () # class names that resolve but are not concrete choices
    aliases: dict = None  # friendlier YAML spelling -> real class name


SECTIONS = {
    "grid": Section(
        key="grid",
        module="safeincave",
        prefix="sf.",
        name_prefixes=("GridHandler",),
    ),
    "equations.momentum": Section(
        key="equations.momentum",
        module="safeincave",
        prefix="sf.",
        bases=("safeincave.LinearMomentumBase",),
        exclude=("LinearMomentumBase",),
        aliases={"momentum_newton": "LinearMomentumNewton"},
    ),
    "equations.heat": Section(
        key="equations.heat",
        module="safeincave",
        prefix="sf.",
        bases=("safeincave.HeatDiffusion",),
    ),
    "material.elastic": Section(
        key="material.elastic",
        module="safeincave.Materials.Constitutive",
        prefix="sf.",
        bases=("safeincave.Materials.Constitutive.Spring.Spring",),
    ),
    "material.non_elastic": Section(
        key="material.non_elastic",
        module="safeincave.Materials.Constitutive",
        prefix="sf.",
        bases=(
            "safeincave.Materials.Constitutive.NonElasticElement.NonElasticElement",
        ),
        exclude=("NonElasticElement",),
    ),
    "material.thermoelastic": Section(
        key="material.thermoelastic",
        module="safeincave.Materials.Constitutive",
        prefix="sf.",
        bases=("safeincave.Materials.Constitutive.Thermoelastic.Thermoelastic",),
    ),
    "time": Section(
        key="time",
        module="safeincave",
        prefix="sf.",
        bases=("safeincave.TimeControllerBase",),
        exclude=("TimeControllerBase",),
        aliases={"adaptive": "TimeControllerAdaptive"},
    ),
    "bcs.momentum": Section(
        key="bcs.momentum",
        module="safeincave.BC.Momentum",
        prefix="momBC.",
        bases=("safeincave.BC.base.GeneralBC",),
        exclude=("GeneralBC",),
    ),
    "bcs.heat": Section(
        key="bcs.heat",
        module="safeincave.BC.Heat",
        prefix="heatBC.",
        bases=("safeincave.BC.base.GeneralBC",),
        exclude=("GeneralBC",),
    ),
    "caverns": Section(
        key="caverns",
        module="safeincave.Cavern",
        prefix="sf.Cavern.",
        bases=("safeincave.Cavern.Cavern",),
        exclude=("Cavern",),
    ),
    "logging": Section(
        key="logging",
        module="safeincave",
        prefix="sf.",
        bases=("safeincave.SimulationLogging",),
    ),
    "simulator": Section(
        key="simulator",
        module="safeincave",
        prefix="sf.",
        bases=("safeincave.Simulation.Simulators.base.Simulator",),
        aliases={
            "newton": "Simulator_MNewton",
            "thermo_mechanical_newton": "Simulator_TMNewton",
            "geostatic": "GeostaticStep",
        },
    ),
}


def _import_attr(dotted_path: str):
    module_path, _, attr = dotted_path.rpartition(".")
    return getattr(importlib.import_module(module_path), attr)


def _accepts(section: Section, name: str, obj) -> bool:
    if not inspect.isclass(obj) or inspect.isabstract(obj):
        return False
    if name in section.exclude:
        return False
    if section.bases:
        bases = tuple(_import_attr(path) for path in section.bases)
        return issubclass(obj, bases)
    if section.name_prefixes:
        return name.startswith(section.name_prefixes)
    return True


def _public_names(module) -> list:
    names = getattr(module, "__all__", None)
    if names is None:
        names = [name for name in dir(module) if not name.startswith("_")]
    return names


def legal_names(section_key: str) -> list:
    """All class names currently legal as ``type:`` in the given section."""
    section = SECTIONS[section_key]
    module = importlib.import_module(section.module)
    names = []
    for name in _public_names(module):
        obj = getattr(module, name, None)
        if _accepts(section, name, obj):
            names.append(name)
    if section.aliases:
        names.extend(section.aliases)
    return sorted(set(names))


def resolve(section_key: str, type_name: str, context: str):
    """Resolve ``type_name`` in the section's namespace, or raise.

    ``type_name`` may be the real class name or one of the section's
    friendlier YAML aliases (see :attr:`Section.aliases`); aliases are
    resolved to their real class name before the usual lookup.

    Raises :class:`TranspileError` with the list of legal names and a
    closest-match suggestion when the name is unknown or of the wrong
    category.
    """
    section = SECTIONS[section_key]
    real_name = (section.aliases or {}).get(type_name, type_name)
    module = importlib.import_module(section.module)
    obj = getattr(module, real_name, None)
    if obj is not None and _accepts(section, real_name, obj):
        return obj

    legal = legal_names(section_key)
    lines = [f"{context}: unknown type '{type_name}'."]
    if obj is not None:
        lines[0] = (
            f"{context}: '{type_name}' exists but is not a valid choice here "
            f"(section '{section.key}')."
        )
    suggestions = closest(type_name, legal, n=3)
    if suggestions:
        lines.append(f"Did you mean: {', '.join(suggestions)}?")
    lines.append(f"Legal types for '{section.key}': {', '.join(legal)}.")
    raise TranspileError("\n".join(lines))
