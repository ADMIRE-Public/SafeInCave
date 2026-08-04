# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Keyword validation against real constructor signatures.

Every YAML block is checked against ``inspect.signature`` of the class it
instantiates: missing required parameters and unknown keywords are reported
with the full real signature, Abaqus-keyword style. There is no schema to
maintain — the code *is* the schema.

The small exceptions are ``KEYWORD_ALIASES``/``TYPE_ALIASES``/``FIELD_ALIASES``
below: a handful of friendlier YAML spellings for real constructor
parameter/class/output-field names, translated to the real name before
validation (or, for output fields, before use) runs. Everything above still
applies to the *translated* name -- these are not a parallel schema, just
spelling sugar for a few names that read poorly verbatim in YAML.
"""

from __future__ import annotations

import inspect

from .errors import TranspileError
from .suggest import closest

# Friendly YAML keyword -> real constructor parameter name.
KEYWORD_ALIASES: dict = {
    "youngs_modulus": "E",
    "poissons_ratio": "nu",
}

# Friendly YAML type name -> real class name.
TYPE_ALIASES: dict = {
    "plastic_drucker_prager": "PlasticDPR",
}

# Friendly YAML output field name -> real field name (e.g. on the momentum
# equation / SaveFields), used as the default label when so aliased.
FIELD_ALIASES: dict = {
    "displacements": "u",
}


def alias_kwargs(node: dict, cls) -> dict:
    """Translate ``KEYWORD_ALIASES`` keys in ``node`` to their real names.

    Only translates a key when ``cls``'s real signature actually uses the
    alias's target name and *not* the alias itself (e.g. ``Spring`` takes
    ``E``, so ``youngs_modulus`` -> ``E``; ``PlasticDPR`` already takes
    ``youngs_modulus`` natively, so it is left alone there) -- otherwise a
    class that happens to use the long name as its real parameter would get
    it silently renamed to something it doesn't accept.
    """
    real_params = init_parameters(cls)
    result = {}
    for key, value in node.items():
        target = KEYWORD_ALIASES.get(key)
        if target is not None and key not in real_params and target in real_params:
            result[target] = value
        else:
            result[key] = value
    return result


def resolve_field_alias(name: str) -> str:
    """Translate ``name`` through ``FIELD_ALIASES`` if it is a known alias."""
    return FIELD_ALIASES.get(name, name)


def resolve_type_alias(name: str) -> str:
    """Translate ``name`` through ``TYPE_ALIASES`` if it is a known alias."""
    return TYPE_ALIASES.get(name, name)


# Parameters supplied by the transpiler itself (wired to objects built from
# other YAML sections). They are not legal YAML keywords.
AUTO_WIRED = {
    "grid",
    "eq",
    "eq_mom",
    "eq_heat",
    "t_control",
    "outputs",
    "caverns",
    "simulation_logger",
}


def init_parameters(cls) -> dict:
    """Ordered ``{name: Parameter}`` of ``cls.__init__`` without ``self``."""
    sig = inspect.signature(cls.__init__)
    params = dict(sig.parameters)
    params.pop("self", None)
    return params


def format_signature(cls) -> str:
    """Human-readable constructor signature, e.g. ``Spring(E, nu, name='spring')``."""
    parts = []
    for name, param in init_parameters(cls).items():
        if param.kind in (param.VAR_POSITIONAL, param.VAR_KEYWORD):
            continue
        if param.default is param.empty:
            parts.append(name)
        else:
            parts.append(f"{name}={param.default!r}")
    return f"{cls.__name__}({', '.join(parts)})"


def _has_var_keyword(cls) -> bool:
    return any(
        p.kind is p.VAR_KEYWORD for p in init_parameters(cls).values()
    )


def validate_kwargs(cls, provided: dict, context: str) -> None:
    """Check the YAML keys of one block against the real ``__init__`` signature."""
    params = init_parameters(cls)
    explicit = {
        name: p
        for name, p in params.items()
        if p.kind in (p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY)
    }
    settable = {name for name in explicit if name not in AUTO_WIRED}

    missing = [
        name
        for name, p in explicit.items()
        if p.default is p.empty and name not in AUTO_WIRED and name not in provided
    ]
    unknown = [] if _has_var_keyword(cls) else [
        key for key in provided if key not in settable
    ]

    if not missing and not unknown:
        return

    lines = [f"{context}: invalid parameters for {cls.__name__}."]
    if missing:
        lines.append(f"Missing required: {', '.join(missing)}.")
    for key in unknown:
        suggestion = closest(key, sorted(settable), n=1)
        hint = f" Did you mean '{suggestion[0]}'?" if suggestion else ""
        if key in AUTO_WIRED:
            hint = " This parameter is wired automatically and must not appear in YAML."
        lines.append(f"Unknown parameter: '{key}'.{hint}")
    lines.append(f"Signature: {format_signature(cls)}")
    raise TranspileError("\n".join(lines))


def tensor_params(cls) -> set:
    """Constructor parameters annotated as torch tensors.

    Values for these parameters are per-element fields in safeincave; the
    generated code broadcasts YAML scalars with ``to.ones(grid.n_elems)`` and
    expands ``{region: value}`` maps via ``grid.region_indices``.
    """
    names = set()
    for name, param in init_parameters(cls).items():
        if "Tensor" in str(param.annotation):
            names.add(name)
    return names


def constitutive_tensor_params(cls) -> set:
    """Tensor parameters of a constitutive element.

    Falls back to "every parameter except ``name``" when the class carries no
    tensor annotations (e.g. extension classes without type hints), since all
    constitutive element inputs are per-element tensors by contract.
    """
    names = tensor_params(cls)
    if names:
        return names
    return {
        name
        for name, param in init_parameters(cls).items()
        if name != "name"
        and param.kind in (param.POSITIONAL_OR_KEYWORD, param.KEYWORD_ONLY)
    }
