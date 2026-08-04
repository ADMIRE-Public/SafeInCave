# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Generic yield tracking for multi-mechanism plasticity models.

Registers a configurable ``Yield`` variable in
:mod:`safeincave.Output.DataExtract`'s variable registry: indicates whether
any active yield surface (Drucker-Prager, Mohr-Coulomb, Rankine, etc.) is
active, regardless of which one.

Usage:
    Import this module to activate generic yield tracking:
    >>> from safeincave.Output import register_generic_yield
    >>> register_generic_yield()

    Then log ``"Yield"`` via :class:`safeincave.Output.SimulationLogging`.
"""

from __future__ import annotations

from safeincave.Output import register_variable

# Global registry of known yield surface variable names
_REGISTERED_YIELD_SURFACES = {"yield_indicator", "yield_indicator_h"}


def register_yield_surface(name: str) -> None:
    """Register a custom yield surface name for auto-detection."""
    _REGISTERED_YIELD_SURFACES.add(name.lower())


def _extract_generic_yield(elem_id: int, **kwargs) -> int:
    """Auto-detect yield from any registered surface."""
    yielding = 0
    # Check all registered yield surfaces
    for surface_name in _REGISTERED_YIELD_SURFACES:
        indicator = kwargs.get(surface_name)
        if indicator is not None and int(indicator[elem_id].item()):
            yielding = 1
            break
    return yielding


def register_generic_yield() -> None:
    """Register the generic 'Yield' variable in SimulationLogging."""
    register_variable("Yield", "Yield", _extract_generic_yield)


# Auto-register on import for backward compatibility
register_generic_yield()
