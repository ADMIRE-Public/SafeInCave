# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Shared utilities for cavern boundary condition updates.

This module provides common functions for updating boundary conditions
based on cavern state. The update_cavern_bcs logic is shared between
momentum and heat equation BC handlers.
"""

import numpy as np


def compute_cavern_pressure_load(cavern, x_coord):
    """
    Compute the pressure load for a cavern boundary condition.

    This function calculates the load applied by a cavern on the boundary,
    including hydrostatic effects. Used by both momentum and heat BC handlers.

    Parameters
    ----------
    cavern : Cavern
        Cavern object with attributes:
        - direction : int
            Coordinate index for hydrostatic variation (0=x, 1=y, 2=z).
        - density : float
            Fluid/solid density ρ.
        - ref_pos : float
            Reference elevation H.
        - P : float
            Current gauge pressure.
        - gravity : float
            Gravitational acceleration.
    x_coord : ufl.core.expr.Expr
        Spatial coordinate vector for hydrostatic calculations.

    Returns
    -------
    ufl.core.expr.Expr
        Load expression to be applied on the boundary.

    Notes
    -----
    The applied load is modeled as:
    ``p(t) + ρ g (H - x[i])``
    where p(t) = -P (negative because pressure acts inward).
    """
    i = cavern.direction
    rho = cavern.density
    H = cavern.ref_pos
    p = -cavern.P
    load = p + rho * cavern.gravity * (H - x_coord[i])
    return load


def compute_cavern_heat_flux(cavern):
    """
    Compute the heat flux parameters for a cavern boundary condition.

    Parameters
    ----------
    cavern : Cavern
        Cavern object with attributes:
        - T : float
            Current cavern temperature.
        - h_conv : float
            Convective heat-transfer coefficient.

    Returns
    -------
    tuple
        (T_inf, h) where:
        - T_inf : float
            Ambient temperature at the boundary.
        - h : float
            Heat transfer coefficient.
    """
    T_inf = cavern.T
    h = cavern.h_conv
    return T_inf, h
