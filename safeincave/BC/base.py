# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Shared base classes for boundary conditions.

This module provides common base classes that are shared between momentum
and heat boundary condition hierarchies. Both MomentumBC and HeatBC import
from here and extend as needed for equation-specific behavior.
"""

from abc import ABC
import numpy as np


class GeneralBC(ABC):
    """
    Base container for time-dependent boundary-condition data.

    This is the shared base class for all boundary condition types.
    Subclasses set the concrete fields and ``type`` identifier.

    Parameters
    ----------
    boundary_name : str
        Name/label of the boundary as defined by mesh tags.
    values : list of float
        Boundary values sampled at the times in ``time_values``.
    time_values : list of float
        Monotonically increasing times associated with ``values``.

    Attributes
    ----------
    boundary_name : str
        Name/label of the boundary as defined by mesh tags.
    type : str or None
        Boundary condition type identifier (e.g., ``"dirichlet"``, ``"neumann"``).
        Set by subclasses.
    values : list of float
        Boundary values sampled at the times in ``time_values``.
    time_values : list of float
        Monotonically increasing times associated with ``values``.
    """

    def __init__(self, boundary_name: str, values: list, time_values: list):
        self.boundary_name = boundary_name
        self.type = None
        self.values = values
        self.time_values = time_values


class DirichletBC(GeneralBC):
    """
    Time-dependent Dirichlet (essential) boundary condition.

    Parameters
    ----------
    boundary_name : str
        Named boundary in the mesh tags.
    values : list of float
        Prescribed values over time (interpolated with :func:`numpy.interp`).
    time_values : list of float
        Times corresponding to ``values`` (must be monotone).

    Attributes
    ----------
    type : str
        Always ``"dirichlet"``.
    """

    def __init__(self, boundary_name: str, values: list, time_values: list):
        super().__init__(boundary_name, values, time_values)
        self.type = "dirichlet"


class NeumannBC(GeneralBC):
    """
    Time-dependent Neumann (traction/flux) boundary condition.

    Parameters
    ----------
    boundary_name : str
        Named boundary in the mesh tags.
    values : list of float
        Flux/intensity values over time.
    time_values : list of float
        Times corresponding to ``values`` (must be monotone).

    Attributes
    ----------
    type : str
        Always ``"neumann"``.
    """

    def __init__(self, boundary_name: str, values: list, time_values: list):
        super().__init__(boundary_name, values, time_values)
        self.type = "neumann"
