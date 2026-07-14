# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from ..base import GeneralBC


class RobinBC(GeneralBC):
    """
    Time-dependent Robin (convective) boundary condition.

    The Robin condition typically has the form
    :math:`h (T - T_\\infty)` on the boundary, where ``h`` is a heat
    transfer coefficient and ``T_∞`` may be time-dependent.

    Parameters
    ----------
    boundary_name : str
        Named boundary in the mesh tags.
    values : list of float
        Ambient values (e.g., ``T_∞``) sampled over time.
    h : float
        Robin/convective coefficient.
    time_values : list of float
        Times corresponding to ``values``.

    Attributes
    ----------
    type : str
        Always ``"robin"``.
    h : float
        Robin coefficient.
    """

    def __init__(self, boundary_name: str, values: list, h: float, time_values: list):
        super().__init__(boundary_name, values, time_values)
        self.type = "robin"
        self.h = h
