# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Simulators for coupled and single-equation problems."""

from .base import Simulator  # noqa: F401
from .thermo_mechanical import Simulator_TM  # noqa: F401
from .mechanical import Simulator_M  # noqa: F401
from .thermal import Simulator_T  # noqa: F401
from .mechanical_out import Simulator_Mout  # noqa: F401
from .mechanical_newton import Simulator_MNewton  # noqa: F401
from .thermo_mechanical_newton import Simulator_TMNewton  # noqa: F401

__all__ = [
    "Simulator",
    "Simulator_TM",
    "Simulator_M",
    "Simulator_T",
    "Simulator_Mout",
    "Simulator_MNewton",
    "Simulator_TMNewton",
]
