# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
import torch as to


class Thermoelastic:
    """
    Thermoelastic contribution producing thermal strain :math:`\\varepsilon_{th}
    = \\alpha\\,\\Delta T\\,I`.

    Parameters
    ----------
    alpha : torch.Tensor
        Linear coefficient of thermal expansion per element, shape (N,).
    name : str, optional
        Identifier for the element, by default "thermoelastic".

    Attributes
    ----------
    alpha : torch.Tensor
        Thermal expansion coefficients, shape (N,).
    n_elems : int
        Number of elements.
    eps_th : torch.Tensor
        Thermal strain tensor per element, shape (N, 3, 3).
    I : torch.Tensor
        Identity tensor (broadcasted to N), shape (N, 3, 3).
    name : str
        Element name.
    """

    def __init__(self, alpha, name="thermoelastic"):
        self.alpha = alpha
        self.name = name
        self.n_elems = self.alpha.shape[0]
        self.eps_th = to.zeros((self.n_elems, 3, 3))
        self.I = to.eye(3, dtype=to.float64).unsqueeze(0).repeat(self.n_elems, 1, 1)

    def compute_eps_th(self, dT_DG_vec):
        """
        Compute thermal strain from a temperature increment.

        Parameters
        ----------
        dT_DG_vec : torch.Tensor
            Temperature increment per element (N,) or broadcastable to (N,).

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.eps_th = alpha * dT * I`.
        """
        self.eps_th = self.alpha[:, None, None] * dT_DG_vec[:, None, None] * self.I
