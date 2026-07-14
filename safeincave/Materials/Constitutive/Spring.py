# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
import torch as to
from ...Utils import dotdot_torch


class Spring:
    """
    Linear isotropic elastic element in Voigt notation.

    Parameters
    ----------
    E : torch.Tensor
        Young's modulus per element, shape (N,).
    nu : torch.Tensor
        Poisson's ratio per element, shape (N,).
    name : str, optional
        Element name, by default "spring".

    Attributes
    ----------
    E, nu : torch.Tensor
        Material parameters, shape (N,).
    n_elems : int
        Number of elements.
    C, C_inv : torch.Tensor
        Stiffness and its inverse in tensorial Voigt form, shape (N, 6, 6).
    C_tilde, C_tilde_inv : torch.Tensor
        Deviatoric stiffness and inverse, shape (N, 6, 6).
    K : torch.Tensor
        Bulk modulus per element, shape (N,).
    eps_e : torch.Tensor
        Elastic strain tensor per element, shape (N, 3, 3).
    name : str
        Element name.
    """

    def __init__(self, E: to.Tensor, nu: to.Tensor, name: str = "spring") -> None:
        self.E = E
        self.nu = nu
        self.name = name
        self.n_elems = self.E.shape[0]
        self.eps_e = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

    def initialize(self) -> None:
        """
        Build stiffness operators and bulk modulus.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `C`, `C_inv`, `C_tilde`, `C_tilde_inv`, and `K`.
        """
        self.C = self.__compute_C(self.n_elems, self.nu, self.E)
        self.C_inv = self.__compute_C_inv(self.C)
        self.C_tilde = self.__compute_C_tilde(self.n_elems, self.nu, self.E)
        self.C_tilde_inv = self.__compute_C_tilde_inv(self.n_elems, self.nu, self.E)
        self.K = self.E / (3 * (1 - 2 * self.nu))

    def compute_eps_e(self, stress: to.Tensor) -> None:
        """
        Compute elastic strain from stress using `C_inv`.

        Parameters
        ----------
        stress : torch.Tensor
            Cauchy stress tensor per element, shape (N, 3, 3).

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.eps_e` (N, 3, 3).
        """
        self.eps_e = dotdot_torch(self.C_inv, stress)

    def __compute_C(self, n_elems: int, nu: to.Tensor, E: to.Tensor) -> to.Tensor:
        """
        Construct isotropic stiffness in tensorial Voigt form.

        Parameters
        ----------
        n_elems : int
            Number of elements.
        nu : torch.Tensor
            Poisson's ratio, shape (N,).
        E : torch.Tensor
            Young's modulus, shape (N,).

        Returns
        -------
        torch.Tensor
            Stiffness matrix `(N, 6, 6)` with shear terms `2G` on the diagonal
            of the shear block (tensorial, not engineering).
        """
        C = to.zeros((n_elems, 6, 6), dtype=to.float64)
        a0 = E / ((1 + nu) * (1 - 2 * nu))
        C[:, 0, 0] = a0 * (1 - nu)
        C[:, 1, 1] = a0 * (1 - nu)
        C[:, 2, 2] = a0 * (1 - nu)
        C[:, 3, 3] = a0 * (1 - 2 * nu)
        C[:, 4, 4] = a0 * (1 - 2 * nu)
        C[:, 5, 5] = a0 * (1 - 2 * nu)
        C[:, 0, 1] = C[:, 1, 0] = C[:, 0, 2] = C[:, 2, 0] = C[:, 2, 1] = C[:, 1, 2] = (
            a0 * nu
        )
        return C

    def __compute_C_inv(self, C: to.Tensor) -> to.Tensor:
        """
        Invert the stiffness matrix per element.

        Parameters
        ----------
        C : torch.Tensor
            Stiffness matrix, shape (N, 6, 6).

        Returns
        -------
        torch.Tensor
            Element-wise inverse, shape (N, 6, 6).
        """
        return to.linalg.inv(self.C)

    def __compute_C_tilde(self, n_elems: int, nu: to.Tensor, E: to.Tensor) -> to.Tensor:
        """
        Construct deviatoric stiffness (2G on all six Voigt diagonals).

        Parameters
        ----------
        n_elems : int
        nu : torch.Tensor
        E : torch.Tensor

        Returns
        -------
        torch.Tensor
            `(N, 6, 6)` with diagonal entries `2G`, zeros elsewhere.
        """
        G = E / (2 * (1 + nu))
        C_tilde = to.zeros((n_elems, 6, 6), dtype=to.float64)
        C_tilde[:, 0, 0] = 2 * G
        C_tilde[:, 1, 1] = 2 * G
        C_tilde[:, 2, 2] = 2 * G
        C_tilde[:, 3, 3] = 2 * G
        C_tilde[:, 4, 4] = 2 * G
        C_tilde[:, 5, 5] = 2 * G
        return C_tilde

    def __compute_C_tilde_inv(self, n_elems: int, nu: to.Tensor, E: to.Tensor) -> to.Tensor:
        G = E / (2 * (1 + nu))
        C_tilde_inv = to.zeros((n_elems, 6, 6), dtype=to.float64)
        C_tilde_inv[:, 0, 0] = 1 / (2 * G)
        C_tilde_inv[:, 1, 1] = 1 / (2 * G)
        C_tilde_inv[:, 2, 2] = 1 / (2 * G)
        C_tilde_inv[:, 3, 3] = 1 / (2 * G)
        C_tilde_inv[:, 4, 4] = 1 / (2 * G)
        C_tilde_inv[:, 5, 5] = 1 / (2 * G)
        return C_tilde_inv
