# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
import torch as to
import numpy as np
from .Utils import MPa


class Material:
    """
    Composite material model that aggregates elastic, thermoelastic,
    and non-elastic (e.g., viscoelastic/viscoplastic) elements.

    The class stores element-wise stiffness operators in Voigt form
    and assembles effective operators used during constitutive updates.

    Parameters
    ----------
    n_elems : int
        Number of finite elements (batch size).

    Attributes
    ----------
    n_elems : int
        Number of elements.
    elems_ne : list[NonElasticElement]
        Collection of non-elastic elements contributing to non-elastic response.
    elems_th : list[Thermoelastic]
        Collection of thermoelastic contributors (thermal strain).
    elems_e : list[Spring]
        Collection of elastic contributors (linear isotropic springs).
    C_inv, C : torch.Tensor
        Element-wise stiffness (and inverse) in tensorial Voigt form, shape (N, 6, 6).
    C_tilde_inv, C_tilde : torch.Tensor
        Element-wise deviatoric stiffness (and inverse), shape (N, 6, 6).
    G, B : torch.Tensor
        Assembled non-elastic tangent-like (G) and state variable (B) operators, shapes
        (N, 6, 6) and (N, 3, 3).
    IT, T : torch.Tensor
        Volumetric coupling tensors assembled from non-elastic elements.
    B_vol, T_vol : torch.Tensor
        Element-wise volumetric parts (scalars), shape (N,).
    G_tilde, B_tilde : torch.Tensor
        Deviatoric parts of G and B, shapes (N, 6, 6) and (N, 3, 3).
    CT, CT_tilde : torch.Tensor
        Effective consistent tangents after time integration, shapes (N, 6, 6).
    density, cp, k, alpha_th : torch.Tensor
        Optional material properties (per element).

    Notes
    -----
    - Voigt ordering is assumed to be `[xx, yy, zz, xy, xz, yz]` with
      **tensorial shear** convention (no engineering factors).
    """

    def __init__(self, n_elems: int):
        self.n_elems = n_elems
        self.elems_ne = []
        self.elems_th = []
        self.elems_e = []

        self.C_inv = to.zeros((n_elems, 6, 6), dtype=to.float64)
        self.C = to.zeros((n_elems, 6, 6), dtype=to.float64)

        self.C_tilde_inv = to.zeros((n_elems, 6, 6), dtype=to.float64)
        self.C_tilde = to.zeros((n_elems, 6, 6), dtype=to.float64)

    def set_density(self, density: to.Tensor) -> None:
        """
        Set mass density per element.

        Parameters
        ----------
        density : torch.Tensor
            1D tensor of shape (N,) with densities.
        """
        self.density = density

    def set_specific_heat_capacity(self, cp: to.Tensor) -> None:
        """
        Set specific heat capacity per element.

        Parameters
        ----------
        cp : torch.Tensor
            1D tensor of shape (N,) with specific heat capacities.
        """
        self.cp = cp

    def set_thermal_conductivity(self, k: to.Tensor) -> None:
        """
        Set thermal conductivity per element.

        Parameters
        ----------
        k : torch.Tensor
            1D tensor of shape (N,) with conductivities.
        """
        self.k = k

    def set_thermal_expansion(self, alpha_th: to.Tensor) -> None:
        """
        Set coefficient of thermal expansion per element.

        Parameters
        ----------
        alpha_th : torch.Tensor
            1D tensor of shape (N,) with linear thermal expansion coefficients.
        """
        self.alpha_th = alpha_th

    def add_to_elastic(self, elem: Spring):
        """
        Add an elastic (linear isotropic) contributor and accumulate stiffness.

        Parameters
        ----------
        elem : Spring
            Elastic element. Its `initialize()` is called inside.

        Side Effects
        ------------
        - Updates `C`, `C_inv`, `C_tilde`, `C_tilde_inv` by addition.
        - Stores `K`, `E`, and shear modulus estimate `ShearMod`.
        - Appends `elem` to `elems_e`.
        """
        elem.initialize()
        self.C_inv += elem.C_inv
        self.C += elem.C
        self.C_tilde_inv += elem.C_tilde_inv
        self.C_tilde += elem.C_tilde
        self.elems_e.append(elem)
        self.K = elem.K
        self.E = elem.E
        self.ShearMod = 3 * self.K * self.E / (9 * self.K - self.E)

    def add_to_non_elastic(self, elem) -> None:
        """
        Add a non-elastic element contributor.

        Parameters
        ----------
        elem : NonElasticElement
            Inelastic mechanism (e.g., creep, viscoplasticity).
        """
        self.elems_ne.append(elem)

    def add_to_thermoelastic(self, elem: Thermoelastic) -> None:
        """
        Add a thermoelastic contributor.

        Parameters
        ----------
        elem : Thermoelastic
            Provides thermal strain contributions.
        """
        self.elems_th.append(elem)

    def compute_G_B(
        self, stress: to.Tensor, dt: float, theta: float, T: to.Tensor
    ) -> None:
        """
        Assemble non-elastic operators G and B over all inelastic elements.

        Parameters
        ----------
        stress : torch.Tensor
            Current Cauchy stress per element, shape (N, 3, 3).
        dt : float
            Time step size.
        theta : float
            Time integration parameter: 0 for fully implicit, 0.5 for Crank-Nicolson, 1 for explicit.
        T : torch.Tensor
            Temperature per element (shape (N,) or compatible).

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.G` (N,6,6) and `self.B` (N,3,3) as sums of element contributions.
        """
        self.G = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
        self.B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        for elem_ne in self.elems_ne:
            elem_ne.compute_G_B(stress, dt, theta, T)
            self.G += elem_ne.G
            self.B += elem_ne.B

    def compute_T_IT(self) -> None:
        """
        Assemble volumetric coupling tensors T and IT from inelastic elements.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.T` (N,3,3) and `self.IT` (N,6,6) as sums of element contributions.
        """
        self.IT = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
        self.T = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        for elem_ne in self.elems_ne:
            elem_ne.compute_T_IT()
            self.IT += elem_ne.IT
            self.T += elem_ne.T

    def compute_Bvol_Tvol(self, stress: to.Tensor, dt: float) -> None:
        """
        Compute volumetric parts of B and T.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3). (Not used directly here.)
        dt : float
            Time step size. (Not used directly here.)

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.B_vol` and `self.T_vol` (shape (N,)) from element contributions.
        """
        self.B_vol = to.zeros(self.n_elems, dtype=to.float64)
        self.T_vol = to.zeros(self.n_elems, dtype=to.float64)
        for elem_ne in self.elems_ne:
            elem_ne.compute_Bvol_Tvol()
            self.B_vol += elem_ne.B_vol
            self.T_vol += elem_ne.T_vol

    def compute_Gtilde_Btilde(self, stress: to.Tensor, dt: float) -> None:
        """
        Compute deviatoric parts of G and B.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3). (Not used directly here.)
        dt : float
            Time step size. (Not used directly here.)

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.G_tilde` and `self.B_tilde` (N,6,6) and (N,3,3).
        """
        self.G_tilde = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
        self.B_tilde = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        for elem_ne in self.elems_ne:
            elem_ne.compute_Gtilde_Btilde()
            self.G_tilde += elem_ne.G_tilde
            self.B_tilde += elem_ne.B_tilde

    def compute_CT(self, dt: float, theta: float) -> None:
        """
        Compute consistent tangent `CT = (C_inv + dt*(1-theta)*G)^{-1}`.

        Parameters
        ----------
        dt : float
            Time step size.
        theta : float
            Time integration parameter: 0 for fully implicit, 0.5 for Crank-Nicolson, 1 for explicit.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.CT` (N, 6, 6).
        """
        self.CT = to.linalg.inv(self.C_inv + dt * (1 - theta) * self.G)

    def compute_CT_tilde(self, dt: float, theta: float) -> None:
        """
        Compute deviatoric consistent tangent `CT_tilde`.

        Parameters
        ----------
        dt : float
            Time step size.
        theta : float
            Time integration parameter in [0, 1].

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.CT_tilde` (N, 6, 6).
        """
        self.CT_tilde = to.linalg.inv(
            self.C_tilde_inv + dt * (1 - theta) * self.G_tilde
        )


# Import all constitutive models from the ConstitutiveModels folder
from .ConstitutiveModels import (
    NonElasticElement,
    Thermoelastic,
    Spring,
    Viscoelastic,
    LinearDashpot,
    DislocationCreep,
    PressureSolutionCreep,
    ViscoplasticDesai,
    MunsonDawsonCreep,
    ModifiedCamClayViscoplastic,
)

__all__ = [
    "Material",
    "NonElasticElement",
    "Thermoelastic",
    "Spring",
    "Viscoelastic",
    "LinearDashpot",
    "DislocationCreep",
    "PressureSolutionCreep",
    "ViscoplasticDesai",
    "MunsonDawsonCreep",
    "ModifiedCamClayViscoplastic",
]
