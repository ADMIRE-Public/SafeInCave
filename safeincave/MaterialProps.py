# Copyright 2025 The safeincave community.
#
# This file is part of safeincave.
#
# Licensed under the GNU GENERAL PUBLIC LICENSE, Version 3 (the "License"); you may not
# use this file except in compliance with the License.  You may obtain a copy
# of the License at
#
#     https://spdx.org/licenses/GPL-3.0-or-later.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations under
# the License.
from __future__ import annotations
from abc import ABC, abstractmethod
import torch as to
import numpy as np
from .Utils import dotdot_torch, MPa
# try:
#     from .Utils import dotdot_torch, MPa
# except:
#     from Utils import dotdot_torch, MPa

class Material():
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
		self.ShearMod = 3*self.K*self.E/(9*self.K - self.E)

	def add_to_non_elastic(self, elem: NonElasticElement) -> None:
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

	def compute_G_B(self, stress: to.Tensor, dt: float, theta: float, T: to.Tensor) -> None:
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
		self.CT = to.linalg.inv(self.C_inv + dt*(1-theta)*self.G)

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
		self.CT_tilde = to.linalg.inv(self.C_tilde_inv + dt*(1-theta)*self.G_tilde)


class Thermoelastic():
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
		self.eps_th = self.alpha[:,None,None]*dT_DG_vec[:,None,None]*self.I


class Spring():
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
	def __init__(self, E, nu, name="spring"):
		self.E = E
		self.nu = nu
		self.name = name
		self.n_elems = self.E.shape[0]
		self.eps_e = to.tensor((self.n_elems, 3, 3), dtype=to.float64)

	def initialize(self):
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
		self.K = self.E/(3*(1 - 2*self.nu))

	def compute_eps_e(self, stress):
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

	def __compute_C(self, n_elems, nu, E):
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
		a0 = E/((1 + nu)*(1 - 2*nu))
		C[:,0,0] = a0*(1 - nu)
		C[:,1,1] = a0*(1 - nu)
		C[:,2,2] = a0*(1 - nu)
		C[:,3,3] = a0*(1 - 2*nu)
		C[:,4,4] = a0*(1 - 2*nu)
		C[:,5,5] = a0*(1 - 2*nu)
		C[:,0,1] = C[:,1,0] = C[:,0,2] = C[:,2,0] = C[:,2,1] = C[:,1,2] = a0*nu
		return C

	def __compute_C_inv(self, C):
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

	def __compute_C_tilde(self, n_elems, nu, E):
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
		G = E/(2*(1 + nu))
		C_tilde = to.zeros((n_elems, 6, 6), dtype=to.float64)
		C_tilde[:,0,0] = 2*G
		C_tilde[:,1,1] = 2*G
		C_tilde[:,2,2] = 2*G
		C_tilde[:,3,3] = 2*G
		C_tilde[:,4,4] = 2*G
		C_tilde[:,5,5] = 2*G
		return C_tilde

	def __compute_C_tilde_inv(self, n_elems, nu, E):
		G = E/(2*(1 + nu))
		C_tilde_inv = to.zeros((n_elems, 6, 6), dtype=to.float64)
		C_tilde_inv[:,0,0] = 1/(2*G)
		C_tilde_inv[:,1,1] = 1/(2*G)
		C_tilde_inv[:,2,2] = 1/(2*G)
		C_tilde_inv[:,3,3] = 1/(2*G)
		C_tilde_inv[:,4,4] = 1/(2*G)
		C_tilde_inv[:,5,5] = 1/(2*G)
		return C_tilde_inv



class NonElasticElement(ABC):
    """
    Abstract base for inelastic mechanisms (e.g., viscoelasticity,
    dislocation creep, viscoplasticity). Provides common storage
    and utility updates for internal variables.

    Parameters
    ----------
    n_elems : int
        Number of elements.

    Attributes
    ----------
    n_elems : int
        Number of elements.
    eps_ne_rate, eps_ne_rate_old : torch.Tensor
        Current and previous non-elastic strain rate, shape (N, 3, 3).
    eps_ne_old, eps_ne_k : torch.Tensor
        Non-elastic strain at old and current time, shape (N, 3, 3).
    B : torch.Tensor
        State variable term (N, 3, 3) assembled in `compute_G_B`.
    G : torch.Tensor
        Tangent-like operator (N, 6, 6) assembled in `compute_G_B`.
    """
    def __init__(self, n_elems):
    	self.n_elems = n_elems
    	self.eps_ne_rate = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
    	self.eps_ne_rate_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
    	self.eps_ne_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
    	self.eps_ne_k = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
    	self.B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
    	self.G = to.zeros((self.n_elems, 6, 6), dtype=to.float64)

    @abstractmethod
    def compute_eps_ne_rate(self, stress_vec: to.Tensor, phi1: float, Temp: to.Tensor, return_eps_ne: bool=False) -> None:
    	pass

    def increment_internal_variables(self, *args) -> None:
    	pass

    def update_internal_variables(self, *args) -> None:
    	pass

    def compute_eps_ne_k(self, phi1: float, phi2: float) -> None:
        """
        Predictor for non-elastic strain at the previous iteration k.

        Parameters
        ----------
        phi1 : float
            Typically `phi1=dt*theta`.
        phi2 : float
            Typically `phi2=dt*(1-theta)`.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.eps_ne_k`.
        """
        self.eps_ne_k = self.eps_ne_old + phi1*self.eps_ne_rate_old + phi2*self.eps_ne_rate

    def update_eps_ne_old(self, stress: to.Tensor, stress_k: to.Tensor, phi2: float) -> None:
        """
        Update non-elastic strain from previous time step.

        Parameters
        ----------
        stress : torch.Tensor
            Stress at current iteration k+1, shape (N, 3, 3).
        stress_k : torch.Tensor
            Stress at previous iteration k, shape (N, 3, 3).
        phi2 : float
            Typically `phi2=dt*(1-theta)`.

        Returns
        -------
        None

        Side Effects
        ------------
        Updates `self.eps_ne_old`.
        """
        self.eps_ne_old = self.eps_ne_k + phi2*dotdot_torch(self.G, stress - stress_k) - phi2*self.B

    def update_eps_ne_rate_old(self) -> None:
        """
        Update the current non-nelastic strain rate to the previous time step.

        Returns
        -------
        None
        """
        self.eps_ne_rate_old = self.eps_ne_rate.clone()

    def compute_E(self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor) -> None:
        """
        Finite-difference approximation of the 6×6 operator E = d(eps_ne)/d(stress).

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).
        dt : float
            Time step.
        theta : float
            Time integration parameter.
        Temp : torch.Tensor
            Temperature per element.

        Returns
        -------
        torch.Tensor
            Operator `E` with shape (N, 6, 6).
        """
        phi1 = dt*theta
        EPSILON = 1e-2
        E = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
        stress_eps = stress.clone()
        c1 = 1.0
        c2 = 2.0
        magic_indexes = [(0,0,0,c1), (1,1,1,c1), (2,2,2,c1), (0,1,3,c2), (0,2,4,c2), (1,2,5,c2)]
        for i, j, k, phi in magic_indexes:
        	stress_eps[:,i,j] += EPSILON
        	eps_A = self.compute_eps_ne_rate(stress_eps, phi1, Temp, return_eps_ne=True)
        	stress_eps[:,i,j] -= EPSILON
        	stress_eps[:,i,j] -= EPSILON
        	eps_B = self.compute_eps_ne_rate(stress_eps, phi1, Temp, return_eps_ne=True)
        	stress_eps[:,i,j] += EPSILON
        	E[:,:,k] = phi*(eps_A[:,[0,1,2,0,0,1],[0,1,2,1,2,2]] - eps_B[:,[0,1,2,0,0,1],[0,1,2,1,2,2]]) / (2*EPSILON)
        return E

    def compute_B_and_H_over_h(self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor) -> None:
        """
        Compute state variable term `B` and linearization term `H/h`.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).
        dt : float
            Time step.
        theta : float
            Time integration parameter.
        Temp : torch.Tensor
            Temperature per element.

        Returns
        -------
        B : torch.Tensor
            Driving term, shape (N, 3, 3).
        H_over_h : torch.Tensor
            Linearization ratio, shape (N, 6, 6).

        Notes
        -----
        Default implementation returns zeros; subclasses override.
        """
        B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        H_over_h = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
        return B, H_over_h

    def compute_G_B(self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor) -> None:
        """
        Assemble `G` and `B` for the element based on `E` and `H/h`.

        Parameters
        ----------
        stress : torch.Tensor
        dt : float
        theta : float
        Temp : torch.Tensor

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.B` and `self.G`.
        """
        self.B, H_over_h = self.compute_B_and_H_over_h(stress, dt, theta, Temp)
        E = self.compute_E(stress, dt, theta, Temp)
        self.G = E.clone() - H_over_h.clone()

    def compute_T_IT(self) -> None:
        """
        Build volumetric coupling tensors `T` (3×3) and `IT` (6×6) from `G`.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.T` and `self.IT`.
        """
        self.T = to.zeros((self.n_elems, 3, 3))
        self.T[:,0,0] = self.G[:,0,0] + self.G[:,1,0] + self.G[:,2,0]
        self.T[:,1,1] = self.G[:,0,1] + self.G[:,1,1] + self.G[:,2,1]
        self.T[:,2,2] = self.G[:,0,2] + self.G[:,1,2] + self.G[:,2,2]
        self.T[:,1,0] = self.T[:,0,1] = (self.G[:,0,3] + self.G[:,1,3] + self.G[:,2,3])/2
        self.T[:,2,0] = self.T[:,0,2] = (self.G[:,0,4] + self.G[:,1,4] + self.G[:,2,4])/2
        self.T[:,2,1] = self.T[:,1,2] = (self.G[:,0,5] + self.G[:,1,5] + self.G[:,2,5])/2

        self.IT = to.zeros((self.n_elems, 6, 6))
        self.IT[:,0,0] = self.T[:,0,0]
        self.IT[:,0,1] = self.T[:,1,1]
        self.IT[:,0,2] = self.T[:,2,2]
        self.IT[:,0,3] = self.T[:,0,1] + self.T[:,1,0]
        self.IT[:,0,4] = self.T[:,0,2] + self.T[:,2,0]
        self.IT[:,0,5] = self.T[:,1,2] + self.T[:,2,1]
        self.IT[:,1,:] = self.IT[:,0,:]
        self.IT[:,2,:] = self.IT[:,0,:]

    def compute_Bvol_Tvol(self) -> None:
        """
        Compute volumetric parts of `B` and `T`.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.B_vol` and `self.T_vol` as traces of `B` and `T`.
        """
        self.T_vol = to.einsum("bii->b", self.T)
        self.B_vol = to.einsum("bii->b", self.B)

    def compute_Gtilde_Btilde(self) -> None:
        """
        Compute deviatoric parts `G_tilde` and `B_tilde`.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.G_tilde` and `self.B_tilde`.
        """
        I = to.eye(3).expand(self.n_elems, -1, -1)
        self.G_tilde = self.G - self.IT/3
        self.B_tilde = self.B - self.B_vol[:,None,None]*I/3

		



class Viscoelastic(NonElasticElement):
    """
    Kelvin–Voigt-type viscoelastic element.

    Parameters
    ----------
    eta : torch.Tensor
        Viscosity parameter per element, shape (N,).
    E : torch.Tensor
        Young's modulus per element, shape (N,).
    nu : torch.Tensor
        Poisson's ratio per element, shape (N,).
    name : str, optional
        Element name, by default "kelvin_voigt".

    Attributes
    ----------
    C1 : torch.Tensor
        Elastic stiffness in Voigt form, shape (N, 6, 6).
    eta, E, nu : torch.Tensor
        Material parameters, shape (N,).
    """
    def __init__(self, eta: to.Tensor, E: to.Tensor, nu: to.Tensor, name: bool="kelvin_voigt"):
        super().__init__(E.shape[0])
        self.eta = eta
        self.E = E
        self.nu = nu
        self.name = name

        # Assemble C1 tensor (n_elems, 6, 6)
        self.C1 = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
        a0 = self.E/((1 + self.nu)*(1 - 2*self.nu))
        self.C1[:,0,0] = a0*(1 - self.nu)
        self.C1[:,1,1] = a0*(1 - self.nu)
        self.C1[:,2,2] = a0*(1 - self.nu)
        self.C1[:,3,3] = a0*(1 - 2*self.nu)
        self.C1[:,4,4] = a0*(1 - 2*self.nu)
        self.C1[:,5,5] = a0*(1 - 2*self.nu)
        self.C1[:,0,1] = self.C1[:,1,0] = self.C1[:,0,2] = self.C1[:,2,0] = self.C1[:,2,1] = self.C1[:,1,2] = a0*self.nu

    def compute_eps_ne_rate(self, stress_vec: to.Tensor, phi1: float, Temp: to.Tensor, return_eps_ne: bool=False):
        """
        Compute viscoelastic strain rate (Kelvin–Voigt form).

        Parameters
        ----------
        stress_vec : torch.Tensor
            Stress tensor per element, shape (N, 3, 3).
        phi1 : float
            Time integration factor (dt*theta).
        Temp : torch.Tensor
            Temperature per element (unused here).
        return_eps_ne : bool, default=False
            If True, return the rate; else store it.

        Returns
        -------
        None or torch.Tensor
            (N, 3, 3) if `return_eps_ne=True`, else `None`.
        """
        eps_ne_rate = dotdot_torch(self.G, stress_vec - dotdot_torch(self.C1, self.eps_ne_old + phi1*self.eps_ne_rate_old))
        if return_eps_ne:
        	return eps_ne_rate.clone()
        else:
        	self.eps_ne_rate = eps_ne_rate.clone()

    def compute_E(self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor) -> to.Tensor:
        """
        Closed-form 6×6 operator for viscoelasticity:
        `E = (eta*I + phi2*C1)^{-1}`.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3). (Unused here.)
        dt : float
            Time step.
        theta : float
            Time integration parameter.
        Temp : torch.Tensor
            Temperature per element (unused).

        Returns
        -------
        torch.Tensor
            `E` with shape (N, 6, 6).
        """
        phi2 = dt*(1 - theta)
        I = to.eye(6, dtype=to.float64).unsqueeze(0).repeat(self.n_elems, 1, 1)
        E = to.linalg.inv(self.eta[:,None,None]*I + phi2*self.C1)
        return E




class DislocationCreep(NonElasticElement):
    """
    Power-law dislocation creep: :math:`\\dot\\varepsilon_{ne}
    = A\\,\\exp(-Q/(RT))\\,\\q^{n-1}\\,\\mathbf{s}`.

    Parameters
    ----------
    A : torch.Tensor
        Pre-exponential factor per element, shape (N,).
    Q : torch.Tensor
        Activation energy per element, shape (N,).
    n : torch.Tensor
        Stress exponent per element, shape (N,).
    name : str, optional
        Element name, by default "creep".

    Attributes
    ----------
    R : float
        Gas constant used (8.32).
    A, Q, n : torch.Tensor
        Material parameters, shape (N,).
    """
    def __init__(self, A: to.Tensor, Q: to.Tensor, n: to.Tensor, name: bool="creep"):
        super().__init__(A.shape[0])
        self.R = 8.32
        self.Q = Q
        self.A = A
        self.n = n
        self.name = name

    def compute_eps_ne_rate(self, stress_vec: to.Tensor, phi1: float, Temp: to.Tensor, return_eps_ne: bool=False):
        """
        Compute creep strain rate from current stress.

        Parameters
        ----------
        stress_vec : torch.Tensor
            Stress tensor per element, shape (N, 3, 3).
        phi1 : float
            Time integration factor (unused here).
        Temp : torch.Tensor
            Temperature per element (N,) or broadcastable.
        return_eps_ne : bool, default=False
            If True, return the rate; else store it.

        Returns
        -------
        None or torch.Tensor
            (N, 3, 3) if `return_eps_ne=True`, else `None`.
        """
        s_xx = stress_vec[:,0,0]
        s_yy = stress_vec[:,1,1]
        s_zz = stress_vec[:,2,2]
        s_xy = stress_vec[:,0,1]
        s_xz = stress_vec[:,0,2]
        s_yz = stress_vec[:,1,2]

        sigma_mean = (s_xx + s_yy + s_zz) / 3
        dev = stress_vec.clone()
        dev[:,0,0] = s_xx - sigma_mean
        dev[:,1,1] = s_yy - sigma_mean
        dev[:,2,2] = s_zz - sigma_mean

        q_vm = to.sqrt( 0.5*( (s_xx - s_yy)**2 + (s_xx - s_zz)**2 + (s_yy - s_zz)**2 + 6*(s_xy**2 + s_xz**2 + s_yz**2) ) )

        A_bar = self.A*to.exp(-self.Q/self.R/Temp)*q_vm**(self.n - 1)
        eps_rate = A_bar[:,None,None]*dev
        if return_eps_ne:
        	return eps_rate
        else:
        	self.eps_ne_rate = eps_rate



class ViscoplasticDesai(NonElasticElement):
    """
    Viscoplastic model of Desai-type with hardening (state) variable `alpha`.

    Parameters
    ----------
    mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0 : torch.Tensor
        Model parameters per element, shape (N,).
    name : str, optional
        Element name, by default "desai".

    Attributes
    ----------
    alpha : torch.Tensor
        Current hardening variable per element, shape (N,).
    Fvp : torch.Tensor
        Current value of the yield function per element, shape (N,).
    qsi, qsi_old : torch.Tensor
        Accumulated viscoplastic strain measure and its previous value, shape (N,).
    P : torch.Tensor
        Sensitivity of the residue to stress, shape (N, 3, 3).
    r, h : torch.Tensor
        Residue and its derivative w.r.t. `alpha`, shapes (N,) and (N,).
    """
    def __init__(self,
						mu_1: to.Tensor,
						N_1: to.Tensor,
						a_1: to.Tensor,
						eta: to.Tensor,
						n: to.Tensor,
						beta_1: to.Tensor,
						beta: to.Tensor,
						m: to.Tensor,
						gamma: to.Tensor,
						sigma_t: to.Tensor,
						alpha_0: to.Tensor,
						name: bool="desai"):
        super().__init__(mu_1.shape[0])
        self.name = name
        self.mu_1 = mu_1
        self.N_1 = N_1
        self.a_1 = a_1
        self.eta = eta
        self.n = n
        self.beta_1 = beta_1
        self.beta = beta
        self.m = m
        self.gamma = gamma
        self.sigma_t = sigma_t
        self.alpha_0 = alpha_0
        self.F_0 = 1.0
        self.n_elems = self.alpha_0.shape[0]
        self.alpha = self.alpha_0.clone()
        self.Fvp = to.zeros(self.n_elems, dtype=to.float64)
        self.qsi = to.zeros(self.n_elems, dtype=to.float64)
        self.qsi_old = to.zeros(self.n_elems, dtype=to.float64)

    def compute_residue(self, eps_rate: to.Tensor, alpha: to.Tensor, dt: float) -> to.Tensor:
        """
        Residue of the implicit hardening equation.

        Parameters
        ----------
        eps_rate : torch.Tensor
            Current inelastic strain rate, shape (N, 3, 3).
        alpha : torch.Tensor
            Hardening variable, shape (N,).
        dt : float
            Time step.

        Returns
        -------
        torch.Tensor
            Residue per element, shape (N,).

        Notes
        -----
        Updates `self.qsi` internally based on `eps_rate`.
        """
        self.qsi = self.qsi_old + to.sum(eps_rate**2, axis=(-2, -1))**0.5*dt
        return alpha - self.a_1 / (((self.a_1/self.alpha_0)**(1/self.eta) + self.qsi)**self.eta)

    def update_internal_variables(self) -> None:
        """
        Commit accumulated measure `qsi`.

        Returns
        -------
        None
        """
        self.qsi_old = self.qsi.clone()

    def increment_internal_variables(self, stress: to.Tensor, stress_k: to.Tensor, dt:float) -> None:
        """
        Increment hardening variable `alpha` using linearization.

        Parameters
        ----------
        stress : torch.Tensor
            End-of-step stress, shape (N, 3, 3).
        stress_k : torch.Tensor
            Intermediate stress, shape (N, 3, 3).
        dt : float
            Time step.

        Returns
        -------
        None

        Side Effects
        ------------
        Updates `self.alpha`.
        """
        delta_alpha = -(self.r + to.einsum('bij,bij->b', self.P, stress - stress_k))/self.h
        self.alpha += delta_alpha

    def compute_stress_invariants(self, s_xx: to.Tensor,
										s_yy: to.Tensor,
										s_zz: to.Tensor,
										s_xy: to.Tensor,
										s_xz: to.Tensor,
										s_yz: to.Tensor) -> tuple[to.Tensor, to.Tensor, to.Tensor, to.Tensor, to.Tensor, to.Tensor, to.Tensor, to.Tensor]:
        """
        Compute invariants (I1, I2, I3, J2, J3) and auxiliary quantities.

        Parameters
        ----------
        s_xx, s_yy, s_zz, s_xy, s_xz, s_yz : torch.Tensor
            Normal and shear components (MPa-scaled as provided), each shape (N,).

        Returns
        -------
        I1, I2, I3, J2, J3, Sr, I1_star, ind_J2_leq_0 : tuple of torch.Tensor
            Invariants and helper arrays; `ind_J2_leq_0` are indices where `J2 <= 0`.
        """
        I1 = s_xx + s_yy + s_zz
        I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
        I3 = s_xx*s_yy*s_zz + 2*s_xy*s_yz*s_xz - s_zz*s_xy**2 - s_xx*s_yz**2 - s_yy*s_xz**2
        J2 = (1/3)*I1**2 - I2
        J3 = (2/27)*I1**3 - (1/3)*I1*I2 + I3
        Sr = -(J3*np.sqrt(27))/(2*J2**1.5)

        # Check where J2 <= 0.0
        ind_J2_leq_0 = to.where(J2 <= 0.0)[0]

        # Sr will be nan if, J2=0.0. So, replace it by 0.0
        Sr[ind_J2_leq_0] = 0.0

        I1_star = I1 + self.sigma_t
        return I1, I2, I3, J2, J3, Sr, I1_star, ind_J2_leq_0

    def extract_stress_components(self, stress: to.Tensor) -> to.Tensor:
        """
        Extract and scale stress components from a 3×3 tensor.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).

        Returns
        -------
        tuple[torch.Tensor, ...]
            `(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)`, each shape (N,).
        """
        stress_vec = -stress
        s_xx = stress_vec[:,0,0]/MPa
        s_yy = stress_vec[:,1,1]/MPa
        s_zz = stress_vec[:,2,2]/MPa
        s_xy = stress_vec[:,0,1]/MPa
        s_xz = stress_vec[:,0,2]/MPa
        s_yz = stress_vec[:,1,2]/MPa
        return s_xx, s_yy, s_zz, s_xy, s_xz, s_yz

    def compute_Fvp(self, alpha, I1, J2, Sr):
        """
        Compute the Desai viscoplastic yield function value `Fvp`.

        Parameters
        ----------
        alpha : torch.Tensor
            Hardening variable per element, shape (N,).
        I1 : torch.Tensor
            First invariant, shape (N,).
        J2 : torch.Tensor
            Second deviatoric invariant, shape (N,).
        Sr : torch.Tensor
            Lode-related parameter, shape (N,).

        Returns
        -------
        torch.Tensor
            `Fvp` per element, shape (N,).
        """
        F1 = (alpha*I1**self.n - self.gamma*I1**2)
        F2 = (to.exp(self.beta_1*I1) - self.beta*Sr)
        Fvp = J2 + F1*F2**self.m
        return Fvp

    def compute_initial_hardening(self, stress: to.Tensor, Fvp_0=0.0) -> None:
        """
        Initialize `alpha` from a target `Fvp_0` and the current stress state.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).
        Fvp_0 : float, default=0.0
            Target initial value for `Fvp`.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.alpha_0`, `self.alpha`, and `self.Fvp`.
        """
        s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.extract_stress_components(stress)
        I1, I2, I3, J2, J3, Sr, I1_star, _ = self.compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)
        self.alpha_0 =  self.gamma*I1_star**(2-self.n) + (Fvp_0 - J2)*I1_star**(-self.n)*(to.exp(self.beta_1*I1_star) - self.beta*Sr)**(-self.m)
        self.alpha = self.alpha_0.clone()

        s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.extract_stress_components(stress)
        I1, I2, I3, J2, J3, Sr, I1_star, ind_J2_leq_0 = self.compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)
        self.Fvp = self.compute_Fvp(self.alpha, I1_star, J2, Sr)


    def compute_eps_ne_rate(self, stress: to.Tensor, phi1: float, Temp: to.Tensor, alpha=None, return_eps_ne=False):
        """
        Compute viscoplastic strain rate and optionally return it.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).
        phi1 : float
            Time integration factor (dt*theta).
        Temp : torch.Tensor
            Temperature per element (unused).
        alpha : torch.Tensor or None, optional
            Hardening variable override; if `None`, use `self.alpha`.
        return_eps_ne : bool, default=False
            If True, return rate; else store it.

        Returns
        -------
        None or torch.Tensor
            (N, 3, 3) if `return_eps_ne=True`, else `None`.

        Notes
        -----
        Also updates `self.Fvp` when `return_eps_ne=False`.
        """
        if alpha == None:
        	alpha = self.alpha

        s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.extract_stress_components(stress)
        I1, I2, I3, J2, J3, Sr, I1_star, ind_J2_leq_0 = self.compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)

        # Compute yield function
        Fvp = self.compute_Fvp(alpha, I1_star, J2, Sr)
        if not return_eps_ne:
        	self.Fvp = Fvp.clone()


        # Compute flow direction, i.e. d(Fvp)/d(stress)
        F1 = (-alpha*I1**self.n + self.gamma*I1**2)
        F2 = (to.exp(self.beta_1*I1) - self.beta*Sr)
        dF1_dI1 = 2*self.gamma*I1 - self.n*alpha*I1**(self.n-1)
        dF2m_dI1 = self.beta_1*self.m*to.exp(self.beta_1*I1)*F2**(self.m-1)
        dF_dI1 = -(dF1_dI1*F2**self.m + F1*dF2m_dI1)

        dF2_dJ2 = -(3*self.beta*J3*27**0.5)/(4*J2**(5/2))
        dF_dJ2 = 1 - F1*self.m*F2**(self.m-1)*dF2_dJ2
        dF_dJ3 = -self.m*F1*self.beta*np.sqrt(27)*F2**(self.m-1)/(2*J2**1.5)


        dI1_dSxx = 1.0
        dI1_dSyy = 1.0
        dI1_dSzz = 1.0
        dI1_dSxy = 0.0
        dI1_dSxz = 0.0
        dI1_dSyz = 0.0

        dI2_dSxx = s_yy + s_zz
        dI2_dSyy = s_xx + s_zz
        dI2_dSzz = s_xx + s_yy
        dI2_dSxy = -2*s_xy
        dI2_dSxz = -2*s_xz
        dI2_dSyz = -2*s_yz

        dI3_dSxx = s_yy*s_zz - s_yz**2
        dI3_dSyy = s_xx*s_zz - s_xz**2
        dI3_dSzz = s_xx*s_yy - s_xy**2
        dI3_dSxy = 2*(s_xz*s_yz - s_zz*s_xy)
        dI3_dSxz = 2*(s_xy*s_yz - s_yy*s_xz)
        dI3_dSyz = 2*(s_xz*s_xy - s_xx*s_yz)

        dJ2_dI1 = (2/3)*I1
        dJ2_dI2 = -1.0

        dJ2_dSxx = dJ2_dI1*dI1_dSxx + dJ2_dI2*dI2_dSxx
        dJ2_dSyy = dJ2_dI1*dI1_dSyy + dJ2_dI2*dI2_dSyy
        dJ2_dSzz = dJ2_dI1*dI1_dSzz + dJ2_dI2*dI2_dSzz
        dJ2_dSxy = dJ2_dI1*dI1_dSxy + dJ2_dI2*dI2_dSxy
        dJ2_dSxz = dJ2_dI1*dI1_dSxz + dJ2_dI2*dI2_dSxz
        dJ2_dSyz = dJ2_dI1*dI1_dSyz + dJ2_dI2*dI2_dSyz

        dJ3_dI1 = (2/9)*I1**2 - (1/3)*I2
        dJ3_dI2 = -(1/3)*I1
        dJ3_dI3 = 1.0

        dJ3_dSxx = dJ3_dI1*dI1_dSxx + dJ3_dI2*dI2_dSxx + dJ3_dI3*dI3_dSxx
        dJ3_dSyy = dJ3_dI1*dI1_dSyy + dJ3_dI2*dI2_dSyy + dJ3_dI3*dI3_dSyy
        dJ3_dSzz = dJ3_dI1*dI1_dSzz + dJ3_dI2*dI2_dSzz + dJ3_dI3*dI3_dSzz
        dJ3_dSxy = dJ3_dI1*dI1_dSxy + dJ3_dI2*dI2_dSxy + dJ3_dI3*dI3_dSxy
        dJ3_dSxz = dJ3_dI1*dI1_dSxz + dJ3_dI2*dI2_dSxz + dJ3_dI3*dI3_dSxz
        dJ3_dSyz = dJ3_dI1*dI1_dSyz + dJ3_dI2*dI2_dSyz + dJ3_dI3*dI3_dSyz

        dQdS_00 = dF_dI1*dI1_dSxx + dF_dJ2*dJ2_dSxx + dF_dJ3*dJ3_dSxx
        dQdS_11 = dF_dI1*dI1_dSyy + dF_dJ2*dJ2_dSyy + dF_dJ3*dJ3_dSyy
        dQdS_22 = dF_dI1*dI1_dSzz + dF_dJ2*dJ2_dSzz + dF_dJ3*dJ3_dSzz
        dQdS_01 = dQdS_10 = dF_dI1*dI1_dSxy + dF_dJ2*dJ2_dSxy + dF_dJ3*dJ3_dSxy
        dQdS_02 = dQdS_20 = dF_dI1*dI1_dSxz + dF_dJ2*dJ2_dSxz + dF_dJ3*dJ3_dSxz
        dQdS_12 = dQdS_21 = dF_dI1*dI1_dSyz + dF_dJ2*dJ2_dSyz + dF_dJ3*dJ3_dSyz

        # Initialize viscoplastic direction
        dQdS = to.zeros_like(stress, dtype=to.float64)
        dQdS[:,0,0] = dQdS_00
        dQdS[:,1,1] = dQdS_11
        dQdS[:,2,2] = dQdS_22
        dQdS[:,1,0] = dQdS[:,0,1] = dQdS_01
        dQdS[:,2,0] = dQdS[:,0,2] = dQdS_02
        dQdS[:,2,1] = dQdS[:,1,2] = dQdS_12

        # Wherever J2=0, make viscoplasticity to be zero
        dQdS[ind_J2_leq_0,:,:] = 0.0

        # Calculate strain rate
        ramp_idx = to.where(Fvp > 0)[0]
        lmbda = to.zeros(self.n_elems, dtype=to.float64)
        if len(ramp_idx) != 0:
        	lmbda[ramp_idx] = self.mu_1[ramp_idx]*(Fvp[ramp_idx]/self.F_0)**self.N_1[ramp_idx]
        eps_vp_rate = -dQdS*lmbda[:, None, None]

        if return_eps_ne:
        	return eps_vp_rate
        else:
            self.eps_ne_rate = eps_vp_rate


    def compute_B_and_H_over_h(self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor) -> tuple[to.Tensor, to.Tensor]:
        """
        Compute `B` and `H/h` via perturbations of `alpha` and stress.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).
        dt : float
            Time step.
        theta : float
            Time integration parameter.
        Temp : torch.Tensor
            Temperature per element (unused).

        Returns
        -------
        B : torch.Tensor
            Driving term, shape (N, 3, 3).
        H_over_h : torch.Tensor
            Linearization ratio, shape (N, 6, 6).

        Notes
        -----
        Uses finite differences to approximate sensitivities.
        """
        # EPSILON_ALPHA = 1e-7
        EPSILON_ALPHA = 0.0001*self.alpha
        EPSILON_STRESS = 1e-1

        alpha_eps = self.alpha + EPSILON_ALPHA
        eps_ne_rate_eps = self.compute_eps_ne_rate(stress, dt*theta, Temp, alpha=alpha_eps, return_eps_ne=True)

        self.r = self.compute_residue(self.eps_ne_rate, self.alpha, dt)
        r_eps = self.compute_residue(eps_ne_rate_eps, alpha_eps, dt)
        self.h = (r_eps - self.r) / EPSILON_ALPHA
        Q = (eps_ne_rate_eps - self.eps_ne_rate) / EPSILON_ALPHA[:,None,None]
        B = (self.r / self.h)[:,None,None] * Q

        self.P = to.zeros_like(stress)
        stress_eps = stress.clone()
        for i, j in [(0,0), (1,1), (2,2), (0,1), (0,2), (1,2)]:
        	stress_eps[:,i,j] += EPSILON_STRESS
        	eps_ne_rate_eps = self.compute_eps_ne_rate(stress_eps, dt*theta, Temp, return_eps_ne=True)
        	r_eps = self.compute_residue(eps_ne_rate_eps, self.alpha, dt)
        	self.P[:,i,j] = (r_eps - self.r) / EPSILON_STRESS
        	self.P[:,j,i] = self.P[:,i,j]
        	stress_eps[:,i,j] -= EPSILON_STRESS

        H = self.compute_H(Q, self.P)
        H_over_h = H/self.h[:,None,None]

        return B, H_over_h


    def compute_H(self, Q: to.Tensor, P: to.Tensor) -> to.Tensor:
        """
        Build the 6×6 matrix `H` from tensors `Q` and `P`.

        Parameters
        ----------
        Q : torch.Tensor
            Sensitivity of rate to `alpha`, shape (N, 3, 3).
        P : torch.Tensor
            Sensitivity of residue to stress, shape (N, 3, 3).

        Returns
        -------
        torch.Tensor
            `H` with shape (N, 6, 6) in tensorial Voigt ordering.
        """
        n_elems, _, _ = P.shape
        H = to.zeros((n_elems, 6, 6), dtype=to.float64)
        H[:,0,0] = Q[:,0,0]*P[:,0,0]
        H[:,0,1] = Q[:,0,0]*P[:,1,1]
        H[:,0,2] = Q[:,0,0]*P[:,2,2]
        H[:,0,3] = 2*Q[:,0,0]*P[:,0,1]
        H[:,0,4] = 2*Q[:,0,0]*P[:,0,2]
        H[:,0,5] = 2*Q[:,0,0]*P[:,1,2]

        H[:,1,0] = Q[:,1,1]*P[:,0,0]
        H[:,1,1] = Q[:,1,1]*P[:,1,1]
        H[:,1,2] = Q[:,1,1]*P[:,2,2]
        H[:,1,3] = 2*Q[:,1,1]*P[:,0,1]
        H[:,1,4] = 2*Q[:,1,1]*P[:,0,2]
        H[:,1,5] = 2*Q[:,1,1]*P[:,1,2]

        H[:,2,0] = Q[:,2,2]*P[:,0,0]
        H[:,2,1] = Q[:,2,2]*P[:,1,1]
        H[:,2,2] = Q[:,2,2]*P[:,2,2]
        H[:,2,3] = 2*Q[:,2,2]*P[:,0,1]
        H[:,2,4] = 2*Q[:,2,2]*P[:,0,2]
        H[:,2,5] = 2*Q[:,2,2]*P[:,1,2]

        H[:,3,0] = Q[:,0,1]*P[:,0,0]
        H[:,3,1] = Q[:,0,1]*P[:,1,1]
        H[:,3,2] = Q[:,0,1]*P[:,2,2]
        H[:,3,3] = 2*Q[:,0,1]*P[:,0,1]
        H[:,3,4] = 2*Q[:,0,1]*P[:,0,2]
        H[:,3,5] = 2*Q[:,0,1]*P[:,1,2]

        H[:,4,0] = Q[:,0,2]*P[:,0,0]
        H[:,4,1] = Q[:,0,2]*P[:,1,1]
        H[:,4,2] = Q[:,0,2]*P[:,2,2]
        H[:,4,3] = 2*Q[:,0,2]*P[:,0,1]
        H[:,4,4] = 2*Q[:,0,2]*P[:,0,2]
        H[:,4,5] = 2*Q[:,0,2]*P[:,1,2]

        H[:,5,0] = Q[:,1,2]*P[:,0,0]
        H[:,5,1] = Q[:,1,2]*P[:,1,1]
        H[:,5,2] = Q[:,1,2]*P[:,2,2]
        H[:,5,3] = 2*Q[:,1,2]*P[:,0,1]
        H[:,5,4] = 2*Q[:,1,2]*P[:,0,2]
        H[:,5,5] = 2*Q[:,1,2]*P[:,1,2]
        return H

















