"""
It defines classes for all types for elements that can be included into
the constitutive model, such as Spring, DislocationCreep, etc.
"""
# Copyright 2024 The safeincave community.
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

import torch as to
import numpy as np
from Utils import dotdot2, MPa, compute_C

class Spring():
	"""
    This class implements the necessary data and methods for a linear spring element.

    Parameters
    ----------
    props : dict
        Dictionary containing the material properties of a linear elastic spring.
    """
	def __init__(self, props):
		self._E = props["E"]
		self._nu = props["nu"]
		try:
			self.n_elems = self._E.shape[0]
		except:
			self.n_elems = 1
		self._eps_e = to.tensor((self.n_elems, 3, 3), dtype=to.float64)

	def initialize(self):
		"""
		Initialize spring data associated to the linear spring, i.e., :math:`\\mathbb{C}_0` and 
		:math:`\\mathbb{C}_0^{-1}` for each element of the grid.

		.. note::

			The quantity :math:`\\mathbb{C}_0` is a 4th-order tensor, but due to Voigt notation it is represented as 6x6 matrix for each grid element.

		"""
		self._C0 = compute_C(self.n_elems, self._nu, self._E)
		self._C0_inv = to.linalg.inv(self._C0)

	def compute_eps_e(self, stress):
		"""
		This function computes the elastic strain as :math:`\\varepsilon_e = \\mathbb{C}_0^{-1} : \\pmb{\\sigma}`.
		"""
		# :math:`\varepsilon_e = \mathbb{C}_0^{-1} : \pmb{\sigma}`
		self._eps_e = dotdot2(self._C0_inv, stress)

	@property
	def E(self):
		"""
		list : This list contains the Young's modulus for each element of the grid.
		"""
		return self._E

	@property
	def nu(self):
		"""
		list : This list contains the Poisson's ratio for each element of the grid.
		"""
		return self._nu

	@property
	def C0(self):
		"""
		torch.Tensor : This is a (nelems, 6, 6) tensor storing the :math:`\\mathbb{C}_0` matrix for each element of the grid.
		"""
		return self._C0
	
	@property
	def C0_inv(self):
		"""
		torch.Tensor : This is a (nelems, 6, 6) tensor storing :math:`\\mathbb{C}_0^{-1}` matrix for each element of the grid.
		"""
		return self._C0_inv

	@property
	def eps_e(self):
		"""
		torch.Tensor : This is a (nelems, 3, 3) tensor storing the elastic strain tensor for all grid elements.
		"""
		return self._eps_e




	


class Viscoelastic():
	"""
    This class implements the necessary data and methods for a linear Kelvin-Voigt element (viscoelastic).

    Parameters
    ----------
    props : dict
        Dictionary containing the material properties of a Kelvin-Voigt element. It includes Young's modulus 
        and Poisson's ratio for the spring and a viscosity for the dashpot.
    """
	def __init__(self, props):

		self._eta = props["eta"]
		self._E = props["E"]
		self._nu = props["nu"]
		try:
			self.n_elems = self._E.shape[0]
		except:
			self.n_elems = 1
		self._eps_ve_rate = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._eps_ve_rate_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._eps_ve_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._eps_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._eps_bar = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self._B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._G = to.zeros((self.n_elems, 6, 6), dtype=to.float64)

		# Assemble C1 tensor (n_elems, 6, 6)
		self._C1 = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		a0 = self._E/((1 + self._nu)*(1 - 2*self._nu))
		self._C1[:,0,0] = a0*(1 - self._nu)
		self._C1[:,1,1] = a0*(1 - self._nu)
		self._C1[:,2,2] = a0*(1 - self._nu)
		self._C1[:,3,3] = a0*(1 - 2*self._nu)
		self._C1[:,4,4] = a0*(1 - 2*self._nu)
		self._C1[:,5,5] = a0*(1 - 2*self._nu)
		self._C1[:,0,1] = self._C1[:,1,0] = self._C1[:,0,2] = self._C1[:,2,0] = self._C1[:,2,1] = self._C1[:,1,2] = a0*self.nu

	def compute_G_B(self, stress_vec, phi2):
		"""
		Compute matrices :math:`\\mathbb{G}_{ve}` and :math:`\\mathbf{B}_{ve}` for the Kelvin-Voigt element.

		.. note::
		
			Since the viscoelastic element is linear, the value of :math:`\\mathbb{G}_{ve}` is constant,
			and it is given by :math:`\\mathbb{G}_{ve} = \\left( \\eta \\mathbf{I} + \\phi_2 \\mathbb{C}_1 \\right)^{-1}`.

		.. note::
		
			For a viscoelastic element, since it does not depend on internal parameters,
			we must have :math:`\\mathbf{B}_{ve} = 0`.

		Parameters
		----------
		stress_vec : torch.Tensor
			This is a (nelems, 3, 3) storing the stress tensor for each grid element.

		phi2 : float
			Quantity associated to the time integration, :math:`\\phi_2 = \\Delta t (1 - \\theta)`.
		"""
		I = to.eye(6, dtype=to.float64).unsqueeze(0).repeat(self.n_elems, 1, 1)
		self._G = to.linalg.inv(self._eta[:,None,None]*I + phi2*self._C1)

	def increment_internal_variables(self, *args):
		"""
		The viscoelastic element does not depend on internal variables, so this function does nothing.
		"""
		pass

	def update_internal_variables(self, *args):
		"""
		The viscoelastic element does not depend on internal variables, so this function does nothing.
		"""
		pass

	def compute_eps_bar(self, phi1, phi2):
		"""
		Computes :math:`\\bar{\\pmb{\\varepsilon}}_{ve}^k = \\pmb{\\varepsilon}_{ve}^t + \\phi_1 \\dot{\\pmb{\\varepsilon}}_{ve}^t + \\phi_2 \\dot{\\pmb{\\varepsilon}}_{ve}^k`.
		
		Parameters
		----------
		phi1 : float
			Quantity associated to time integration, :math:`\\phi_1 = \\Delta t \\theta`. 

		phi2 : float
			Quantity associated to time integration, :math:`\\phi_2 = \\Delta t (1 - \\theta)`. 

		"""
		self._eps_bar = self._eps_ve_old + phi1*self._eps_ve_rate_old + phi2*self._eps_ve_rate

	def compute_eps_ve(self, stress, stress_k, phi2):
		"""
		Computes the viscoelastic strain as:

		.. math::

			\\pmb{\\varepsilon}_{ve} = \\bar{\\pmb{\\varepsilon}}_{ve}^k + \\phi_2 \\left[ \\mathbb{G}_{ve} : \\left( \\pmb{\\sigma}^{k+1} - \\pmb{\\sigma}^k \\right) - \\mathbf{B}_{ve} \\right]
		
		Parameters
		----------
		stress : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress of the **current** iteration level 
			(i.e. :math:`\\pmb{\\sigma}^{k+1}`).

		stress_k : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress of the **previous** iteration level 
			(i.e. :math:`\\pmb{\\sigma}^{k}`).

		phi2 : float
			Quantity associated to time integration, :math:`\\phi_2 = \\Delta t (1 - \\theta)`.
		"""
		self._eps_ve = self._eps_bar + phi2*dotdot2(self._G, stress - stress_k) - phi2*self._B

	def update_eps_ve_old(self):
		"""
		Updates viscoelastic strain at time :math:`t` with the value of time :math:`t + \\Delta t`,
		that is, :math:`\\pmb{\\varepsilon}_{ve}^t \\leftarrow \\pmb{\\varepsilon}_{ve}^{t+\\Delta t}`
		"""
		self._eps_ve_old = self._eps_ve.clone()

	def update_eps_ve_rate_old(self):
		"""
		Updates viscoelastic strain rate at time :math:`t` with the value of time :math:`t + \\Delta t`,
		that is, :math:`\\dot{\\pmb{\\varepsilon}}_{ve}^t \\leftarrow \\dot{\\pmb{\\varepsilon}}_{ve}^{t+\\Delta t}`
		"""
		self._eps_ve_rate_old = self._eps_ve_rate.clone()

	def compute_eps_ve_rate(self, stress_vec, phi1, return_eps_ve=False):
		"""
		Computes the viscoelastic strain rate as

		.. math::
			\\dot{\\pmb{\\varepsilon}}_{ve} = \\mathbb{G}_{ve} : \\left[ \\pmb{\\sigma} - \\mathbb{C}_1 : \\left( \\pmb{\\varepsilon}_{ve}^t + \\phi_1 \\dot{\\pmb{\\varepsilon}}_{ve}^t \\right) \\right]

		"""
		eps_ve_rate = dotdot2(self._G, stress_vec - dotdot2(self._C1, self.eps_ve_old + phi1*self._eps_ve_rate_old))
		if return_eps_ve:
			return eps_ve_rate.clone()
		else:
			self._eps_ve_rate = eps_ve_rate.clone()

	@property
	def eta(self):
		"""
		list : This list contains the viscosities associated to the Kelvin-Voigt dashpot for each element of the grid.
		"""
		return self._eta

	@property
	def E(self):
		"""
		list : This list contains the Young's modulus for each element of the grid.
		"""
		return self._E

	@property
	def nu(self):
		"""
		list : This list contains the Poisson's ratio for each element of the grid.
		"""
		return self._nu

	@property
	def eps_ve_rate(self):
		"""
		torch.Tensor : A (n_elems, 3, 3) tensor storing :math:`\\dot{\\pmb{\\varepsilon}}_{ve}` for all grid elements.
		"""
		return self._eps_ve_rate

	@property
	def eps_ve_rate_old(self):
		"""
		torch.Tensor : A (n_elems, 3, 3) tensor storing :math:`\\dot{\\pmb{\\varepsilon}}_{ve}^t` (i.e. at the 
		**previous** time level) for all grid elements.
		"""
		return self._eps_ve_rate_old

	@property
	def eps_ve_old(self):
		"""
		torch.Tensor : A (n_elems, 3, 3) tensor storing :math:`\\pmb{\\varepsilon}_{ve}^t` (i.e. at the 
		**previous** time level) for all grid elements.
		"""
		return self._eps_ve_old

	@property
	def eps_ve(self):
		"""
		torch.Tensor : A (nelems, 3, 3) tensor storing :math:`\\pmb{\\varepsilon}_{ve}` (i.e. at the 
		**current** time level) for all grid elements.
		"""
		return self._eps_ve

	@property
	def eps_bar(self):
		"""
		torch.Tensor : A (nelems, 3, 3) tensor storing :math:`\\bar{\\pmb{\\varepsilon}}_{ve}` for all grid elements.
		"""
		return self._eps_bar

	@property
	def B(self):
		"""
		torch.Tensor : A (nelems, 3, 3) tensor storing matrix :math:`\\textbf{B}_{ve}` for all grid elements.
		"""
		return self._B

	@property
	def G(self):
		"""
		torch.Tensor : A (nelems, 6, 6) tensor storing :math:`\\mathbb{G}_{ve}` for all grid elements.
		"""
		return self._G

	@property
	def C1(self):
		"""
		torch.Tensor : A (nelems, 6, 6) tensor storing :math:`\\mathbb{C}_{1}` for all grid elements.
		"""
		return self._C1
	



class DislocationCreep():
	"""
    This class implements the necessary data and methods for the dislocation creep element (non-linear dashpot).

    Parameters
    ----------
    props : dict
        Dictionary containing the material properties of a dislocation creep element. It includes
        the following parameters:

        - *A*: [:math:`\\text{Pa}^{-1} \\text{s}^{-1}`].
        - *n*: [:math:`-`].
        - *Q*: Activation energy [:math:`\\text{J}/\\text{mol}`].
        - *R*: Universal gas constant [8.32 :math:`\\text{JK}^{-1}\\text{mol}^{-1}`].
        - *T*: Temperature [:math:`\\text{K}`].
    """
	def __init__(self, props):
		R = 8.32
		self._Aexp = props["A"]*to.exp(-props["Q"]/R/props["T"])
		self._n = props["n"]
		try:
			self.n_elems = self._Aexp.shape[0]
		except:
			self.n_elems = 1
		self._eps_ie_rate = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._eps_ie_rate_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._eps_ie_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._eps_bar = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self._B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._G = to.zeros((self.n_elems, 6, 6), dtype=to.float64)

	def increment_internal_variables(self, *args):
		"""
		The disclocation creep element does not depend on internal variables, so this function does nothing.
		"""
		pass

	def update_internal_variables(self, *args):
		"""
		The disclocation creep element does not depend on internal variables, so this function does nothing.
		"""
		pass

	def compute_eps_ie(self, stress, stress_k, phi2):
		"""
		Computes the dislocation creep strain as:

		.. math::

			\\pmb{\\varepsilon}_{cr} = \\bar{\\pmb{\\varepsilon}}_{cr}^k + \\phi_2 \\left[ \\mathbb{G}_{cr} : \\left( \\pmb{\\sigma}^{k+1} - \\pmb{\\sigma}^k \\right) - \\mathbf{B}_{cr} \\right]
		
		Parameters
		----------
		stress : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress of the **current** iteration level 
			(i.e. :math:`\\pmb{\\sigma}^{k+1}`).

		stress_k : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress of the **previous** iteration level 
			(i.e. :math:`\\pmb{\\sigma}^{k}`).

		phi2 : float
			Quantity associated to time integration, :math:`\\phi_2 = \\Delta t (1 - \\theta)`.
		"""
		self._eps_ie = self._eps_bar + phi2*dotdot2(self._G, stress - stress_k) - phi2*self._B

	def update_eps_ie_old(self):
		"""
		Updates dislocation creep strain at time :math:`t` with the value of time :math:`t + \\Delta t`,
		that is, :math:`\\pmb{\\varepsilon}_{cr}^t \\leftarrow \\pmb{\\varepsilon}_{cr}^{t+\\Delta t}`
		"""
		self._eps_ie_old = self._eps_ie.clone()

	def update_eps_ie_rate_old(self):
		"""
		Updates dislocation creep strain rate at time :math:`t` with the value of time :math:`t + \\Delta t`,
		that is, :math:`\\dot{\\pmb{\\varepsilon}}_{cr}^t \\leftarrow \\dot{\\pmb{\\varepsilon}}_{cr}^{t+\\Delta t}`
		"""
		self._eps_ie_rate_old = self._eps_ie_rate.clone()

	def compute_eps_bar(self, phi1, phi2):
		"""
		Computes :math:`\\bar{\\pmb{\\varepsilon}}_{cr}^k = \\pmb{\\varepsilon}_{cr}^t + \\phi_1 \\dot{\\pmb{\\varepsilon}}_{cr}^t + \\phi_2 \\dot{\\pmb{\\varepsilon}}_{cr}^k`.
		
		Parameters
		----------
		phi1 : float
			Quantity associated to time integration, :math:`\\phi_1 = \\Delta t \\theta`. 

		phi2 : float
			Quantity associated to time integration, :math:`\\phi_2 = \\Delta t (1 - \\theta)`. 

		"""
		self._eps_bar = self._eps_ie_old + phi1*self._eps_ie_rate_old + phi2*self._eps_ie_rate

	def compute_eps_ie_rate(self, stress_vec, return_eps_ie=False, *args):
		"""
		Computes the dislocation creep strain rate as

		.. math::
			\\dot{\\pmb{\\varepsilon}}_{cr} = A \\exp \\left( -\\frac{Q}{RT} \\right) q^{n-1} \\mathbf{s}

		where :math:`\\mathbf{s} = \\pmb{\\sigma} - \\frac{1}{3} \\text{tr}(\\pmb{\\sigma}) \\mathbf{I}` is
		the deviatoric stress and :math:`q = \\sqrt{\\frac{3}{2} \\mathbf{s} : \\mathbf{s}}` denotes the
		von Mises stress.

		Parameters
		----------
		stress_vec : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress tensor for all grid elements.

		return_eps_ie : bool
			If *False* it stores :math:`\\dot{\\pmb{\\varepsilon}}_{ie}` to *eps_ie_rate*, otherwise the function returns as a value.
			This is used when computing :math:`\\frac{\\partial \\dot{\\pmb{\\varepsilon}}_{ie}}{\\partial \\pmb{\\sigma}}`
			by finite differences.

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

		A_bar = self._Aexp*q_vm**(self._n - 1)
		if return_eps_ie:
			return A_bar[:,None,None]*dev
		else:
			self._eps_ie_rate = A_bar[:,None,None]*dev

	# def compute_E(self, stress_vec):
	# 	"""
	# 	This function computes :math:`\\frac{\\partial \\dot{\\pmb{\\varepsilon}}_{ie}}{\\partial \\pmb{\\sigma}}`
	# 	by finite differences.

	# 	Parameters
	# 	----------
	# 	stress_vec : torch.Tensor
	# 		A (nelems, 3, 3) pytorch tensor storing the stress tensor for all grid elements.

	# 	"""
	# 	EPSILON = 1e-2
	# 	self._E = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
	# 	stress_eps = stress_vec.clone()
	# 	c1 = 1.0
	# 	c2 = 2.0
	# 	magic_indexes = [(0,0,0,c1), (1,1,1,c1), (2,2,2,c1), (0,1,3,c2), (0,2,4,c2), (1,2,5,c2)]
	# 	for i, j, k, phi in magic_indexes:
	# 		stress_eps[:,i,j] += EPSILON
	# 		eps_ie_rate_eps = self.compute_eps_ie_rate(stress_eps, return_eps_ie=True)
	# 		self._E[:,:,k] = phi*(eps_ie_rate_eps[:,[0,1,2,0,0,1],[0,1,2,1,2,2]] - self._eps_ie_rate[:,[0,1,2,0,0,1],[0,1,2,1,2,2]]) / EPSILON
	# 		stress_eps[:,i,j] -= EPSILON

	def compute_E(self, stress_vec):
		"""
		This function computes :math:`\\frac{\\partial \\dot{\\pmb{\\varepsilon}}_{ie}}{\\partial \\pmb{\\sigma}}`
		by finite differences.

		Parameters
		----------
		stress_vec : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress tensor for all grid elements.

		"""
		EPSILON = 1e-2
		self._E = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		stress_eps = stress_vec.clone()
		c1 = 1.0
		c2 = 1.0
		magic_indexes = [(0,0,0,c1), (1,1,1,c1), (2,2,2,c1), (0,1,3,c2), (0,2,4,c2), (1,2,5,c2)]
		for i, j, k, phi in magic_indexes:
			stress_eps[:,i,j] += EPSILON
			eps_ie_rate_eps = self.compute_eps_ie_rate(stress_eps, return_eps_ie=True)

			if k <= 2:
				self._E[:,:,k] += phi*(eps_ie_rate_eps[:,[0,1,2,0,0,1],[0,1,2,1,2,2]] - self._eps_ie_rate[:,[0,1,2,0,0,1],[0,1,2,1,2,2]]) / EPSILON

			else:
				self._E[:,:,k] += (eps_ie_rate_eps[:,[0,1,2,0,0,1],[0,1,2,1,2,2]] - self._eps_ie_rate[:,[0,1,2,0,0,1],[0,1,2,1,2,2]]) / EPSILON
				self._E[:,:,k] += (eps_ie_rate_eps[:,[0,1,2,1,2,2], [0,1,2,0,0,1]] - self._eps_ie_rate[:,[0,1,2,1,2,2], [0,1,2,0,0,1]]) / EPSILON

			stress_eps[:,i,j] -= EPSILON

	def compute_G_B(self, stress_vec, *args):
		"""
		This function computes matrices :math:`\\mathbb{G}_{ie}` and :math:`\\mathbf{B}_{ie}`.

		.. note::

			In this case, :math:`\\mathbf{B}_{ie} = 0` since the dislocation creep element does not depend
			on internal parameters.

		Parameters
		----------
		stress_vec : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress tensor for all grid elements.
		
		"""
		self.compute_E(stress_vec)
		self._G = self._E.clone()

	@property
	def Aexp(self):
		"""
		list : List containing the term :math:`A_\\text{exp} = \\exp \\left( -Q/RT \\right)`
		for all grid elements.
		"""
		return self._Aexp

	@property
	def n(self):
		"""
		list : List containing the material property *n* (:math:`-`) for all grid elements.
		"""
		return self._n

	@property
	def eps_ie_rate(self):
		"""
		torch.Tensor : A (n_elems, 3, 3) tensor storing :math:`\\dot{\\pmb{\\varepsilon}}_{cr}` for all grid elements.
		"""
		return self._eps_ie_rate

	@property
	def eps_ie_rate_old(self):
		"""
		torch.Tensor : A (n_elems, 3, 3) tensor storing :math:`\\dot{\\pmb{\\varepsilon}}_{cr}^t` (i.e. at the 
		**previous** time level) for all grid elements.
		"""
		return self._eps_ie_rate_old

	@property
	def eps_ie_old(self):
		"""
		torch.Tensor : A (n_elems, 3, 3) tensor storing :math:`\\pmb{\\varepsilon}_{cr}^t` (i.e. at the 
		**previous** time level) for all grid elements.
		"""
		return self._eps_ie_old

	@property
	def eps_ie(self):
		"""
		torch.Tensor : A (nelems, 3, 3) tensor storing :math:`\\pmb{\\varepsilon}_{ie}` (i.e. at the 
		**current** time level) for all grid elements.
		"""
		return self._eps_ie
	

	@property
	def eps_bar(self):
		"""
		torch.Tensor : A (nelems, 3, 3) tensor storing :math:`\\bar{\\pmb{\\varepsilon}}_{cr}` for all grid elements.
		"""
		return self._eps_bar

	@property
	def B(self):
		"""
		torch.Tensor : A (nelems, 3, 3) tensor storing matrix :math:`\\textbf{B}_{cr}` for all grid elements.
		"""
		return self._B

	@property
	def G(self):
		"""
		torch.Tensor : A (nelems, 6, 6) tensor storing :math:`\\mathbb{G}_{cr}` for all grid elements.
		"""
		return self._G



class ViscoplasticDesai():
	"""
    This class implements the necessary data and methods for the viscoplastic model proposed by Desai (1987).

    Parameters
    ----------
    props : dict
        Dictionary containing the material properties of viscoplastic model. It includes the 
        following parameters:

        - :math:`\\mu_1`: [:math:`\\text{s}^{-1}`].
        - :math:`N_1`: [:math:`-`].
        - :math:`a_1`: [:math:`\\text{MPa}^{2-n}`].
        - :math:`\\eta`: [:math:`-`].
        - :math:`n`: [:math:`-`].
        - :math:`\\beta_1`: [:math:`\\text{MPa}^{-1}`].
        - :math:`\\beta`: [:math:`-`].
        - :math:`m`: [:math:`-`].
        - :math:`\\gamma`: [:math:`-`].
        - :math:`\\sigma_t`: Tensile strength [:math:`\\text{MPa}`].
        - :math:`\\alpha_0`: Initial hardening parameter [:math:`-`].

    """
	def __init__(self, props):
		self.F_0 = 1.0
		self._mu_1 = props["mu_1"]
		self._N_1 = props["N_1"]
		self._a_1 = props["a_1"]
		self._eta = props["eta"]
		self._n = props["n"]
		self._beta_1 = props["beta_1"]
		self._beta = props["beta"]
		self._m = props["m"]
		self._gamma = props["gamma"]
		self._sigma_t = props["sigma_t"]
		self._alpha_0 = props["alpha_0"]
		try:
			self._alpha = self._alpha_0.clone()
		except:
			self._alpha = self._alpha_0

		try:
			self.n_elems = self._alpha_0.shape[0]
		except:
			self.n_elems = 1

		self._Fvp = to.zeros(self.n_elems, dtype=to.float64)
		self._qsi = to.zeros(self.n_elems, dtype=to.float64)
		self._qsi_old = to.zeros(self.n_elems, dtype=to.float64)
		self._r = to.zeros(self.n_elems, dtype=to.float64)
		self._eps_ie_rate = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._eps_ie_rate_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._eps_ie_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self._eps_bar = to.zeros((self.n_elems, 3, 3), dtype=to.float64)


	def compute_G_B(self, stress_vec, dt):
		"""
		This function computes matrices :math:`\\mathbb{G}_{vp}` and :math:`\\mathbf{B}_{vp}`.
		Since the viscoplastic model depends on a internal parameter (:math:`\\alpha`), these
		matrices are given by

		.. math::

			\\mathbb{G}_{i} = \\frac{\\partial \\dot{\\pmb{\\varepsilon}}_{i}}{\\partial \\pmb{\\sigma}} - \\frac{1}{h_i} \\frac{\\partial \\dot{\\pmb{\\varepsilon}}_{i}}{\\partial \\alpha_i} \\frac{\\partial r_i}{\\partial \\pmb{\\sigma}}
			\\quad \\text{and} \\quad
			\\mathbf{B}_{i} = \\frac{r_i}{h_i} \\frac{\\partial \\dot{\\pmb{\\varepsilon}}_{i}}{\\partial \\alpha_i}

		where :math:`i = vp`.

		.. note::

			Remember that :math:`\\mathbb{G}_{ie}` is a 4th-order tensor represented as matrix
			by using Voigt notation.

		Parameters
		----------
		stress_vec : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress tensor for all grid elements.

		dt : float
			Time step size.
		
		"""
		# EPSILON_ALPHA = 1e-7
		EPSILON_ALPHA = 0.0001*self._alpha
		EPSILON_STRESS = 1e-1

		alpha_eps = self._alpha + EPSILON_ALPHA
		eps_ie_rate_eps = self.compute_eps_ie_rate(stress_vec, alpha=alpha_eps, return_eps_ie=True)

		self._r = self.compute_residue(self._eps_ie_rate, self._alpha, dt)
		r_eps = self.compute_residue(eps_ie_rate_eps, alpha_eps, dt)
		self._h = (r_eps - self._r) / EPSILON_ALPHA
		self._Q = (eps_ie_rate_eps - self._eps_ie_rate) / EPSILON_ALPHA[:,None,None]
		self._B = (self._r / self._h)[:,None,None] * self._Q

		self._P = to.zeros_like(stress_vec)
		stress_eps = stress_vec.clone()
		for i, j in [(0,0), (1,1), (2,2), (0,1), (0,2), (1,2)]:
			stress_eps[:,i,j] += EPSILON_STRESS
			eps_ie_rate_eps = self.compute_eps_ie_rate(stress_eps, return_eps_ie=True)
			r_eps = self.compute_residue(eps_ie_rate_eps, self._alpha, dt)
			self._P[:,i,j] = (r_eps - self._r) / EPSILON_STRESS
			self._P[:,j,i] = self._P[:,i,j]
			stress_eps[:,i,j] -= EPSILON_STRESS

		self._H = self.compute_H(self._Q, self._P)

		self.compute_E(stress_vec)
		self._G = self._E - self._H/self._h[:,None,None]

	# def compute_E(self, stress_vec):
	# 	"""
	# 	This function computes :math:`\\frac{\\partial \\dot{\\pmb{\\varepsilon}}_{vp}}{\\partial \\pmb{\\sigma}}`
	# 	by finite differences.

	# 	Parameters
	# 	----------
	# 	stress_vec : torch.Tensor
	# 		A (nelems, 3, 3) pytorch tensor storing the stress tensor for all grid elements.

	# 	"""
	# 	EPSILON_STRESS = 1e-1
	# 	self._E = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
	# 	stress_eps = stress_vec.clone()
	# 	c1 = 1.0
	# 	c2 = 2.0
	# 	magic_indexes = [(0,0,0,c1), (1,1,1,c1), (2,2,2,c1), (0,1,3,c2), (0,2,4,c2), (1,2,5,c2)]
	# 	for i, j, k, phi in magic_indexes:
	# 		stress_eps[:,i,j] += EPSILON_STRESS
	# 		eps_ie_rate_eps = self.compute_eps_ie_rate(stress_eps, return_eps_ie=True)
	# 		self._E[:,:,k] = phi*(eps_ie_rate_eps[:,[0,1,2,0,0,1],[0,1,2,1,2,2]] - self._eps_ie_rate[:,[0,1,2,0,0,1],[0,1,2,1,2,2]]) / EPSILON_STRESS
	# 		stress_eps[:,i,j] -= EPSILON_STRESS

	def compute_E(self, stress_vec):
		"""
		This function computes :math:`\\frac{\\partial \\dot{\\pmb{\\varepsilon}}_{ie}}{\\partial \\pmb{\\sigma}}`
		by finite differences.

		Parameters
		----------
		stress_vec : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress tensor for all grid elements.

		"""
		EPSILON = 1e-2
		self._E = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		stress_eps = stress_vec.clone()
		c1 = 1.0
		c2 = 1.0
		magic_indexes = [(0,0,0,c1), (1,1,1,c1), (2,2,2,c1), (0,1,3,c2), (0,2,4,c2), (1,2,5,c2)]
		for i, j, k, phi in magic_indexes:
			stress_eps[:,i,j] += EPSILON
			eps_ie_rate_eps = self.compute_eps_ie_rate(stress_eps, return_eps_ie=True)

			if k <= 2:
				self._E[:,:,k] += phi*(eps_ie_rate_eps[:,[0,1,2,0,0,1],[0,1,2,1,2,2]] - self._eps_ie_rate[:,[0,1,2,0,0,1],[0,1,2,1,2,2]]) / EPSILON

			else:
				self._E[:,:,k] += (eps_ie_rate_eps[:,[0,1,2,0,0,1],[0,1,2,1,2,2]] - self._eps_ie_rate[:,[0,1,2,0,0,1],[0,1,2,1,2,2]]) / EPSILON
				self._E[:,:,k] += (eps_ie_rate_eps[:,[0,1,2,1,2,2], [0,1,2,0,0,1]] - self._eps_ie_rate[:,[0,1,2,1,2,2], [0,1,2,0,0,1]]) / EPSILON

			stress_eps[:,i,j] -= EPSILON

	def compute_residue(self, eps_rate, alpha, dt):
		"""
		Computes the residue of the hardening rule. In this case,

		.. math::

			r^k_{vp} = \\alpha^k - a_1 \\left[ \\left( \\frac{a_1}{\\alpha_0} \\right)^{1/\\eta} + \\xi^k \\right]^{-\\eta}, \\quad \\text{where} \\quad \\xi^k = \\int_{t_0}^t \\sqrt{ \\dot{\\pmb{\\varepsilon}}^k_{vp} : \\dot{\\pmb{\\varepsilon}}^k_{vp} } \\mathrm{dt}
		
		Parameters
		----------
		eps_rate : torch.Tensor
			A (nelems, 3, 3) tensor storing :math:`\\dot{\\varepsilon}_{vp}` for all grid elements.

		alpha : list
			List containing the hardening parameter values for all grid elements.

		dt : float
			Time step size.
		"""
		self._qsi = self._qsi_old + to.sum(eps_rate**2, axis=(-2, -1))**0.5*dt
		return alpha - self._a_1 / (((self._a_1/self._alpha_0)**(1/self._eta) + self._qsi)**self._eta)

	def compute_eps_bar(self, phi1, phi2):
		"""
		Computes :math:`\\bar{\\pmb{\\varepsilon}}_{cr}^k = \\pmb{\\varepsilon}_{cr}^t + \\phi_1 \\dot{\\pmb{\\varepsilon}}_{cr}^t + \\phi_2 \\dot{\\pmb{\\varepsilon}}_{cr}^k`.
		
		Parameters
		----------
		phi1 : float
			Quantity associated to time integration, :math:`\\phi_1 = \\Delta t \\theta`. 

		phi2 : float
			Quantity associated to time integration, :math:`\\phi_2 = \\Delta t (1 - \\theta)`. 

		"""
		self._eps_bar = self._eps_ie_old + phi1*self._eps_ie_rate_old + phi2*self._eps_ie_rate

	def compute_eps_ie(self, stress, stress_k, phi2):
		"""
		Computes the dislocation creep strain as:

		.. math::

			\\pmb{\\varepsilon}_{vp} = \\bar{\\pmb{\\varepsilon}}_{vp}^k + \\phi_2 \\left[ \\mathbb{G}_{vp} : \\left( \\pmb{\\sigma}^{k+1} - \\pmb{\\sigma}^k \\right) - \\mathbf{B}_{vp} \\right]
		
		Parameters
		----------
		stress : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress of the **current** iteration level 
			(i.e. :math:`\\pmb{\\sigma}^{k+1}`).

		stress_k : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress of the **previous** iteration level 
			(i.e. :math:`\\pmb{\\sigma}^{k}`).

		phi2 : float
			Quantity associated to time integration, :math:`\\phi_2 = \\Delta t (1 - \\theta)`.
		"""
		self._eps_ie = self._eps_bar + phi2*dotdot2(self.G, stress - stress_k) - phi2*self._B

	def update_eps_ie_old(self):
		"""
		Updates viscoplastic strain at time :math:`t` with the value of time :math:`t + \\Delta t`,
		that is, :math:`\\pmb{\\varepsilon}_{vp}^t \\leftarrow \\pmb{\\varepsilon}_{vp}^{t+\\Delta t}`
		"""
		self._eps_ie_old = self.eps_ie.clone()

	def update_eps_ie_rate_old(self):
		"""
		Updates viscoplastic strain rate at time :math:`t` with the value of time :math:`t + \\Delta t`,
		that is, :math:`\\dot{\\pmb{\\varepsilon}}_{vp}^t \\leftarrow \\dot{\\pmb{\\varepsilon}}_{vp}^{t+\\Delta t}`
		"""
		self._eps_ie_rate_old = self._eps_ie_rate.clone()

	def update_internal_variables(self):
		"""
		Updates the previous value of the accumulated viscoplastic strain with the current one, that is
		:math:`\\xi^t \\leftarrow \\xi^{t+\\Delta t}`.
		"""
		self._qsi_old = self._qsi.clone()

	def increment_internal_variables(self, stress, stress_k, dt):
		"""
		Increment the hardening parameter by the following expression

		.. math::

			\\alpha^{k+1}_i \\leftarrow \\alpha^k_i + \\delta \\alpha_i
			\\quad \\text{where} \\quad
			\\delta \\alpha_{i} = - \\frac{1}{h_i} \\left( r^k_{i} + \\frac{\\partial r_i^k}{\\partial \\pmb{\\sigma}} : \\delta \\pmb{\\sigma} \\right).
		
		where :math:`i=vp`.
		"""
		delta_alpha = -(self._r + to.einsum('bij,bij->b', self._P, stress - stress_k))/self._h
		self._alpha += delta_alpha

	def compute_stress_invariants(self, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz):
		"""
		This function compute the following stress invariants

		.. math::

			I_1 = \\sigma_{xx} + \\sigma_{yy} + \\sigma_{zz}

		.. math::

			I_2 = \\sigma_{xx}\\sigma_{yy} + \\sigma_{yy}\\sigma_{zz} + \\sigma_{xx}\\sigma_{zz} - \\sigma_{xy}^2 - \\sigma_{yz}^2 - \\sigma_{xz}^2

		.. math::

			I_3 = \\sigma_{xx}\\sigma_{yy}\\sigma_{zz} + 2\\sigma_{xy}\\sigma_{yz}\\sigma_{xz} - \\sigma_{zz}\\sigma_{xy}^2 - \\sigma_{xx}\\sigma_{yz}^2 - \\sigma_{yy}\\sigma_{xz}^2

		.. math::

			J_2 = \\frac{1}{3}*I_1^2 - I_2

		.. math::

			J_3 = \\frac{2}{27} I_1^3 - \\frac{1}{3} I_1 I_2 + I_3

		.. math::

			 S_r = -\\frac{J_3\\sqrt{27}}{2 J_2^{1.5}}

		.. math::

			I_1^* = I_1 + \\sigma_{t}

		Parameters
		----------
		s_xx : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{xx}` for all grid elements.

		s_yy : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{yy}` for all grid elements.
			
		s_zz : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{zz}` for all grid elements.
			
		s_xy : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{xy}` for all grid elements.
			
		s_xz : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{xz}` for all grid elements.
			
		s_yz : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{yz}` for all grid elements.

		Returns
		-------
		I1 : torch.Tensor
			A (nelems,) array storing :math:`I_1` for all grid elements.

		I2 : torch.Tensor
			A (nelems,) array storing :math:`I_2` for all grid elements.

		I3 : torch.Tensor
			A (nelems,) array storing :math:`I_3` for all grid elements.

		J2 : torch.Tensor
			A (nelems,) array storing :math:`J_2` for all grid elements.

		J3 : torch.Tensor
			A (nelems,) array storing :math:`J_3` for all grid elements.

		Sr : torch.Tensor
			A (nelems,) array storing :math:`S_r` for all grid elements.

		I1_star : torch.Tensor
			A (nelems,) array storing :math:`I_1^*` for all grid elements.
		"""
		I1 = s_xx + s_yy + s_zz
		I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
		I3 = s_xx*s_yy*s_zz + 2*s_xy*s_yz*s_xz - s_zz*s_xy**2 - s_xx*s_yz**2 - s_yy*s_xz**2
		J2 = (1/3)*I1**2 - I2
		J3 = (2/27)*I1**3 - (1/3)*I1*I2 + I3
		Sr = -(J3*np.sqrt(27))/(2*J2**1.5)
		I1_star = I1 + self._sigma_t
		# print("I1", I1.max())
		# print("I2", I2.max())
		# print("J2", J2.max())
		# print()
		return I1, I2, I3, J2, J3, Sr, I1_star

	def extract_stress_components(self, stress):
		"""
		This is a convenient function to extract the stress components.

		Parameters
		----------
		stress : torch.Tensor
			A (nelems, 3, 3) storing the stress tensor for all grid elements.

		Returns
		-------
		s_xx : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{xx}` for all grid elements.

		s_yy : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{yy}` for all grid elements.
			
		s_zz : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{zz}` for all grid elements.
			
		s_xy : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{xy}` for all grid elements.
			
		s_xz : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{xz}` for all grid elements.
			
		s_yz : torch.Tensor
			A (nelems,) array storing :math:`\\sigma_{yz}` for all grid elements.
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
		Computes the yield function value according to

		.. math::

			F_{vp}(\\pmb{\\sigma}, \\alpha) = J_2 - (-\\alpha I_1^{n} + \\gamma I_1^2) \\left[ \\exp{(\\beta_1 I_1)} - \\beta S_r \\right]^m

		Parameters
		----------
		alpha : torch.Tensor
			A (nelems,) storing the hardening parameters for all grid elements.
		I1 : torch.Tensor
			A (nelems,) storing :math:`I_1` for all grid elements.
		J2 : torch.Tensor
			A (nelems,) storing :math:`J_2` for all grid elements.
		Sr : torch.Tensor
			A (nelems,) storing :math:`S_r` for all grid elements.

		Returns
		-------
		Fvp : torch.Tensor
			A (nelems,) storing the yield function values :math:`F_{vp}` for all grid elements.
		"""
		F1 = (alpha*I1**self._n - self._gamma*I1**2)
		F2 = (to.exp(self._beta_1*I1) - self._beta*Sr)
		# print(F1[0])
		# print("J2", J2.max())
		# print("I1", I1.max())
		# print()
		Fvp = J2 + F1*F2**self._m
		return Fvp

	def compute_initial_hardening(self, stress, Fvp_0=0.0):
		"""
		Computes the initial hardening parameter values such that :math:`F_{vp} = F_{vp,0}`.
		This is done through the following expression:

		.. math::

			\\alpha_0 = \\gamma I_1^{2-n} - I_1^{-n} J_2 \\left[ \\exp\\left( \\beta_1 I_1 \\right) - \\beta S_r \\right]^{-m}

		.. note::

			Specifying :math:`F_{vp} = 0` means that all the elements are on the onset of viscoplastic deformation.

		Parameters
		----------
		stress : torch.Tensor
			A (nelems, 3, 3) storing the stress tensor for all grid elements.
		Fvp_0 : float
			Specified value for the yield function value in the entire domain.
		"""
		s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.extract_stress_components(stress)
		I1, I2, I3, J2, J3, Sr, I1_star = self.compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)
		self._alpha_0 =  self._gamma*I1_star**(2-self._n) + (Fvp_0 - J2)*I1_star**(-self._n)*(to.exp(self._beta_1*I1_star) - self._beta*Sr)**(-self._m)
		self._alpha = self._alpha_0.clone()

	def compute_eps_ie_rate(self, stress, alpha=None, return_eps_ie=False):
		"""
		Computes the viscoplastic strain rate as

		.. math::
			\\dot{\\pmb{\\varepsilon}}_{vp} = \\mu_1 \\left\\langle \\dfrac{ F_{vp} }{F_0} \\right\\rangle^{N_1} \\dfrac{\\partial F_{vp}}{\\partial \\pmb{\\sigma}}

		where :math:`\\langle \\cdot \\rangle` is the ramp function.

		Parameters
		----------
		stress_vec : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing the stress tensor for all grid elements.

		alpha : torch.Tensor
			A (nelems,) pytorch tensor storing the hardening parameters for all grid elements.

		return_eps_ie : bool
			If *False* it stores :math:`\\dot{\\pmb{\\varepsilon}}_{ie}` to *eps_ie_rate*, otherwise the function returns as a value.
			This is used when computing :math:`\\frac{\\partial \\dot{\\pmb{\\varepsilon}}_{ie}}{\\partial \\pmb{\\sigma}}`
			by finite differences.

		"""
		if alpha == None:
			alpha = self._alpha

		s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.extract_stress_components(stress)
		I1, I2, I3, J2, J3, Sr, I1_star = self.compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)

		# Compute yield function
		Fvp = self.compute_Fvp(alpha, I1_star, J2, Sr)
		if not return_eps_ie:
			self._Fvp = Fvp.clone()

		# Compute flow direction, i.e. d(Fvp)/d(stress)
		F1 = (-alpha*I1**self._n + self._gamma*I1**2)
		F2 = (to.exp(self._beta_1*I1) - self._beta*Sr)
		dF1_dI1 = 2*self._gamma*I1 - self._n*alpha*I1**(self.n-1)
		dF2m_dI1 = self._beta_1*self._m*to.exp(self._beta_1*I1)*F2**(self._m-1)
		dF_dI1 = -(dF1_dI1*F2**self._m + F1*dF2m_dI1)

		dF2_dJ2 = -(3*self._beta*J3*27**0.5)/(4*J2**(5/2))
		dF_dJ2 = 1 - F1*self._m*F2**(self._m-1)*dF2_dJ2
		dF_dJ3 = -self._m*F1*self._beta*np.sqrt(27)*F2**(self._m-1)/(2*J2**1.5)

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

		dQdS = to.zeros_like(stress, dtype=to.float64)
		dQdS[:,0,0] = dQdS_00
		dQdS[:,1,1] = dQdS_11
		dQdS[:,2,2] = dQdS_22
		dQdS[:,1,0] = dQdS[:,0,1] = dQdS_01
		dQdS[:,2,0] = dQdS[:,0,2] = dQdS_02
		dQdS[:,2,1] = dQdS[:,1,2] = dQdS_12

		ramp_idx = to.where(Fvp > 0)[0]
		lmbda = to.zeros(self.n_elems, dtype=to.float64)
		lmbda[ramp_idx] = self._mu_1[ramp_idx]*(Fvp[ramp_idx]/self.F_0)**self._N_1[ramp_idx]
		# lmbda[ramp_idx] = self._mu_1[ramp_idx]*(Fvp[ramp_idx]/self.F_0[ramp_idx])**self._N_1[ramp_idx]
		eps_vp_rate = -dQdS*lmbda[:, None, None]

		if return_eps_ie:
			return eps_vp_rate
		else:
			self._eps_ie_rate = eps_vp_rate

	def compute_H(self, Q, P):
		"""
		This function computes the tensor product between **Q** and **P**, both rank-2 tensors,
		which implies that it results in a rank-4 tensor. That is,

		.. math::
			\\mathbb{H} = \\mathbf{Q} \\otimes \\mathbf{P}
			\\quad \\rightarrow \\quad
			H_{ijkl} = Q_{ij} P_{kl}

		This operation, however, is performed considering the symmetry of **Q** and **P**,
		which allows us to represent :math:`\\mathbb{H}` using Voigt notation.

		Parameters
		----------
		Q : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing **Q** for all grid elements.
		P : torch.Tensor
			A (nelems, 3, 3) pytorch tensor storing **P** for all grid elements.

		Returns
		-------
		H : torch.Tensor
			A (nelems, 6, 6) pytorch tensor storing :math:`\\mathbb{H}` for all grid elements. 
			Note that Voigt notation is employed here. 
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

	@property
	def mu_1(self):
		return self._mu_1

	@property
	def N_1(self):
		return self._N_1

	@property
	def a_1(self):
		return self._a_1

	@property
	def eta(self):
		return self._eta

	@property
	def n(self):
		return self._n

	@property
	def beta_1(self):
		return self._beta_1

	@property
	def beta(self):
		return self._beta

	@property
	def m(self):
		return self._m

	@property
	def gamma(self):
		return self._gamma

	@property
	def sigma_t(self):
		return self._sigma_t

	@property
	def alpha_0(self):
		return self._alpha_0

	@property
	def alpha(self):
		return self._alpha

	@property
	def Fvp(self):
		return self._Fvp

	@property
	def qsi(self):
		return self._qsi

	@property
	def qsi_old(self):
		return self._qsi_old

	@property
	def r(self):
		return self._r

	@property
	def eps_ie(self):
		return self._eps_ie

	@property
	def eps_ie_rate(self):
		return self._eps_ie_rate

	@property
	def eps_ie_rate_old(self):
		return self._eps_ie_rate_old

	@property
	def eps_ie_old(self):
		return self._eps_ie_old

	@property
	def eps_bar(self):
		return self._eps_bar

	@property
	def h(self):
		return self._h

	@property
	def Q(self):
		return self._Q

	@property
	def P(self):
		return self._P

	@property
	def B(self):
		return self._B

	@property
	def G(self):
		return self._G

	@property
	def H(self):
		return self._H

