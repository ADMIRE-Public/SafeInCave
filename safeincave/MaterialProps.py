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
from abc import ABC, abstractmethod
import torch as to
import numpy as np
from Utils import dotdot2, MPa

class Material():
	def __init__(self, n_elems):
		self.n_elems = n_elems
		self.elems_ne = []
		self.elems_th = []
		self.elems_e = []

		self.C_inv = to.zeros((n_elems, 6, 6), dtype=to.float64)
		self.C = to.zeros((n_elems, 6, 6), dtype=to.float64)

		self.C_tilde_inv = to.zeros((n_elems, 6, 6), dtype=to.float64)
		self.C_tilde = to.zeros((n_elems, 6, 6), dtype=to.float64)

	def set_density(self, density : to.Tensor):
		self.density = density

	def set_specific_heat_capacity(self, cp : to.Tensor):
		self.cp = cp

	def set_thermal_conductivity(self, k : to.Tensor):
		self.k = k

	def set_thermal_expansion(self, alpha_th : to.Tensor):
		self.alpha_th = alpha_th


	def add_to_elastic(self, elem):
		elem.initialize()
		self.C_inv += elem.C_inv
		self.C += elem.C
		self.C_tilde_inv += elem.C_tilde_inv
		self.C_tilde += elem.C_tilde
		self.elems_e.append(elem)
		self.K = elem.K
		self.E = elem.E
		self.ShearMod = 3*self.K*self.E/(9*self.K - self.E)

	def add_to_non_elastic(self, elem):
		self.elems_ne.append(elem)

	def add_to_thermoelastic(self, elem):
		self.elems_th.append(elem)

	def compute_G_B(self, stress, dt, theta, T):
		self.G = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ne in self.elems_ne:
			elem_ne.compute_G_B(stress, dt, theta, T)
			self.G += elem_ne.G
			self.B += elem_ne.B

	def compute_T_IT(self):
		self.IT = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.T = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ne in self.elems_ne:
			elem_ne.compute_T_IT()
			self.IT += elem_ne.IT
			self.T += elem_ne.T

	def compute_Bvol_Tvol(self, stress, dt):
		self.B_vol = to.zeros(self.n_elems, dtype=to.float64)
		self.T_vol = to.zeros(self.n_elems, dtype=to.float64)
		for elem_ne in self.elems_ne:
			elem_ne.compute_Bvol_Tvol()
			self.B_vol += elem_ne.B_vol
			self.T_vol += elem_ne.T_vol

	def compute_Gtilde_Btilde(self, stress, dt):
		self.G_tilde = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.B_tilde = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ne in self.elems_ne:
			elem_ne.compute_Gtilde_Btilde()
			self.G_tilde += elem_ne.G_tilde
			self.B_tilde += elem_ne.B_tilde

	def compute_CT(self, dt, theta):
		self.CT = to.linalg.inv(self.C_inv + dt*(1-theta)*self.G)

	def compute_CT_tilde(self, dt, theta):
		self.CT_tilde = to.linalg.inv(self.C_tilde_inv + dt*(1-theta)*self.G_tilde)


class Thermoelastic():
	def __init__(self, alpha, name="thermoelastic"):
		self.alpha = alpha
		self.name = name
		self.n_elems = self.alpha.shape[0]
		self.eps_th = to.zeros((self.n_elems, 3, 3))
		self.I = to.eye(3, dtype=to.float64).unsqueeze(0).repeat(self.n_elems, 1, 1)

	def compute_eps_th(self, dT_DG_vec):
		self.eps_th = self.alpha[:,None,None]*dT_DG_vec[:,None,None]*self.I


class Spring():
	def __init__(self, E, nu, name="spring"):
		self.E = E
		self.nu = nu
		self.name = name
		self.n_elems = self.E.shape[0]
		self.eps_e = to.tensor((self.n_elems, 3, 3), dtype=to.float64)

	def initialize(self):
		self.C = self.__compute_C(self.n_elems, self.nu, self.E)
		self.C_inv = self.__compute_C_inv(self.C)
		self.C_tilde = self.__compute_C_tilde(self.n_elems, self.nu, self.E)
		self.C_tilde_inv = self.__compute_C_tilde_inv(self.n_elems, self.nu, self.E)
		self.K = self.E/(3*(1 - 2*self.nu))

	def compute_eps_e(self, stress):
		self.eps_e = dotdot2(self.C_inv, stress)

	def __compute_C(self, n_elems, nu, E):
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
		return to.linalg.inv(self.C)

	def __compute_C_tilde(self, n_elems, nu, E):
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
	def __init__(self, n_elems):
		self.n_elems = n_elems
		self.eps_ne_rate = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ne_rate_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ne_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ne_k = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.G = to.zeros((self.n_elems, 6, 6), dtype=to.float64)

	@abstractmethod
	def compute_eps_ne_rate(self, stress_vec, phi1, Temp, return_eps_ne=False):
		pass

	def increment_internal_variables(self, *args):
		pass

	def update_internal_variables(self, *args):
		pass

	def compute_eps_ne_k(self, phi1, phi2):
		self.eps_ne_k = self.eps_ne_old + phi1*self.eps_ne_rate_old + phi2*self.eps_ne_rate

	def update_eps_ne_old(self, stress, stress_k, phi2):
		self.eps_ne_old = self.eps_ne_k + phi2*dotdot2(self.G, stress - stress_k) - phi2*self.B

	def update_eps_ne_rate_old(self):
		self.eps_ne_rate_old = self.eps_ne_rate.clone()

	def compute_E(self, stress, dt, theta, Temp):
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

	def compute_B_and_H_over_h(self, stress, dt, theta, Temp):
		B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		H_over_h = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		return B, H_over_h

	def compute_G_B(self, stress, dt, theta, Temp):
		self.B, H_over_h = self.compute_B_and_H_over_h(stress, dt, theta, Temp)
		E = self.compute_E(stress, dt, theta, Temp)
		self.G = E.clone() - H_over_h.clone()

	def compute_T_IT(self):
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

	def compute_Bvol_Tvol(self):
		self.T_vol = to.einsum("bii->b", self.T)
		self.B_vol = to.einsum("bii->b", self.B)

	def compute_Gtilde_Btilde(self):
		I = to.eye(3).expand(self.n_elems, -1, -1)
		self.G_tilde = self.G - self.IT/3
		self.B_tilde = self.B - self.B_vol[:,None,None]*I/3

		



class Viscoelastic(NonElasticElement):
	def __init__(self, eta, E, nu, name="kelvin_voigt"):
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

	def compute_eps_ne_rate(self, stress_vec, phi1, Temp, return_eps_ne=False):
		eps_ne_rate = dotdot2(self.G, stress_vec - dotdot2(self.C1, self.eps_ne_old + phi1*self.eps_ne_rate_old))
		if return_eps_ne:
			return eps_ne_rate.clone()
		else:
			self.eps_ne_rate = eps_ne_rate.clone()

	def compute_E(self, stress, dt, theta, Temp):
		phi2 = dt*(1 - theta)
		I = to.eye(6, dtype=to.float64).unsqueeze(0).repeat(self.n_elems, 1, 1)
		E = to.linalg.inv(self.eta[:,None,None]*I + phi2*self.C1)
		return E




class DislocationCreep(NonElasticElement):
	def __init__(self, A, Q, n, name="creep"):
		super().__init__(A.shape[0])
		self.R = 8.32
		self.Q = Q
		self.A = A
		self.n = n
		self.name = name

	def compute_eps_ne_rate(self, stress_vec, phi1, Temp, return_eps_ne=False):
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
	def __init__(self, mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0, name="desai"):
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

	def compute_residue(self, eps_rate, alpha, dt):
		self.qsi = self.qsi_old + to.sum(eps_rate**2, axis=(-2, -1))**0.5*dt
		return alpha - self.a_1 / (((self.a_1/self.alpha_0)**(1/self.eta) + self.qsi)**self.eta)

	def update_internal_variables(self):
		self.qsi_old = self.qsi.clone()

	def increment_internal_variables(self, stress, stress_k, dt):
		delta_alpha = -(self.r + to.einsum('bij,bij->b', self.P, stress - stress_k))/self.h
		self.alpha += delta_alpha

	def __compute_stress_invariants(self, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz):
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

	def __extract_stress_components(self, stress):
		stress_vec = -stress
		s_xx = stress_vec[:,0,0]/MPa
		s_yy = stress_vec[:,1,1]/MPa
		s_zz = stress_vec[:,2,2]/MPa
		s_xy = stress_vec[:,0,1]/MPa
		s_xz = stress_vec[:,0,2]/MPa
		s_yz = stress_vec[:,1,2]/MPa
		return s_xx, s_yy, s_zz, s_xy, s_xz, s_yz

	def __compute_Fvp(self, alpha, I1, J2, Sr):
		F1 = (alpha*I1**self.n - self.gamma*I1**2)
		F2 = (to.exp(self.beta_1*I1) - self.beta*Sr)
		Fvp = J2 + F1*F2**self.m
		return Fvp

	def compute_initial_hardening(self, stress, Fvp_0=0.0):
		s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.__extract_stress_components(stress)
		I1, I2, I3, J2, J3, Sr, I1_star, _ = self.__compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)
		self.alpha_0 =  self.gamma*I1_star**(2-self.n) + (Fvp_0 - J2)*I1_star**(-self.n)*(to.exp(self.beta_1*I1_star) - self.beta*Sr)**(-self.m)
		self.alpha = self.alpha_0.clone()

	def compute_eps_ne_rate(self, stress, phi1, Temp, alpha=None, return_eps_ne=False):
		if alpha == None:
			alpha = self.alpha

		s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.__extract_stress_components(stress)
		I1, I2, I3, J2, J3, Sr, I1_star, ind_J2_leq_0 = self.__compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)

		# Compute yield function
		Fvp = self.__compute_Fvp(alpha, I1_star, J2, Sr)
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
		ramp_idx = list(to.where(Fvp > 0)[0])
		lmbda = to.zeros(self.n_elems, dtype=to.float64)
		if len(ramp_idx) != 0:
			lmbda[ramp_idx] = self.mu_1[ramp_idx]*(Fvp[ramp_idx]/self.F_0)**self.N_1[ramp_idx]
		eps_vp_rate = -dQdS*lmbda[:, None, None]

		if return_eps_ne:
			return eps_vp_rate
		else:
			self.eps_ne_rate = eps_vp_rate


	def compute_B_and_H_over_h(self, stress, dt, theta, Temp):
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


	def compute_H(self, Q, P):
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

















