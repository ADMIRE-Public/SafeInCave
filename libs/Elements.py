import torch as to
import numpy as np
from Utils import dotdot2, MPa

class Spring():
	def __init__(self, props_e, *args):
		self.E = props_e["E"]
		self.nu = props_e["nu"]
		try:
			self.n_elems = self.E.shape[0]
		except:
			self.n_elems = 1
		self.eps_e = to.tensor((self.n_elems, 3, 3), dtype=to.float64)


	def initialize(self):
		self.compute_C0()
		self.compute_C0_inv()

	def compute_C0(self):
		self.C0 = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		a0 = self.E/((1 + self.nu)*(1 - 2*self.nu))
		self.C0[:,0,0] = a0*(1 - self.nu)
		self.C0[:,1,1] = a0*(1 - self.nu)
		self.C0[:,2,2] = a0*(1 - self.nu)
		self.C0[:,3,3] = a0*(1 - 2*self.nu)
		self.C0[:,4,4] = a0*(1 - 2*self.nu)
		self.C0[:,5,5] = a0*(1 - 2*self.nu)
		self.C0[:,0,1] = self.C0[:,1,0] = self.C0[:,0,2] = self.C0[:,2,0] = self.C0[:,2,1] = self.C0[:,1,2] = a0*self.nu

	def compute_C0_inv(self):
		self.C0_inv = to.linalg.inv(self.C0)

	def compute_eps_e(self, stress):
		self.eps_e = dotdot2(self.C0_inv, stress)


class Viscoelastic():
	def __init__(self, props, psi_1=1.0, psi_2=1.0):
		self.psi_1 = psi_1
		self.psi_2 = psi_2
		self.eta = props["eta"]
		self.E = props["E"]
		self.nu = props["nu"]
		try:
			self.n_elems = self.E.shape[0]
		except:
			self.n_elems = 1
		self.eps_ve_rate = to.zeros((self.n_elems, 3, 3))
		self.eps_ve_rate_old = to.zeros((self.n_elems, 3, 3))
		self.eps_ve_old = to.zeros((self.n_elems, 3, 3))
		self.eps_ve = to.zeros((self.n_elems, 3, 3))
		self.eps_t = to.zeros((self.n_elems, 3, 3))

		self.B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.G = to.zeros((self.n_elems, 6, 6), dtype=to.float64)

	def initialize(self):
		self.compute_C1()

	def compute_C1(self):
		self.C1 = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		a0 = self.E/((1 + self.nu)*(1 - 2*self.nu))
		self.C1[:,0,0] = a0*(1 - self.nu)
		self.C1[:,1,1] = a0*(1 - self.nu)
		self.C1[:,2,2] = a0*(1 - self.nu)
		self.C1[:,3,3] = a0*(1 - 2*self.nu)
		self.C1[:,4,4] = a0*(1 - 2*self.nu)
		self.C1[:,5,5] = a0*(1 - 2*self.nu)
		self.C1[:,0,1] = self.C1[:,1,0] = self.C1[:,0,2] = self.C1[:,2,0] = self.C1[:,2,1] = self.C1[:,1,2] = a0*self.nu

	def compute_G_B(self, stress_vec, phi2, *args):
		I = to.eye(6, dtype=to.float64).unsqueeze(0).repeat(self.n_elems, 1, 1)
		self.G = to.linalg.inv(self.eta[:,None,None]*I + self.psi_1*phi2*self.C1)

	def increment_internal_variables(self, *args):
		pass

	def update_internal_variables(self, *args):
		pass

	def compute_eps_t(self, phi1, phi2):
		self.eps_t = self.eps_ve_old + phi1*self.eps_ve_rate_old + phi2*self.eps_ve_rate

	def compute_eps_ve(self, stress, stress_k, phi2):
		self.eps_ve = self.eps_t + phi2*dotdot2(self.G, stress - stress_k) - phi2*self.B

	def update_eps_ve_old(self):
		self.eps_ve_old = self.eps_ve.clone()

	def update_eps_ve_rate_old(self):
		self.eps_ve_rate_old = self.eps_ve_rate.clone()

	def compute_eps_ve_rate(self, stress_vec, phi1, return_eps_ve=False):
		eps_ve_rate = dotdot2(self.G, stress_vec - dotdot2(self.C1, self.eps_ve_old + phi1*self.eps_ve_rate_old))
		if return_eps_ve:
			return eps_ve_rate.clone()
		else:
			self.eps_ve_rate = eps_ve_rate.clone()



class DislocationCreep():
	def __init__(self, props, psi_1=1.0, psi_2=1.0):
		self.psi_1 = psi_1
		self.psi_2 = psi_2
		self.A = props["A"]*to.exp(-props["Q"]/props["R"]/props["T"])
		self.n = props["n"]
		try:
			self.n_elems = self.A.shape[0]
		except:
			self.n_elems = 1
		self.eps_ie_rate = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ie_rate_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ie_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_t = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.G = to.zeros((self.n_elems, 6, 6), dtype=to.float64)

	def increment_internal_variables(self, *args):
		pass

	def update_internal_variables(self, *args):
		pass

	def compute_eps_ie(self, stress, stress_k, phi2):
		self.eps_ie = self.eps_t + phi2*dotdot2(self.G, stress - stress_k) - phi2*self.B

	def update_eps_ie_old(self):
		self.eps_ie_old = self.eps_ie.clone()

	def update_eps_ie_rate_old(self):
		self.eps_ie_rate_old = self.eps_ie_rate.clone()

	def compute_eps_t(self, phi1, phi2):
		self.eps_t = self.eps_ie_old + phi1*self.eps_ie_rate_old + phi2*self.eps_ie_rate

	def compute_eps_ie_rate(self, stress_vec, return_eps_ie=False, *args):
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

		A_bar = self.A*q_vm**(self.n - 1)
		if return_eps_ie:
			return A_bar[:,None,None]*dev
		else:
			self.eps_ie_rate = A_bar[:,None,None]*dev

	def compute_E(self, stress_vec):
		EPSILON = 1e-2
		self.E = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		stress_eps = stress_vec.clone()
		c1 = 1.0
		c2 = 2.0
		magic_indexes = [(0,0,0,c1), (1,1,1,c1), (2,2,2,c1), (0,1,3,c2), (0,2,4,c2), (1,2,5,c2)]
		for i, j, k, phi in magic_indexes:
			stress_eps[:,i,j] += EPSILON
			eps_ie_rate_eps = self.compute_eps_ie_rate(stress_eps, return_eps_ie=True)
			self.E[:,:,k] = phi*(eps_ie_rate_eps[:,[0,1,2,0,0,1],[0,1,2,1,2,2]] - self.eps_ie_rate[:,[0,1,2,0,0,1],[0,1,2,1,2,2]]) / EPSILON
			stress_eps[:,i,j] -= EPSILON

	def compute_G_B(self, stress_vec, *args):
		if self.psi_1 == 1.0:
			self.compute_E(stress_vec)
			self.G = self.E.clone()
			# print(self.G)


class ViscoplasticDesai():
	def __init__(self, props, psi_1=1.0, psi_2=1.0):
		self.F_0 = 1.0
		# self.F_0 = props["F_0"]
		self.mu_1 = props["mu_1"]
		self.N_1 = props["N_1"]
		self.a_1 = props["a_1"]
		self.eta = props["eta"]
		self.n = props["n"]
		self.beta_1 = props["beta_1"]
		self.beta = props["beta"]
		self.m = props["m"]
		self.gamma = props["gamma"]
		self.sigma_t = props["sigma_t"]
		self.alpha_0 = props["alpha_0"]
		try:
			self.alpha = self.alpha_0.clone()
		except:
			self.alpha = self.alpha_0

		self.psi_1 = psi_1
		self.psi_2 = psi_2

		try:
			self.n_elems = self.alpha_0.shape[0]
		except:
			self.n_elems = 1

		self.Fvp = to.zeros(self.n_elems, dtype=to.float64)
		self.qsi = to.zeros(self.n_elems, dtype=to.float64)
		self.qsi_old = to.zeros(self.n_elems, dtype=to.float64)
		self.r = to.zeros(self.n_elems, dtype=to.float64)
		self.eps_ie_rate = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ie_rate_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ie_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_t = to.zeros((self.n_elems, 3, 3), dtype=to.float64)


	def compute_G_B(self, stress_vec, dt):
		# EPSILON_ALPHA = 1e-7
		EPSILON_ALPHA = 0.0001*self.alpha
		EPSILON_STRESS = 1e-1

		alpha_eps = self.alpha + EPSILON_ALPHA
		eps_ie_rate_eps = self.compute_eps_ie_rate(stress_vec, alpha=alpha_eps, return_eps_ie=True)

		self.r = self.compute_residue(self.eps_ie_rate, self.alpha, dt)
		r_eps = self.compute_residue(eps_ie_rate_eps, alpha_eps, dt)
		self.h = (r_eps - self.r) / EPSILON_ALPHA
		self.Q = (eps_ie_rate_eps - self.eps_ie_rate) / EPSILON_ALPHA[:,None,None]
		self.B = self.psi_1*(self.r / self.h)[:,None,None] * self.Q

		self.P = to.zeros_like(stress_vec)
		stress_eps = stress_vec.clone()
		for i, j in [(0,0), (1,1), (2,2), (0,1), (0,2), (1,2)]:
			stress_eps[:,i,j] += EPSILON_STRESS
			eps_ie_rate_eps = self.compute_eps_ie_rate(stress_eps, return_eps_ie=True)
			r_eps = self.compute_residue(eps_ie_rate_eps, self.alpha, dt)
			self.P[:,i,j] = (r_eps - self.r) / EPSILON_STRESS
			self.P[:,j,i] = self.P[:,i,j]
			stress_eps[:,i,j] -= EPSILON_STRESS

		self.QP = self.compute_QP(self.Q, self.P)

		self.compute_E(stress_vec)
		self.G = self.psi_1*(self.E - self.psi_2*self.QP/self.h[:,None,None])

		# for i, value in enumerate(self.h):
		# 	print(i, value)

	def compute_E(self, stress_vec):
		EPSILON_STRESS = 1e-1
		self.E = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		stress_eps = stress_vec.clone()
		c1 = 1.0
		c2 = 2.0
		magic_indexes = [(0,0,0,c1), (1,1,1,c1), (2,2,2,c1), (0,1,3,c2), (0,2,4,c2), (1,2,5,c2)]
		for i, j, k, phi in magic_indexes:
			stress_eps[:,i,j] += EPSILON_STRESS
			eps_ie_rate_eps = self.compute_eps_ie_rate(stress_eps, return_eps_ie=True)
			self.E[:,:,k] = phi*(eps_ie_rate_eps[:,[0,1,2,0,0,1],[0,1,2,1,2,2]] - self.eps_ie_rate[:,[0,1,2,0,0,1],[0,1,2,1,2,2]]) / EPSILON_STRESS
			stress_eps[:,i,j] -= EPSILON_STRESS
		self.E = self.E

	def compute_residue(self, eps_rate, alpha, dt):
		self.qsi = self.qsi_old + to.sum(eps_rate**2, axis=(-2, -1))**0.5*dt
		return alpha - self.a_1 / (((self.a_1/self.alpha_0)**(1/self.eta) + self.qsi)**self.eta)

	def compute_eps_t(self, phi1, phi2):
		self.eps_t = self.eps_ie_old + phi1*self.eps_ie_rate_old + phi2*self.eps_ie_rate

	def compute_eps_ie(self, stress, stress_k, phi2):
		self.eps_ie = self.eps_t + phi2*dotdot2(self.G, stress - stress_k) - phi2*self.B

	def update_eps_ie_old(self):
		self.eps_ie_old = self.eps_ie.clone()

	def update_eps_ie_rate_old(self):
		self.eps_ie_rate_old = self.eps_ie_rate.clone()

	def update_internal_variables(self):
		self.qsi_old = self.qsi.clone()

	def increment_internal_variables(self, stress, stress_k, dt):
		delta_alpha = -(self.r + self.psi_2*to.einsum('bij,bij->b', self.P, stress - stress_k))/self.h
		self.alpha += delta_alpha
		# self.alpha_q = self.alpha + self.k_v*(self.alpha_0 - self.alpha)*(1 - self.qsi_v/self.qsi)

	def compute_stress_invariants(self, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz):
		I1 = s_xx + s_yy + s_zz
		I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
		I3 = s_xx*s_yy*s_zz + 2*s_xy*s_yz*s_xz - s_zz*s_xy**2 - s_xx*s_yz**2 - s_yy*s_xz**2
		J2 = (1/3)*I1**2 - I2
		J3 = (2/27)*I1**3 - (1/3)*I1*I2 + I3
		Sr = -(J3*np.sqrt(27))/(2*J2**1.5)
		I1_star = I1 + self.sigma_t
		return I1, I2, I3, J2, J3, Sr, I1_star

	def extract_stress_components(self, stress):
		stress_vec = -stress
		s_xx = stress_vec[:,0,0]/MPa
		s_yy = stress_vec[:,1,1]/MPa
		s_zz = stress_vec[:,2,2]/MPa
		s_xy = stress_vec[:,0,1]/MPa
		s_xz = stress_vec[:,0,2]/MPa
		s_yz = stress_vec[:,1,2]/MPa
		return s_xx, s_yy, s_zz, s_xy, s_xz, s_yz

	def compute_Fvp(self, alpha, I1, J2, Sr):
		F1 = (alpha*I1**self.n - self.gamma*I1**2)
		F2 = (to.exp(self.beta_1*I1) - self.beta*Sr)
		Fvp = J2 + F1*F2**self.m
		return Fvp

	def compute_initial_hardening(self, stress, Fvp_0=-15.0):
		s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.extract_stress_components(stress)
		I1, I2, I3, J2, J3, Sr, I1_star = self.compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)
		self.alpha_0 =  self.gamma*I1_star**(2-self.n) + (Fvp_0 - J2)*I1_star**(-self.n)*(to.exp(self.beta_1*I1_star) - self.beta*Sr)**(-self.m)
		self.alpha = self.alpha_0.clone()

	def compute_eps_ie_rate(self, stress, alpha=None, return_eps_ie=False):
		if alpha == None:
			alpha = self.alpha

		s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.extract_stress_components(stress)
		I1, I2, I3, J2, J3, Sr, I1_star = self.compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)

		# Compute yield function
		Fvp = self.compute_Fvp(alpha, I1_star, J2, Sr)
		if not return_eps_ie:
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

		dQdS = to.zeros_like(stress, dtype=to.float64)
		dQdS[:,0,0] = dQdS_00
		dQdS[:,1,1] = dQdS_11
		dQdS[:,2,2] = dQdS_22
		dQdS[:,1,0] = dQdS[:,0,1] = dQdS_01
		dQdS[:,2,0] = dQdS[:,0,2] = dQdS_02
		dQdS[:,2,1] = dQdS[:,1,2] = dQdS_12

		ramp_idx = to.where(Fvp > 0)[0]
		lmbda = to.zeros(self.n_elems, dtype=to.float64)
		lmbda[ramp_idx] = self.mu_1[ramp_idx]*(Fvp[ramp_idx]/self.F_0)**self.N_1[ramp_idx]
		# lmbda[ramp_idx] = self.mu_1[ramp_idx]*(Fvp[ramp_idx]/self.F_0[ramp_idx])**self.N_1[ramp_idx]
		eps_vp_rate = -dQdS*lmbda[:, None, None]

		if return_eps_ie:
			return eps_vp_rate
		else:
			self.eps_ie_rate = eps_vp_rate

	def compute_QP(self, Q, P):
		n_elems, _, _ = P.shape
		QP = to.zeros((n_elems, 6, 6), dtype=to.float64)
		QP[:,0,0] = Q[:,0,0]*P[:,0,0]
		QP[:,0,1] = Q[:,0,0]*P[:,1,1]
		QP[:,0,2] = Q[:,0,0]*P[:,2,2]
		QP[:,0,3] = 2*Q[:,0,0]*P[:,0,1]
		QP[:,0,4] = 2*Q[:,0,0]*P[:,0,2]
		QP[:,0,5] = 2*Q[:,0,0]*P[:,1,2]

		QP[:,1,0] = Q[:,1,1]*P[:,0,0]
		QP[:,1,1] = Q[:,1,1]*P[:,1,1]
		QP[:,1,2] = Q[:,1,1]*P[:,2,2]
		QP[:,1,3] = 2*Q[:,1,1]*P[:,0,1]
		QP[:,1,4] = 2*Q[:,1,1]*P[:,0,2]
		QP[:,1,5] = 2*Q[:,1,1]*P[:,1,2]

		QP[:,2,0] = Q[:,2,2]*P[:,0,0]
		QP[:,2,1] = Q[:,2,2]*P[:,1,1]
		QP[:,2,2] = Q[:,2,2]*P[:,2,2]
		QP[:,2,3] = 2*Q[:,2,2]*P[:,0,1]
		QP[:,2,4] = 2*Q[:,2,2]*P[:,0,2]
		QP[:,2,5] = 2*Q[:,2,2]*P[:,1,2]

		QP[:,3,0] = Q[:,0,1]*P[:,0,0]
		QP[:,3,1] = Q[:,0,1]*P[:,1,1]
		QP[:,3,2] = Q[:,0,1]*P[:,2,2]
		QP[:,3,3] = 2*Q[:,0,1]*P[:,0,1]
		QP[:,3,4] = 2*Q[:,0,1]*P[:,0,2]
		QP[:,3,5] = 2*Q[:,0,1]*P[:,1,2]

		QP[:,4,0] = Q[:,0,2]*P[:,0,0]
		QP[:,4,1] = Q[:,0,2]*P[:,1,1]
		QP[:,4,2] = Q[:,0,2]*P[:,2,2]
		QP[:,4,3] = 2*Q[:,0,2]*P[:,0,1]
		QP[:,4,4] = 2*Q[:,0,2]*P[:,0,2]
		QP[:,4,5] = 2*Q[:,0,2]*P[:,1,2]

		QP[:,5,0] = Q[:,1,2]*P[:,0,0]
		QP[:,5,1] = Q[:,1,2]*P[:,1,1]
		QP[:,5,2] = Q[:,1,2]*P[:,2,2]
		QP[:,5,3] = 2*Q[:,1,2]*P[:,0,1]
		QP[:,5,4] = 2*Q[:,1,2]*P[:,0,2]
		QP[:,5,5] = 2*Q[:,1,2]*P[:,1,2]
		return QP
