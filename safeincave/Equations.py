import torch as to
from Utils import *

class LinearMomentum():
	def __init__(self, m, theta):
		"""
			m: object of ConstitutiveModel class
			theta: float number defining the time integration method 
						-> 0.0 for fully implicit time integration
						-> 0.5 for Crank-Nicholson time integration
						-> 1.0 for explicit time integration
		"""
		self.m = m
		self.n_elems = m.n_elems
		self.theta = theta

		self.eps_e = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_t_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_t_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.GT_ve = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.GT_ie = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.BT_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.CT = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.stress = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.stress_k = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

	def compute_GT_BT_ve(self, dt):
		self.GT_ve = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		phi2 = dt*(1 - self.theta)
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_G_B(self.stress, phi2)
			self.GT_ve += elem_ve.G
			self.BT_ve += elem_ve.B

	def compute_GT_BT_ie(self, dt):
		self.GT_ie = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_G_B(self.stress, dt)
			self.GT_ie += elem_ie.G
			self.BT_ie += elem_ie.B

	def compute_CT(self, dt):
		GT = self.GT_ie + self.GT_ve
		self.CT = to.linalg.inv(self.m.C0_inv + dt*(1-self.theta)*GT)

	def compute_eps_ve(self, dt):
		self.eps_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_eps_ve(self.stress, self.stress_k, dt*(1-self.theta))
			self.eps_ve += elem_ve.eps_ve

	def compute_eps_ie(self, dt):
		self.eps_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_eps_ie(self.stress, self.stress_k, dt*(1-self.theta))
			self.eps_ie += elem_ie.eps_ie

	def compute_eps_e(self):
		self.eps_e = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_e in self.m.elems_e:
			elem_e.compute_eps_e(self.stress)
			self.eps_e += elem_e.eps_e

	def compute_eps_t_ie(self, dt):
		self.eps_t_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_eps_t(dt*self.theta, dt*(1 - self.theta))
			self.eps_t_ie += elem_ie.eps_t

	def compute_eps_t_ve(self, dt):
		self.eps_t_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_eps_t(dt*self.theta, dt*(1 - self.theta))
			self.eps_t_ve += elem_ve.eps_t

	def compute_eps_t(self, dt):
		self.compute_eps_t_ie(dt)
		self.compute_eps_t_ve(dt)
		self.eps_t = self.eps_t_ve + self.eps_t_ie

	def compute_eps_rhs(self, dt, *args):
		self.compute_eps_t(dt)
		BT = self.BT_ie + self.BT_ve
		GT = self.GT_ie + self.GT_ve
		self.eps_rhs = self.eps_t - dt*(1-self.theta)*(BT + dotdot2(GT, self.stress_k))

	def compute_stress_C0(self, eps_e):
		self.stress = dotdot2(self.m.C0, eps_e)

	def compute_stress(self, eps_tot, dt):
		self.compute_eps_t(dt)
		GT = self.GT_ie + self.GT_ve
		BT = self.BT_ie + self.BT_ve
		self.stress = dotdot2(self.CT, eps_tot - self.eps_t + dt*(1-self.theta)*(BT + dotdot2(GT, self.stress_k)))

	def update_stress(self):
		self.stress_k = self.stress.clone()

	def compute_eps_ie_rate(self):
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_eps_ie_rate(self.stress, return_eps_ie=False)

	def update_eps_ie_rate_old(self):
		for elem_ie in self.m.elems_ie:
			elem_ie.update_eps_ie_rate_old()

	def update_eps_ie_old(self):
		for elem_ie in self.m.elems_ie:
			elem_ie.update_eps_ie_old()

	def compute_eps_ve_rate(self, dt):
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_eps_ve_rate(self.stress, dt*self.theta, return_eps_ve=False)

	def update_eps_ve_rate_old(self):
		for elem_ve in self.m.elems_ve:
			elem_ve.update_eps_ve_rate_old()

	def update_eps_ve_old(self):
		for elem_ve in self.m.elems_ve:
			elem_ve.update_eps_ve_old()

	def increment_internal_variables(self, dt):
		for elem_ie in self.m.elems_ie:
			elem_ie.increment_internal_variables(self.stress, self.stress_k, dt)

	def update_internal_variables(self):
		for elem_ie in self.m.elems_ie:
			elem_ie.update_internal_variables()









