import torch as to
from Utils import *

class LinearMomentum():
	"""
	This class is intended to solve the linear momemum balance equation considering
	a general non-linear constitutive models.

	Parameters
	----------
	m : :class:`safeincave.ConstitutiveModel.ConstitutiveModel`
		Object containing the constitutive model data.
	"""
	def __init__(self, m, theta):
		self.m = m
		self.n_elems = m.n_elems
		self.theta = theta

		self.eps_e = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_bar_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_bar_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.GT_ve = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.GT_ie = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.BT_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.CT = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.stress = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.stress_k = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

	def compute_GT_BT_ve(self, dt):
		"""
		Computes matrices :math:`\\mathbb{G}_T` and :math:`\\mathbf{B}_T` for the
		viscoelastic elements. That is

		.. math::

			\\mathbb{G}_T = \\sum_{i=1}^{N_{ve}} \\mathbb{G}_i
			\\quad \\text{and} \\quad
			\\mathbf{B}_T = \\sum_{i=1}^{N_{ve}} \\mathbf{B}_i

		where :math:`N_{ve}` denotes the number of viscoelastic (Kelvin-Voigt)
		elements present in the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
		None 
		"""
		self.GT_ve = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		phi2 = dt*(1 - self.theta)
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_G_B(self.stress, phi2)
			self.GT_ve += elem_ve.G
			self.BT_ve += elem_ve.B

	def compute_GT_BT_ie(self, dt):
		"""
		Computes matrices :math:`\\mathbb{G}_{ie,T}` and :math:`\\mathbf{B}_{ie,T}` for the
		inelastic elements. That is

		.. math::

			\\mathbb{G}_{ie,T} = \\sum_{i=1}^{N_{ie}} \\mathbb{G}_i
			\\quad \\text{and} \\quad
			\\mathbf{B}_{ie,T} = \\sum_{i=1}^{N_{ie}} \\mathbf{B}_i

		where :math:`N_{ie}` denotes the number of inelastic elements present
		in the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
		None 
		"""
		self.GT_ie = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_G_B(self.stress, dt)
			self.GT_ie += elem_ie.G
			self.BT_ie += elem_ie.B

	def compute_CT(self, dt):
		"""
		Computes the consistent tangent matrix :math:`\\mathbb{C}_{T}`, which is given by

		.. math::

			\\mathbb{C}_{T} = \\left( 
				\\mathbb{C}_0^{-1} + \\Delta t (1 - \\theta) \\mathbb{G}_{T}
			\\right)^{-1}

		where :math:`\\mathbb{C}_0` is the stiffness matrix associated to the linear spring,
		and :math:`\\mathbb{G}_{T} = \\mathbb{G}_{ve,T} + \\mathbb{G}_{ie,T}`.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
		None 
		"""
		GT = self.GT_ie + self.GT_ve
		self.CT = to.linalg.inv(self.m.C0_inv + dt*(1-self.theta)*GT)

	def compute_eps_ve(self, dt):
		"""
		Computes the strains of all **viscoelastic** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.eps_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_eps_ve(self.stress, self.stress_k, dt*(1-self.theta))
			self.eps_ve += elem_ve.eps_ve

	def compute_eps_ie(self, dt):
		"""
		Computes the strains of all **inelastic** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.eps_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_eps_ie(self.stress, self.stress_k, dt*(1-self.theta))
			self.eps_ie += elem_ie.eps_ie

	def compute_eps_e(self):
		"""
		Computes the strains of all **elastic** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.eps_e = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_e in self.m.elems_e:
			elem_e.compute_eps_e(self.stress)
			self.eps_e += elem_e.eps_e

	def compute_eps_bar_ie(self, dt):
		"""
		Computes the eps_bar of all **inelastic** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.eps_bar_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_eps_bar(dt*self.theta, dt*(1 - self.theta))
			self.eps_bar_ie += elem_ie.eps_bar

	def compute_eps_bar_ve(self, dt):
		"""
		Computes the eps_bar of all **viscoelastic** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.eps_bar_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_eps_bar(dt*self.theta, dt*(1 - self.theta))
			self.eps_bar_ve += elem_ve.eps_bar

	def compute_eps_bar(self, dt):
		"""
		Computes the eps_bar of **all** elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.compute_eps_bar_ie(dt)
		self.compute_eps_bar_ve(dt)
		self.eps_bar = self.eps_bar_ve + self.eps_bar_ie

	def compute_eps_rhs(self, dt, *args):
		"""
		Computes the term :math:`\\pmb{\\varepsilon}_{rhs}^k` of the linearized momentum equation, that is,

		.. math::

			\\pmb{\\varepsilon}_\\text{rhs}^k = \\bar{\\pmb{\\varepsilon}}_{ne}^k - \\Delta t \\left( 1 - \\theta \\right) \\left( \\mathbb{G}_{ne} : \\pmb{\\sigma}^k + \\mathbf{B}_{ne} \\right)

		Parameters
		----------
		dt : float
			Time step size.

		Returns
		-------
			None
		"""
		self.compute_eps_bar(dt)
		BT = self.BT_ie + self.BT_ve
		GT = self.GT_ie + self.GT_ve
		self.eps_rhs = self.eps_bar - dt*(1-self.theta)*(BT + dotdot2(GT, self.stress_k))

	def compute_stress_C0(self, eps_e):
		"""
		Compute stress tensor as

		.. math::

			\\pmb{\\sigma} = \\mathbb{C}_0 : \\pmb{\\varepsilon}_e

		.. note::
			This operation considers that :math:`\\mathbb{C}_0` is represented by Voigt notation.

		Parameters
		----------
		eps_e : torch.Tensor
			A (nelems, 3, 3) storing the elastic strain for all grid elements.
		"""
		self.stress = dotdot2(self.m.C0, eps_e)

	def compute_stress(self, eps_tot, dt):
		"""
		Compute stress tensor as

		.. math::

			\\pmb{\\sigma}^{k+1} = \\mathbb{C}_T :
										    \\left[
										        \\pmb{\\varepsilon}^{k+1}
										        - \\bar{\\pmb{\\varepsilon}}_{ne}^k
										        + \\Delta t (1 - \\theta)
										            \\left( 
										               \\mathbf{B}_{ne}
										               + \\mathbb{G}_{ne} : \\pmb{\\sigma}^k
										            \\right)
										    \\right]
		
		Parameters
		----------
		eps_tot : torch.Tensor
			A (nelems, 3, 3) storing the total strain for all grid elements.
		"""
		self.compute_eps_bar(dt)
		GT = self.GT_ie + self.GT_ve
		BT = self.BT_ie + self.BT_ve
		self.stress = dotdot2(self.CT, eps_tot - self.eps_bar + dt*(1-self.theta)*(BT + dotdot2(GT, self.stress_k)))

	def update_stress(self):
		"""
		Updates stress of the previous interation, that is :math:`\\pmb{\\sigma}^k \\leftarrow \\pmb{\\sigma}^{k+1}`.
		"""
		self.stress_k = self.stress.clone()

	def compute_eps_ie_rate(self):
		"""
		Computes strain rates of all inelastic elements of the constitutive model.
		"""
		for elem_ie in self.m.elems_ie:
			elem_ie.compute_eps_ie_rate(self.stress, return_eps_ie=False)

	def update_eps_ie_rate_old(self):
		"""
		Updates inelastic strain rates of the previous time level,
		that is :math:`\\dot{\\pmb{\\varepsilon}}^t_{ie} \\leftarrow \\dot{\\pmb{\\varepsilon}}^{t + \\Delta t}_{ie}`.
		"""
		for elem_ie in self.m.elems_ie:
			elem_ie.update_eps_ie_rate_old()

	def update_eps_ie_old(self):
		"""
		Updates inelastic strains of the previous time level,
		that is :math:`\\pmb{\\varepsilon}^t_{ie} \\leftarrow \\pmb{\\varepsilon}^{t + \\Delta t}_{ie}`.
		"""
		for elem_ie in self.m.elems_ie:
			elem_ie.update_eps_ie_old()

	def compute_eps_ve_rate(self, dt):
		"""
		Computes the strain rates of the viscoelastic elements of the constitutive model.

		Parameters
		----------
		dt : float
			Time step size.
		"""
		for elem_ve in self.m.elems_ve:
			elem_ve.compute_eps_ve_rate(self.stress, dt*self.theta, return_eps_ve=False)

	def update_eps_ve_rate_old(self):
		"""
		Updates viscoelastic strain rates of the previous time level,
		that is :math:`\\dot{\\pmb{\\varepsilon}}^t_{ve} \\leftarrow \\dot{\\pmb{\\varepsilon}}^{t + \\Delta t}_{ve}`.
		"""
		for elem_ve in self.m.elems_ve:
			elem_ve.update_eps_ve_rate_old()

	def update_eps_ve_old(self):
		"""
		Updates viscoelastic strains of the previous time level,
		that is :math:`\\pmb{\\varepsilon}^t_{ve} \\leftarrow \\pmb{\\varepsilon}^{t + \\Delta t}_{ve}`.
		"""
		for elem_ve in self.m.elems_ve:
			elem_ve.update_eps_ve_old()

	def increment_internal_variables(self, dt):
		"""
		Increment internal variables.

		Parameters
		----------
		dt : float
			Time step size.
		"""
		for elem_ie in self.m.elems_ie:
			elem_ie.increment_internal_variables(self.stress, self.stress_k, dt)

	def update_internal_variables(self):
		"""
		Update internal variables with values of previous iteration, that is
		:math:`\\alpha_i^k \\leftarrow \\alpha_i^{k+1}.`
		"""
		for elem_ie in self.m.elems_ie:
			elem_ie.update_internal_variables()









