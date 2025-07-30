from abc import ABC
import numpy as np
import dolfinx as do
import ufl

class GeneralBC(ABC):
	def __init__(self):
		self.boundary_name = None
		self.type = None
		self.values = None
		self.time_values = None


class DirichletBC(GeneralBC):
	def __init__(self, boundary_name : str, component : int, values : list, time_values : list):
		self.boundary_name = boundary_name
		self.type = "dirichlet"
		self.values = values
		self.time_values = time_values
		self.component = component

class NeumannBC(GeneralBC):
	def __init__(self, boundary_name : str, direction : int, density : float, ref_pos : float, values : list, time_values : list, g=-9.91):
		self.boundary_name = boundary_name
		self.type = "neumann"
		self.values = values
		self.time_values = time_values
		self.direction = direction
		self.density = density
		self.ref_pos = ref_pos
		self.gravity = g


class BcHandler():
	def __init__(self, equation):
		self.eq = equation
		self.dirichlet_boundaries = []
		self.neumann_boundaries = []
		self.x = ufl.SpatialCoordinate(self.eq.grid.mesh)

	def reset_boundary_conditions(self):
		self.dirichlet_boundaries = []
		self.neumann_boundaries = []

	def add_boundary_condition(self, bc : GeneralBC):
		if bc.type == "dirichlet":
			self.dirichlet_boundaries.append(bc)
		elif bc.type == "neumann":
			self.neumann_boundaries.append(bc)
		else:
			raise Exception(f"Boundary type {bc.type} not supported.")

	def update_dirichlet(self, t):
		self.dirichlet_bcs = []
		for bc in self.dirichlet_boundaries:
			value = np.interp(t, bc.time_values, bc.values)
			dofs = do.fem.locate_dofs_topological(
				self.eq.get_uV().sub(bc.component),
				self.eq.grid.boundary_dim,
				self.eq.grid.get_boundary_tags(bc.boundary_name)
			)
			self.dirichlet_bcs.append(
				do.fem.dirichletbc(
					do.default_scalar_type(value),
					dofs,
					self.eq.get_uV().sub(bc.component)
				)
			)

	def update_neumann(self, t):
		self.neumann_bcs = []
		for bc in self.neumann_boundaries:
			i = bc.direction
			rho = bc.density
			H = bc.ref_pos
			p = -np.interp(t, bc.time_values, bc.values)
			value_neumann = p + rho*bc.gravity*(H - self.x[i])
			self.neumann_bcs.append(value_neumann*self.eq.normal*self.eq.ds(self.eq.grid.get_boundary_tag(bc.boundary_name)))

