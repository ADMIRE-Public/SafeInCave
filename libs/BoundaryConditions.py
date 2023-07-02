from fenics import Constant, assemble, DirichletBC
import abc

class BoundaryConditions(metaclass=abc.ABCMeta):
	@abc.abstractmethod
	def __init__(self, fem_handler, settings):
		self.fem_handler = fem_handler
		self.boundaryConditions = settings["BoundaryConditions"]

	@abc.abstractmethod
	def update_Dirichlet_BC(self):
		pass

	@abc.abstractmethod
	def update_Neumann_BC(self):
		pass

	@abc.abstractmethod
	def update_BCs(self):
		pass

class MechanicsBoundaryConditions(BoundaryConditions):
	def __init__(self, fem_handler, settings):
		super().__init__(fem_handler, settings)

	def update_Dirichlet_BC(self, time_handler):
		self.bcs = []
		index_dict = {"u_x": 0, "u_y": 1, "u_z": 2}
		for key_0, value_0 in self.boundaryConditions.items():
			index_u = index_dict[key_0]
			for BOUNDARY_NAME, VALUES in value_0.items():
				if VALUES["type"] == "DIRICHLET":
					self.bcs.append(
									DirichletBC(
												self.fem_handler.V.sub(index_u),
												Constant(VALUES["value"][time_handler.idx]),
												self.fem_handler.grid.get_boundaries(),
												self.fem_handler.grid.get_boundary_tags(BOUNDARY_NAME)
									)
					)

	def update_Neumann_BC(self, time_handler):
		L_bc = 0
		for key_0, value_0 in self.boundaryConditions.items():
			for BOUNDARY_NAME, VALUES in value_0.items():
				if VALUES["type"] == "NEUMANN":
					load = Constant(VALUES["value"][time_handler.idx])
					L_bc += load*self.fem_handler.n*self.fem_handler.ds(BOUNDARY_NAME)
		self.b = assemble(L_bc)

	def update_BCs(self, time_handler):
		self.update_Neumann_BC(time_handler)
		self.update_Dirichlet_BC(time_handler)