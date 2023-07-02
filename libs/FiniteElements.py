from fenics import Measure, VectorFunctionSpace, TensorFunctionSpace, FunctionSpace, TrialFunction, TestFunction, dot, FacetNormal

class FemHandler():
	def __init__(self, grid):
		self.grid = grid
		self.dx_ = Measure("dx", domain=self.grid.mesh, subdomain_data=self.grid.subdomains)
		self.ds_ = Measure("ds", domain=self.grid.mesh, subdomain_data=self.grid.boundaries)
		self.V = VectorFunctionSpace(self.grid.mesh, "CG", 1)
		self.TS = TensorFunctionSpace(self.grid.mesh, "DG", 0)#, symmetry=True)
		self.P0 = FunctionSpace(self.grid.mesh, "DG", 0)
		self.du = TrialFunction(self.V)
		self.v = TestFunction(self.V)
		self.n = dot(self.v, FacetNormal(self.grid.mesh))

	def dx(self, DOMAIN_NAME=None):
		if DOMAIN_NAME == None:
			return self.dx_
		else:
			return self.dx_(self.grid.get_domain_tags(DOMAIN_NAME))

	def ds(self, BOUNDARY_NAME=None):
		return self.ds_(self.grid.get_boundary_tags(BOUNDARY_NAME))
		if BOUNDARY_NAME == None:
			return self.ds_
		else:
			return self.ds_(self.grid.get_boundary_tags(BOUNDARY_NAME))