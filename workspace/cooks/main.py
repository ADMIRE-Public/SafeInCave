import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
from Grid import GridHandlerGMSH
# import dolfin as do
# import numpy as np
# import pandas as pd
# from time import time
# from itertools import product
# import json

kPa = 1e3
MPa = 1e6

def local_projection(tensor, V):
    dv = do.TrialFunction(V)
    v_ = do.TestFunction(V)
    a_proj = do.inner(dv, v_)*do.dx
    b_proj = do.inner(tensor, v_)*do.dx
    solver = do.LocalSolver(a_proj, b_proj)
    solver.factorize()
    u = do.Function(V)
    solver.solve_local_rhs(u)
    return u

def save_json(data, file_name):
	with open(file_name, "w") as f:
	    json.dump(data, f, indent=4)

def read_json(file_name):
    with open(file_name, "r") as j_file:
        data = json.load(j_file)
    return data



def solve_mixed_FEM(grid, E, nu, h):
	E = do.Constant(E)
	nu = do.Constant(nu)
	K = do.Constant(E/(3*(1 - 2*nu)))
	lame = do.Constant(E*nu/((1+nu)*(1-2*nu)))
	G = do.Constant(E/2/(1+nu))
	Q = (G + lame)/(G*K)

	DG_3x3 = do.TensorFunctionSpace(grid.mesh, "DG", 0)
	DG_1x1 = do.FunctionSpace(grid.mesh, "DG", 0)

	# Define variational problem
	Vue = do.VectorElement('CG', grid.mesh.ufl_cell(), 1) # displacement finite element
	Vte = do.FiniteElement('CG', grid.mesh.ufl_cell(), 1) # pressure finite element
	V = do.FunctionSpace(grid.mesh, Vue*Vte)
	U_ = do.TestFunction(V)
	dU = do.TrialFunction(V)

	ds = do.Measure("ds", domain=grid.mesh, subdomain_data=grid.get_boundaries())
	dx = do.Measure("dx", domain=grid.mesh, subdomain_data=grid.get_subdomains())

	(v_u, v_T) = do.split(U_)
	normal = do.dot(v_u, do.FacetNormal(grid.mesh))

	bc1 = do.DirichletBC(V.sub(0).sub(2), do.Constant(0.), grid.get_boundaries(), grid.get_boundary_tags("BOTTOM"))
	bc2 = do.DirichletBC(V.sub(0).sub(2), do.Constant(0.), grid.get_boundaries(), grid.get_boundary_tags("TOP"))
	bc3 = do.DirichletBC(V.sub(0).sub(0), do.Constant(0.), grid.get_boundaries(), grid.get_boundary_tags("WEST"))
	bc4 = do.DirichletBC(V.sub(0).sub(1), do.Constant(0.), grid.get_boundaries(), grid.get_boundary_tags("WEST"))
	bcs = [bc1, bc2, bc3, bc4]

	shear_load = 10.0*kPa
	T = do.Constant((0, shear_load, 0))
	integral_neumann = []
	integral_neumann.append(do.dot(T, v_u)*ds(grid.get_boundary_tags("EAST")))

	(u_, p_) = do.split(U_)
	(du, dp) = do.split(dU)
	Uold = do.Function(V)
	(uold, pold) = do.split(Uold)

	def eps(v):
	    return do.sym(do.grad(v))

	def sigma_up(v, p):
	    return 2*G*(eps(v) - (1/3)*do.tr(eps(v))*do.Identity(3)) + p*do.Identity(3)

	# Calculate characteristic lenght, h
	h_cell_2 = do.Function(DG_1x1)
	h_cell_2.vector()[:] = h**2
	stab = (1/(3*K) + 1/G)*h_cell_2


	form = do.inner(sigma_up(du, dp), eps(u_))*do.dx
	form += (dp*p_ - K*do.tr(eps(du))*p_ + do.dot(K*stab*do.grad(dp), do.grad(p_)) )*do.dx
	form += -sum(integral_neumann)

	U = do.Function(V)
	do.solve(
			do.lhs(form) == do.rhs(form), 
			U, 
			bcs, 
			solver_parameters={
             "linear_solver": "mumps",
             "preconditioner": "lu"}
	)
	(u, P) = U.split(deepcopy=True)
	n_unknowns = U.vector()[:].size

	# Calculate stresses
	stress = local_projection(sigma_up(u, P), DG_3x3)

	return u, stress, n_unknowns

def run_simulation(grid, h, nu, output_folder):
	start = time()
	print(f"Output folder: {output_folder}")

	# Read grid
	n_elems = grid.mesh.num_cells()
	n_nodes = grid.mesh.num_vertices()

	# Define elastic properties
	E = 170e3

	# Solve problem
	u, stress, n_unknowns = solve_mixed_FEM(grid, E, nu, h)
	u.rename("Displacement", "m")

	# Define function spaces
	CG_1x1 = do.FunctionSpace(grid.mesh, "CG", 1)
	DG_1x1 = do.FunctionSpace(grid.mesh, "DG", 0)

	# Calculate mean and deviatoric stresses
	stress_np = stress.vector()[:].reshape((n_elems,3,3))
	s_xx = stress_np[:,0,0]
	s_yy = stress_np[:,1,1]
	s_zz = stress_np[:,2,2]
	s_xy = stress_np[:,0,1]
	s_xz = stress_np[:,0,2]
	s_yz = stress_np[:,1,2]

	# Save results at cells
	q_cell_vtk = do.File(os.path.join(output_folder, "q_cell", "q_cell.pvd"))
	sigv_cell_vtk = do.File(os.path.join(output_folder, "sigv_cell", "sigv_cell.pvd"))
	stress_cell_vtk = do.File(os.path.join(output_folder, "stress_cell", "stress_cell.pvd"))

	sigv_cell = do.project(do.tr(stress), DG_1x1)
	sigv_cell.rename("Mean Stress", "Pa")
	sigv_cell_vtk << sigv_cell

	q_cell = do.Function(DG_1x1)
	q_cell.vector()[:] = np.sqrt( 0.5*( (s_xx - s_yy)**2 + (s_xx - s_zz)**2 + (s_yy - s_zz)**2 + 6*(s_xy**2 + s_xz**2 + s_yz**2) ) )
	q_cell.rename("Von Mises Stress", "Pa")
	q_cell_vtk << q_cell

	stress_cell_vtk << stress

	# Save results at nodes
	u_node_vtk = do.File(os.path.join(output_folder, "u_node", "u_node.pvd"))
	q_node_vtk = do.File(os.path.join(output_folder, "q_node", "q_node.pvd"))
	sigv_node_vtk = do.File(os.path.join(output_folder, "sigv_node", "sigv_node.pvd"))
	sxx_node_vtk = do.File(os.path.join(output_folder, "sxx_node", "sxx_node.pvd"))
	syy_node_vtk = do.File(os.path.join(output_folder, "syy_node", "syy_node.pvd"))
	szz_node_vtk = do.File(os.path.join(output_folder, "szz_node", "szz_node.pvd"))
	sxy_node_vtk = do.File(os.path.join(output_folder, "sxy_node", "sxy_node.pvd"))
	sxz_node_vtk = do.File(os.path.join(output_folder, "sxz_node", "sxz_node.pvd"))
	syz_node_vtk = do.File(os.path.join(output_folder, "syz_node", "syz_node.pvd"))

	u_node_vtk << u
	q_node = do.project(q_cell, CG_1x1)
	q_node.rename("Von Mises Stress", "MPa")
	q_node_vtk << q_node
	sigv_node = do.project(sigv_cell, CG_1x1)
	sigv_node.rename("Mean Stress", "MPa")
	sigv_node_vtk << sigv_node
	sxx_node_vtk << do.project(stress[0,0], CG_1x1)
	syy_node_vtk << do.project(stress[1,1], CG_1x1)
	szz_node_vtk << do.project(stress[2,2], CG_1x1)
	sxy_node_vtk << do.project(stress[0,1], CG_1x1)
	sxz_node_vtk << do.project(stress[0,2], CG_1x1)
	syz_node_vtk << do.project(stress[1,2], CG_1x1)

	finish = time()
	print(f"CPU Time: {finish - start} seconds.\n")

	data_time = {
		"Time": finish - start,
		"n_unknowns": n_unknowns,
		"n_elems": n_elems,
		"n_nodes": n_nodes
	}
	save_json(data_time, os.path.join(output_folder, "data_time.json"))


def main_1():
	nu = 0.4999

	grid_name = "grid_1"
	grid_path = os.path.join("grids", grid_name)
	grid = GridHandlerGMSH("geom", grid_path)

	print(grid.mesh.num_cells())
	print(grid.mesh.num_vertices())

	# # Read h
	# h_dict = read_json(os.path.join(grid_path, "h.json"))
	# h0 = np.zeros_like(np.array(h_dict["h_opt"]))
	# h_dict["h0"] = h0
	# print(h_dict.keys())
	# df = pd.DataFrame(h_dict)
	# print(df)

	# # for h_name in h_dict.keys():
	# # # for h_name in ["h_ML"]:
	# # 	combinations = list(product([0, 1], repeat=2))
	# # 	print(combinations)
	# # 	for combo in combinations:
	# # 		a2, a3 = combo
	# # 		output_folder = os.path.join("output", "stab", f"nu_{nu}", f"P1P1{h_name}", f"stab_{a2}_{a3}")
	# # 		h = np.array(h_dict[h_name])
	# # 		run_simulation(grid, h, nu, output_folder, deg_u=1, deg_p=1, a2=a2, a3=a3)

	# # output_folder = os.path.join("output", "stab", f"nu_{nu}", f"P1")
	# # run_simulation(grid, h0, nu, output_folder, deg_u=1, deg_p=-1)

	# output_folder = os.path.join("output", "stab", f"nu_{nu}", f"P2")
	# run_simulation(grid, h0, nu, output_folder, deg_u=2, deg_p=-1)

	# output_folder = os.path.join("output", "stab", f"nu_{nu}", f"P2P1")
	# run_simulation(grid, h0, nu, output_folder, deg_u=2, deg_p=1)

if __name__ == '__main__':
	main_1()