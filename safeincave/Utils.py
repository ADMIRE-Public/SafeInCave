"""
Useful functions used in safeincave.
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
import dolfinx as do
from dolfinx.fem.petsc import LinearProblem
import ufl
import json

GPa = 1e9
MPa = 1e6
kPa = 1e3
minute = 60
hour = 60*minute
day = 24*hour
year = 365*day

def read_json(file_name):
	with open(file_name, "r") as j_file:
		data = json.load(j_file)
	return data

def save_json(data, file_name):
	with open(file_name, "w") as f:
	    json.dump(data, f, indent=4)

def local_projection_old(tensor, V):
    u = do.fem.Function(V)
    dv = ufl.TrialFunction(V)
    v_ = ufl.TestFunction(V)
    a_proj = ufl.inner(dv, v_)*ufl.dx
    b_proj = ufl.inner(tensor, v_)*ufl.dx
    problem = LinearProblem(a_proj, b_proj, u=u)
    problem.solve()
    return u

def project(tensor_ufl, V):
	tensor_expr = do.fem.Expression(tensor_ufl, V.element.interpolation_points())
	tensor = do.fem.Function(V)
	tensor.interpolate(tensor_expr)
	return tensor

def epsilon(u):
	grad_u = ufl.sym(ufl.grad(u))
	return grad_u

def dotdot(C, eps):
	tensor = voigt2tensor(ufl.dot(C, tensor2voigt(eps)))
	return tensor

def tensor2voigt(e):
	e_voigt = ufl.as_vector([e[0,0], e[1,1], e[2,2], e[0,1], e[0,2], e[1,2]])
	return e_voigt

def voigt2tensor(s):
	s_tensor = ufl.as_matrix([[s[0], s[3], s[4]],
							  [s[3], s[1], s[5]],
							  [s[4], s[5], s[2]]])
	return s_tensor

def numpy2torch(numpy_array):
	torch_array = to.tensor(numpy_array, dtype=to.float64)
	return torch_array

def dotdot2(C_voigt, eps_tensor):
	n_elems = C_voigt.shape[0]
	eps_voigt = to.zeros((n_elems, 6), dtype=to.float64)
	eps_voigt[:,0] = eps_tensor[:,0,0]
	eps_voigt[:,1] = eps_tensor[:,1,1]
	eps_voigt[:,2] = eps_tensor[:,2,2]
	eps_voigt[:,3] = eps_tensor[:,0,1]
	eps_voigt[:,4] = eps_tensor[:,0,2]
	eps_voigt[:,5] = eps_tensor[:,1,2]
	stress_voigt = to.bmm(C_voigt, eps_voigt.unsqueeze(2)).squeeze(2)
	stress_torch = to.zeros_like(eps_tensor, dtype=to.float64)
	stress_torch[:,0,0] = stress_voigt[:,0]
	stress_torch[:,1,1] = stress_voigt[:,1]
	stress_torch[:,2,2] = stress_voigt[:,2]
	stress_torch[:,0,1] = stress_torch[:,1,0] = stress_voigt[:,3]
	stress_torch[:,0,2] = stress_torch[:,2,0] = stress_voigt[:,4]
	stress_torch[:,1,2] = stress_torch[:,2,1] = stress_voigt[:,5]
	return stress_torch

def dotdot3(eps_tensor, C_voigt):
	Q = C_voigt.clone()
	Q[:,[3,4,5]] /= 2
	Q = Q.transpose(1,2)
	Q[:,[3,4,5]] *= 2
	n_elems = Q.shape[0]
	eps_voigt = to.zeros((n_elems, 6), dtype=to.float64)
	eps_voigt[:,0] = eps_tensor[:,0,0]
	eps_voigt[:,1] = eps_tensor[:,1,1]
	eps_voigt[:,2] = eps_tensor[:,2,2]
	eps_voigt[:,3] = eps_tensor[:,0,1]
	eps_voigt[:,4] = eps_tensor[:,0,2]
	eps_voigt[:,5] = eps_tensor[:,1,2]
	stress_voigt = to.bmm(Q, eps_voigt.unsqueeze(2)).squeeze(2)
	stress_torch = to.zeros_like(eps_tensor, dtype=to.float64)
	stress_torch[:,0,0] = stress_voigt[:,0]
	stress_torch[:,1,1] = stress_voigt[:,1]
	stress_torch[:,2,2] = stress_voigt[:,2]
	stress_torch[:,0,1] = stress_torch[:,1,0] = stress_voigt[:,3]
	stress_torch[:,0,2] = stress_torch[:,2,0] = stress_voigt[:,4]
	stress_torch[:,1,2] = stress_torch[:,2,1] = stress_voigt[:,5]
	return stress_torch

def create_field_nodes(grid, fun):
	coordinates = grid.mesh.geometry.x
	field = to.zeros(grid.n_nodes, dtype=to.float64)
	for i, coord in enumerate(coordinates):
		x, y, z = coord
		field[i] = fun(x, y, z)
	return field

def create_field_elems(grid, fun):
	field = to.zeros(grid.n_elems, dtype=to.float64)
	coordinates = grid.mesh.geometry.x
	conn_aux = grid.mesh.topology.connectivity(3, 0)
	conn = conn_aux.array.reshape((grid.n_elems, 4))
	for i in range(grid.n_elems):
		cell_vertices = conn[i]
		x = sum(coordinates[v] for v in cell_vertices) / len(cell_vertices)
		field[i] = fun(x[0], x[1], x[2])
	return field