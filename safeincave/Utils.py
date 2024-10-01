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
from dolfin import *
import json

MPa = 1e6
minute = 60
hour = 60*minute
day = 24*hour
year = 365*day

def read_json(file_name):
	"""
	This is just a convenient function to read json files.

	Parameters
	----------
	file_name : str
		Full path to json file, including file name, e.g. to/folder/file.json.

	Returns
	-------
	data : dict
		A dictionary containing the json data.
	"""
	with open(file_name, "r") as j_file:
		data = json.load(j_file)
	return data

def save_json(data, file_name):
	"""
	This is just a convenient function to save dictionaries into json files.

	Parameters
	----------
	data : dict
		Dictionary containing the data to be saved in the json file.

	file_name : str
		Full path to where the json file is supposed to be saved,including
		file name, e.g. to/folder/file.json.

	"""
	with open(file_name, "w") as f:
	    json.dump(data, f, indent=4)

def local_projection(tensor, V):
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dx
    b_proj = inner(tensor, v_)*dx
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    u = Function(V)
    solver.solve_local_rhs(u)
    return u

def epsilon(u):
	"""
	It computes the strain tensor based on the displacement field **u**, that is,

	.. math::
		\\pmb{\\varepsilon} (\\textbf{u}) = \\frac{1}{2} \\left( \\nabla \\textbf{u} + \\nabla \\textbf{u}^T \\right)

	Parameters
	----------
	u : dolfin.function.function.Function
		Displacement field.

	Returns
	-------
	grad_u : ufl.tensoralgebra.Sym
		Symmetric gradient of **u**.

	"""
	grad_u = sym(grad(u))
	return grad_u

def dotdot(C, eps):
	"""
	Computes the double dot product between **C** (4th-order tensor in Voigt notation) and
	**eps** (2nd-order tensor in tensor notation). It first converts **eps** to Voigt
	notation, then it performs the *dot* product with **C**, and finally converts the
	resulting quantity to tensor notation.

	Parameters
	----------
	C : dolfin.function.function.Function
		4th-order tensor in Voigt notation.

	eps : dolfin.function.function.Function
		2nd-order tensor in tensor notation.

	Returns
	-------
	tensor : ufl.tensors.ListTensor
		Double dot product between **C** and **eps**.

	"""
	tensor = voigt2tensor(dot(C, tensor2voigt(eps)))
	return tensor

def tensor2voigt(e):
	"""
	Converts tensor notation to Voigt notation.

	Parameters
	----------
	e : dolfin.function.function.Function
		A 2nd-order tensor in tensor notation.

	Returns
	-------
	e_voigt : ufl.tensors.ListTensor
		A 2nd-order tensor in Voigt notation.
	"""
	e_voigt = as_vector([e[0,0], e[1,1], e[2,2], e[0,1], e[0,2], e[1,2]])
	return e_voigt

def voigt2tensor(s):
	"""
	Converts a tensor from Voigt notation to tensor notation.

	Parameters
	----------
	s : ufl.tensors.ListTensor
		A 2nd-order tensor in Voigt notation.

	Returns
	-------
	s_tensor : ufl.tensors.ListTensor
		A 2nd-order tensor in tensor notation.
	"""
	s_tensor = as_matrix([[s[0], s[3], s[4]],
						  [s[3], s[1], s[5]],
						  [s[4], s[5], s[2]]])
	return s_tensor

def numpy2torch(numpy_array):
	"""
	It properly converts a numpy array into a pytorch tensor.

	Parameters
	----------
	numpy_array : numpy.ndarray
		Numpy array to be converted.

	Returns
	-------
	torch_array : torch.tensor
		A pytorch tensor with the same dimension as *numpy_array*.
	"""
	torch_array = to.tensor(numpy_array, dtype=to.float64)
	return torch_array

def dotdot2(C_voigt, eps_tensor):
	"""
	This function performs the double dot product between a 4th-order tensor represented 
	in Voigt notation and a 2nd-order tensor (without Voigt notation). The operation is
	performed by first applying the Voigt notation to the 2nd-order tensor, perform
	matrix-vector products to obtain the 2nd-order tensor in Voigt notation, and finally
	transform the resulting 2nd-order tensor to tensor notation.
	
	Parameters
	----------
	C_voigt : torch.Tensor
		This is a pytorch tensor storing the 4th-order tensors for all grid elements.
		Since Voigt notation is used, C_voigt has dimension (n_elems, 6, 6).

	eps_tensor : torch.Tensor
		This is a pytorch tensor storing 2nd-order tensors for all grid elements using
		tensor notation. Therefore, its dimensions are (n_elems, 3, 3).

	Returns
	-------
	stress_torch : torch.Tensor
		A tensor of dimensions (n_elems, 3, 3) resulting from the double dot product.


	"""
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

def compute_C(n_elems, nu, E):
	"""
	Assemble the :math:`\\mathbb{C}` matrices for each element of the grid.

	Parameters
	----------
	n_elems : int
		Number of grid elements.

	nu : list
		List containing Poisson's ratio values for each grid element.

	E : list
		List containing the value of Young's modulus for each grid element.

	Returns
	-------
	C : torch.Tensor
		Returns a (n_elems, 6, 6) containing matrices :math:`\\mathbb{C}` for all grid element.
	
	"""
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