# Copyright 2025 The safeincave community.
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

from numpy.typing import NDArray
from typing import Callable
from .Grid import GridHandlerGMSH
import numpy as np
import torch as to
import dolfinx as do
import ufl
import json

# Type aliases
UFLVector3 = ufl.core.expr.Expr  # Shape (3,)
UFLVector6 = ufl.core.expr.Expr  # Shape (6,)
UFLTensor3x3 = ufl.core.expr.Expr  # Shape (3,3)
UFLMatrix6x6 = ufl.core.expr.Expr  # Shape (6,6)
Fn = Callable[[float, float, float], float]

# Useful units
GPa = 1e9
MPa = 1e6
kPa = 1e3
minute = 60
hour = 60*minute
day = 24*hour
year = 365*day

def read_json(file_name: str) -> dict:
	"""
	Read a JSON file into a Python dictionary.

	Parameters
	----------
	file_name : str
	    Path to the JSON file.

	Returns
	-------
	dict
	    Parsed JSON content.
	"""
	with open(file_name, "r") as j_file:
		data = json.load(j_file)
	return data

def save_json(data: dict, file_name: str) -> None:
	"""
	Save a Python dictionary to a JSON file.

	Parameters
	----------
	data : dict
	    Data to serialize as JSON.
	file_name : str
	    Output file path.

	Returns
	-------
	None

	Notes
	-----
	Uses an indentation of 4 spaces for readability. Overwrites the
	file if it already exists.
	"""
	with open(file_name, "w") as f:
		json.dump(data, f, indent=4)

def project(tensor_ufl: ufl.core.expr.Expr, V: do.fem.FunctionSpace) -> do.fem.Function:
	"""
	Interpolate a UFL expression into a DOLFINx function space.

	Parameters
	----------
	tensor_ufl : ufl.core.expr.Expr
	    UFL expression to be interpolated (e.g., rank-2 tensor field).
	V : dolfinx.fem.FunctionSpace
	    Target function space. Its element must be compatible with the
	    UFL expression shape.

	Returns
	-------
	dolfinx.fem.Function
	    Function in `V` obtained by interpolation of `tensor_ufl`.

	Notes
	-----
	This uses `dolfinx.fem.Expression` with `V.element.interpolation_points()`
	and then calls `Function.interpolate(...)`. For non-interpolatory elements,
	prefer projection/solve-based approaches.
	"""
	tensor_expr = do.fem.Expression(tensor_ufl, V.element.interpolation_points())
	tensor = do.fem.Function(V)
	tensor.interpolate(tensor_expr)
	return tensor

def epsilon(u: UFLVector3) -> UFLTensor3x3:
	"""
	Small-strain tensor :math:`\\varepsilon(u) = \\mathrm{sym}(\\nabla u)`.

	Parameters
	----------
	u : UFLVector3
	    Displacement field as a UFL vector expression with shape (3,).

	Returns
	-------
	UFLTensor3x3
	    Symmetric gradient of `u`, a (3, 3) UFL tensor.

	Raises
	------
	AssertionError
	    If `u` is not rank-1 with shape (3,).

	Notes
	-----
	UFLVector3 = ufl.core.expr.Expr
	"""
	assert ufl.rank(u) == 1 and u.ufl_shape == (3,)
	grad_u = ufl.sym(ufl.grad(u))
	return grad_u

def dotdot_ufl(C: UFLMatrix6x6, eps: UFLTensor3x3) -> UFLTensor3x3:
	"""
	Compute stress :math:`\\sigma = C : \\varepsilon` using Voigt mapping in UFL.

	Parameters
	----------
	C : UFLMatrix6x6
	    Fourth-order stiffness tensor stored in Voigt form as a (6, 6) UFL object.
	    Assumes **tensorial Voigt** convention (no engineering shear factors).
	eps : UFLTensor3x3
	    Small-strain tensor as a (3, 3) UFL expression.

	Returns
	-------
	UFLTensor3x3
	    Cauchy stress tensor as a (3, 3) UFL expression.

	Raises
	------
	AssertionError
	    If shapes are not `(6, 6)` for `C` or `(3, 3)` for `eps`.

	Notes
	-----
	Internally computes `voigt2tensor(ufl.dot(C, tensor2voigt(eps)))`.
	UFLMatrix6x6 = ufl.core.expr.Expr
	UFLTensor3x3 = ufl.core.expr.Expr
	"""
	assert ufl.rank(C) == 2 and C.ufl_shape == (6, 6)
	assert ufl.rank(eps) == 2 and eps.ufl_shape == (3, 3)
	tensor = voigt2tensor(ufl.dot(C, tensor2voigt(eps)))
	return tensor

def tensor2voigt(e: UFLTensor3x3) -> UFLVector6:
	"""
    Map a 3×3 symmetric tensor to Voigt vector (tensorial Voigt).

    Parameters
    ----------
    e : UFLTensor3x3
        Tensor expression with shape (3, 3). Symmetry is assumed but not enforced.

    Returns
    -------
    UFLVector6
        Voigt vector `[e_xx, e_yy, e_zz, e_xy, e_xz, e_yz]`.

    Raises
    ------
    AssertionError
        If `e` is not rank-2 with shape (3, 3).

    Notes
    -----
    multiply the shear components by 0.5 when mapping back to tensor form.
    UFLTensor3x3 = ufl.core.expr.Expr
    """
	assert ufl.rank(e) == 2 and e.ufl_shape == (3, 3)
	e_voigt = ufl.as_vector([e[0,0], e[1,1], e[2,2], e[0,1], e[0,2], e[1,2]])
	return e_voigt

def voigt2tensor(s: UFLVector6) -> UFLTensor3x3:
	"""
    Map a Voigt vector (length 6) to a symmetric 3×3 tensor (tensorial Voigt).

    Parameters
    ----------
    s : UFLVector6
        Voigt vector `[xx, yy, zz, xy, xz, yz]` (no engineering shear factors).

    Returns
    -------
    UFLTensor3x3
        Symmetric (3, 3) tensor reconstructed from `s`.

    Raises
    ------
    AssertionError
        If the UFL shape of `s` is not `(6,)`.

    Notes
    -----
    UFLVector6 = ufl.core.expr.Expr
    UFLTensor3x3 = ufl.core.expr.Expr
    """
	assert s.ufl_shape == (6,)
	s_tensor = ufl.as_matrix([[s[0], s[3], s[4]],
							  [s[3], s[1], s[5]],
							  [s[4], s[5], s[2]]])
	return s_tensor

def numpy2torch(numpy_array: NDArray[np.float64]) -> to.Tensor:
	"""
    Convert a NumPy array to a `torch.Tensor` with dtype `torch.float64`.

    Parameters
    ----------
    numpy_array : numpy.typing.NDArray[numpy.float64]
        Input NumPy array of any shape.

    Returns
    -------
    torch.Tensor
        Tensor with the same shape as `numpy_array`, cast to `torch.float64`.

    Notes
    -----
    Uses `torch.tensor`, which **copies** data. For zero-copy semantics,
    prefer `torch.from_numpy(numpy_array)` when the dtype is already `float64`.
    """
	torch_array = to.tensor(numpy_array, dtype=to.float64)
	return torch_array

def dotdot_torch(C_voigt: to.Tensor, eps_tensor: to.Tensor) -> to.Tensor:
	"""
    σ = C : ε using Voigt notation (tensorial Voigt, no engineering factors).

    Parameters
    ----------
    C_voigt : (N, 6, 6) torch.Tensor (float64)
        Stiffness matrix per element in tensorial Voigt order [xx, yy, zz, xy, xz, yz].
    eps_tensor : (N, 3, 3) torch.Tensor (float64)
        Small-strain tensor per element.

    Returns
    -------
    (N, 3, 3) torch.Tensor (float64)
        Cauchy stress tensor per element.
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

def create_field_nodes(grid: GridHandlerGMSH, fun: Fn) -> to.Tensor:
	"""
    Sample a scalar field at mesh nodes using a Python callable.

    Parameters
    ----------
    grid : GridHandlerGMSH
        Grid handler with attributes `mesh` and `n_nodes`. Assumes
        `grid.mesh.geometry.x` provides an array of node coordinates of shape (N, 3).
    fun : Callable[[float, float, float], float]
        Function evaluated as `fun(x, y, z)` at each node.

    Returns
    -------
    torch.Tensor
        1D tensor of length `grid.n_nodes` with dtype `torch.float64`.

    Notes
    -----
    This evaluates `fun` in a Python loop; it is not vectorized.
    """
	coordinates = grid.mesh.geometry.x
	field = to.zeros(grid.n_nodes, dtype=to.float64)
	for i, coord in enumerate(coordinates):
		x, y, z = coord
		field[i] = fun(x, y, z)
	return field

def create_field_elems(grid: GridHandlerGMSH, fun: Fn) -> to.Tensor:
	"""
    Sample a scalar field at element centroids using a Python callable.

    Parameters
    ----------
    grid : GridHandlerGMSH
        Grid handler with attributes `mesh` and `n_elems`. Assumes a tetrahedral
        mesh where the 3→0 connectivity returns 4 vertex indices per cell.
    fun : Callable[[float, float, float], float]
        Function evaluated as `fun(x, y, z)` at each cell centroid.

    Returns
    -------
    torch.Tensor
        1D tensor of length `grid.n_elems` with dtype `torch.float64`.

    Notes
    -----
    The centroid is computed as the arithmetic mean of the 4 vertex coordinates
    for each tetrahedral cell. This evaluates `fun` in a Python loop.
    """
	field = to.zeros(grid.n_elems, dtype=to.float64)
	coordinates = grid.mesh.geometry.x
	conn_aux = grid.mesh.topology.connectivity(3, 0)
	conn = conn_aux.array.reshape((grid.n_elems, 4))
	for i in range(grid.n_elems):
		cell_vertices = conn[i]
		x = sum(coordinates[v] for v in cell_vertices) / len(cell_vertices)
		field[i] = fun(x[0], x[1], x[2])
	return field