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

import meshio as ms
import numpy as np
import os
from scipy.sparse import csr_matrix


def build_smoother(points: np.ndarray, conn: np.ndarray) -> np.ndarray:
    """
    Build a node–cell smoothing operator for tetrahedral meshes.

    Given node coordinates and tetrahedral connectivities, this constructs a
    (sparse) smoothing matrix that averages cell-centered quantities to nodes
    and back to cells via a volume-weighted scheme.

    Parameters
    ----------
    points : (n_nodes, 3) ndarray of float
        Node coordinates (x, y, z).
    conn : (n_elems, 4) ndarray of int
        Tetrahedral connectivity (node indices per element).

    Returns
    -------
    smoother : (n_elems, n_elems) scipy.sparse.csr_matrix
        Sparse smoothing operator mapping cell fields to smoothed cell fields.
        (Note: although the annotation says ``np.ndarray``, the returned
        object is a CSR sparse matrix.)

    Notes
    -----
    The operator is assembled as ``B * A`` where:
    - ``A`` distributes cell values to nodes using volume weights of incident
      cells.
    - ``B`` averages nodal values back to cells with uniform weights over the
      element's nodes.

    See Also
    --------
    scipy.sparse.csr_matrix
    """
    def tetrahedron_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
        """Compute signed volume of a tetrahedron given its four vertices."""
        volume = abs((1/6) * ((x2 - x1) * ((y3 - y1)*(z4 - z1) - (z3 - z1)*(y4 - y1)) + 
                     (y2 - y1) * ((z3 - z1)*(x4 - x1) - (x3 - x1)*(z4 - z1)) + 
                     (z2 - z1) * ((x3 - x1)*(y4 - y1) - (y3 - y1)*(x4 - x1))))
        return volume

    def build_node_elem_stencil(conn, coord):
        """Build node-to-incident-element adjacency lists."""
        stencil = [[] for i in range(n_nodes)]
        for elem, elem_conn in enumerate(conn):
            for node in elem_conn:
                if elem not in stencil[node]:
                    stencil[node].append(elem)
        return stencil

    # Initialize
    n_elems = conn.shape[0]
    n_nodes = points.shape[0]

    # Calculate volumes of all tetrahedra
    volumes = np.zeros(n_elems)
    for i in range(n_elems):
        nodes = conn[i]
        x1, y1, z1 = points[nodes[0], 0], points[nodes[0], 1], points[nodes[0], 2]
        x2, y2, z2 = points[nodes[1], 0], points[nodes[1], 1], points[nodes[1], 2]
        x3, y3, z3 = points[nodes[2], 0], points[nodes[2], 1], points[nodes[2], 2]
        x4, y4, z4 = points[nodes[3], 0], points[nodes[3], 1], points[nodes[3], 2]
        volumes[i] = tetrahedron_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

    stencil = build_node_elem_stencil(conn, points)

    A_row, A_col, A_data = [], [], []
    for node in range(n_nodes):
        vol = volumes[stencil[node]].sum()
        for elem in stencil[node]:
            A_row.append(node)
            A_col.append(elem)
            A_data.append(volumes[elem]/vol)
    A_csr = csr_matrix((A_data, (A_row, A_col)), shape=(n_nodes, n_elems))

    B_row, B_col, B_data = [], [], []
    for elem, nodes in enumerate(conn):
        for node in nodes:
            B_row.append(elem)
            B_col.append(node)
            B_data.append(1/len(nodes))
    B_csr = csr_matrix((B_data, (B_row, B_col)), shape=(n_elems, n_nodes))
    smoother = B_csr.dot(A_csr)

    return smoother

def build_mapping(nodes_xdmf: np.ndarray, nodes_msh: np.ndarray) -> list[int]:
    """
    Build an index mapping from XDMF node order to MSH node order.

    For each coordinate triplet in ``nodes_xdmf``, finds the row index in
    ``nodes_msh`` with the exact same coordinates and returns the list of
    corresponding indices.

    Parameters
    ----------
    nodes_xdmf : (n_nodes, 3) ndarray of float
        Node coordinates as read from an XDMF file.
    nodes_msh : (n_nodes, 3) ndarray of float
        Node coordinates as read from a .msh file.

    Returns
    -------
    mapping : list of int
        For each row in ``nodes_xdmf``, the index of the identical row in
        ``nodes_msh``.

    Notes
    -----
    This uses exact floating-point equality. If the two sources differ by
    round-off, consider a tolerance-based nearest matching instead.
    """
    return [np.where((nodes_msh == row).all(axis=1))[0][0] for row in nodes_xdmf]

def find_closest_point(target_point: np.ndarray, points: np.ndarray) -> int:
    """
    Find the index of the closest point in a set to a target point.

    Parameters
    ----------
    target_point : (3,) ndarray of float
        Query point (x, y, z).
    points : (n_points, 3) ndarray of float
        Candidate points.

    Returns
    -------
    idx : int
        Index of the closest point in ``points`` (Euclidean distance).
    """
    x_p, y_p, z_p = target_point
    d = np.sqrt(  (points[:,0] - x_p)**2
                + (points[:,1] - y_p)**2
                + (points[:,2] - z_p)**2 )
    p_idx = d.argmin()
    return p_idx


def compute_cell_centroids(cells: np.ndarray, points: np.ndarray) -> np.ndarray:
    """
    Compute centroids of tetrahedral cells.

    Parameters
    ----------
    cells : (n_cells, 4) ndarray of int
        Tetrahedral connectivity (node indices per cell).
    points : (n_nodes, 3) ndarray of float
        Node coordinates (x, y, z).

    Returns
    -------
    centroids : (n_cells, 3) ndarray of float
        Centroid coordinates for each cell, computed as the arithmetic mean
        of its four vertex coordinates.
    """
    n_cells = cells.shape[0]
    centroids = np.zeros((n_cells, 3))
    for i, cell in enumerate(cells):
        p0 = points[cell[0]]
        p1 = points[cell[1]]
        p2 = points[cell[2]]
        p3 = points[cell[3]]
        x = (p0[0] + p1[0] + p2[0] + p3[0])/4
        y = (p0[1] + p1[1] + p2[1] + p3[1])/4
        z = (p0[2] + p1[2] + p2[2] + p3[2])/4
        centroids[i,:] = np.array([x, y, z])
    return centroids


def read_cell_tensor(xdmf_field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a time series of cell-centered 3x3 tensor fields from an XDMF file.

    Parameters
    ----------
    xdmf_field_path : str
        Path to the XDMF file containing cell data (``cells['tetra']``).

    Returns
    -------
    centroids : (n_cells, 3) ndarray of float
        Centroid coordinates of the tetrahedral cells.
    time_list : (n_steps,) ndarray of float
        Time values for each time step.
    tensor_field : (n_steps, n_cells, 3, 3) ndarray of float
        Tensor values per time step and cell.

    Notes
    -----
    The function assumes a single tensor field is present under
    ``cell_data['tetra']`` at each time step, and reshapes it to (3, 3) per cell.
    """
    reader = ms.xdmf.TimeSeriesReader(xdmf_field_path)
    points, cells = reader.read_points_cells()
    n_cells = cells["tetra"].shape[0]
    n_steps = reader.num_steps

    centroids = compute_cell_centroids(cells["tetra"], points)
    tensor_field = np.zeros((n_steps, n_cells, 3, 3))
    time_list = np.zeros(n_steps)

    for k in range(reader.num_steps):
        # Read data
        time, point_data, cell_data = reader.read_data(k)

        # Add time
        time_list[k] = time

        # Add tensor
        field_name = list(cell_data["tetra"].keys())[0]
        tensor_field[k,:,:] = cell_data["tetra"][field_name].reshape((n_cells, 3, 3))

    return centroids, time_list, tensor_field


def read_cell_scalar(xdmf_field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a time series of cell-centered scalar fields from an XDMF file.

    Parameters
    ----------
    xdmf_field_path : str
        Path to the XDMF file containing cell data (``cells['tetra']``).

    Returns
    -------
    centroids : (n_cells, 3) ndarray of float
        Centroid coordinates of the tetrahedral cells.
    time_list : (n_steps,) ndarray of float
        Time values for each time step.
    scalar_field : (n_steps, n_cells) ndarray of float
        Scalar values per time step and cell.

    Notes
    -----
    The function assumes a single scalar field is present under
    ``cell_data['tetra']`` at each time step.
    """
    reader = ms.xdmf.TimeSeriesReader(xdmf_field_path)

    points, cells = reader.read_points_cells()
    n_cells = cells["tetra"].shape[0]
    n_steps = reader.num_steps

    centroids = compute_cell_centroids(cells["tetra"], points)
    scalar_field = np.zeros((n_steps, n_cells))
    time_list = np.zeros(n_steps)

    for k in range(reader.num_steps):
        # Read data
        time, point_data, cell_data = reader.read_data(k)

        # Add time
        time_list[k] = time

        # Add scalar
        field_name = list(cell_data["tetra"].keys())[0]
        scalar_field[k] = cell_data["tetra"][field_name].flatten()

    return centroids, time_list, scalar_field


def read_node_scalar(xdmf_field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a time series of node-based scalar fields from an XDMF file.

    Parameters
    ----------
    xdmf_field_path : str
        Path to the XDMF file containing point data.

    Returns
    -------
    points : (n_nodes, 3) ndarray of float
        Node coordinates (x, y, z).
    time_list : (n_steps,) ndarray of float
        Time values for each time step.
    scalar_field : (n_steps, n_nodes) ndarray of float
        Scalar values at nodes for each time step.

    Notes
    -----
    The function assumes a single scalar field exists in ``point_data`` at
    each time step and flattens it to 1D per step.
    """
    reader = ms.xdmf.TimeSeriesReader(xdmf_field_path)

    points, cells = reader.read_points_cells()
    n_nodes = points.shape[0]
    n_steps = reader.num_steps

    scalar_field = np.zeros((n_steps, n_nodes))
    time_list = np.zeros(n_steps)

    for k in range(reader.num_steps):
        # Read data
        time, point_data, cell_data = reader.read_data(k)

        # Add time
        time_list[k] = time

        # Add scalar
        field_name = list(point_data.keys())[0]
        scalar_field[k] = point_data[field_name].flatten()

    return points, time_list, scalar_field


def read_node_vector(xdmf_field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a time series of node-based 3D vector fields from an XDMF file.

    Parameters
    ----------
    xdmf_field_path : str
        Path to the XDMF file containing point data.

    Returns
    -------
    points : (n_nodes, 3) ndarray of float
        Node coordinates (x, y, z).
    time_list : (n_steps,) ndarray of float
        Time values for each time step.
    vector_field : (n_steps, n_nodes, 3) ndarray of float
        Vector values (vx, vy, vz) at nodes for each time step.

    Notes
    -----
    The function assumes a single vector field exists in ``point_data`` at
    each time step with shape ``(n_nodes, 3)``.
    """
    reader = ms.xdmf.TimeSeriesReader(xdmf_field_path)

    points, cells = reader.read_points_cells()
    n_nodes = points.shape[0]
    n_steps = reader.num_steps

    vector_field = np.zeros((n_steps, n_nodes, 3))
    time_list = np.zeros(n_steps)

    for k in range(reader.num_steps):
        # Read data
        time, point_data, cell_data = reader.read_data(k)

        # Add time
        time_list[k] = time

        # Add scalar
        field_name = list(point_data.keys())[0]
        vector_field[k,:,:] = point_data[field_name]

    return points, time_list, vector_field