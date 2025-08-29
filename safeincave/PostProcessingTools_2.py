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


def build_smoother(points, conn):
    def tetrahedron_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
        volume = abs((1/6) * ((x2 - x1) * ((y3 - y1)*(z4 - z1) - (z3 - z1)*(y4 - y1)) + 
                     (y2 - y1) * ((z3 - z1)*(x4 - x1) - (x3 - x1)*(z4 - z1)) + 
                     (z2 - z1) * ((x3 - x1)*(y4 - y1) - (y3 - y1)*(x4 - x1))))
        return volume

    def compute_volumes():
        conn = mesh.topology.connectivity(3, 0).array.reshape(-1, 4)
        points = mesh.geometry.x
        volumes = np.zeros(n_elems)
        for i in range(n_elems):
            nodes = conn[i]
            x1, y1, z1 = coord[nodes[0], 0], coord[nodes[0], 1], coord[nodes[0], 2]
            x2, y2, z2 = coord[nodes[1], 0], coord[nodes[1], 1], coord[nodes[1], 2]
            x3, y3, z3 = coord[nodes[2], 0], coord[nodes[2], 1], coord[nodes[2], 2]
            x4, y4, z4 = coord[nodes[3], 0], coord[nodes[3], 1], coord[nodes[3], 2]
            volumes[i] = tetrahedron_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)

    def build_node_elem_stencil(conn, coord):
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

def build_mapping(nodes_xdmf, nodes_msh):
    return [np.where((nodes_msh == row).all(axis=1))[0][0] for row in nodes_xdmf]

def find_closest_point(target_point, points) -> int:
    x_p, y_p, z_p = target_point
    d = np.sqrt(  (points[:,0] - x_p)**2
                + (points[:,1] - y_p)**2
                + (points[:,2] - z_p)**2 )
    p_idx = d.argmin()
    return p_idx


def compute_cell_centroids(cells, points):
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


def read_cell_tensor(xdmf_field_path):
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


def read_cell_scalar(xdmf_field_path):
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


def read_node_scalar(xdmf_field_path):
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


def read_node_vector(xdmf_field_path):
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