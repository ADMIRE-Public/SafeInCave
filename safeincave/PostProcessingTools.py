"""
Useful to read xdmf files and post-process results.
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

import meshio as ms
import pandas as pd
import numpy as np
import os

def find_point_mapping(original_points, new_points):
    tol = 1e-10
    point_mapping = np.empty(len(new_points), dtype=int)
    for i, node in enumerate(new_points):
        match = np.where(np.all(np.abs(original_points - node) < tol, axis=1))[0]
        if len(match) == 1:
            point_mapping[i] = match[0]
        else:
            raise ValueError(f"Node {i} not found uniquely in the new mesh.")
    return point_mapping

def find_mapping(msh_points, msh_cells, xdmf_file):
    with ms.xdmf.TimeSeriesReader(xdmf_file) as reader:
        points, cells = reader.read_points_cells()
        mesh = ms.Mesh(points=points, cells=cells)
        x = mesh.points[:,0]
        y = mesh.points[:,1]
        z = mesh.points[:,2]
        xdmf_points = pd.DataFrame({'x': x, 'y': y, 'z': z})
        p1 = mesh.cells["tetra"][:,0]
        p2 = mesh.cells["tetra"][:,1]
        p3 = mesh.cells["tetra"][:,2]
        p4 = mesh.cells["tetra"][:,3]
        xdmf_cells = pd.DataFrame({'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4})
    mapping = find_point_mapping(xdmf_points.values, msh_points.values)
    return mapping


def read_msh_as_pandas(file_name):
    msh = ms.read(file_name)
    # print(msh.cells)
    # for k in range(len(msh.cells)):
    #     m = msh.cells["tetra"][k].data.shape[1]
    #     if m == 4:
    #         break
    # df_points = pd.DataFrame(msh.points, columns=["x", "y", "z"])
    # df_cells = pd.DataFrame(msh.cells[k].data, columns=["p1", "p2", "p3", "p4"])
    df_points = pd.DataFrame(msh.points, columns=["x", "y", "z"])
    df_cells = pd.DataFrame(msh.cells["tetra"], columns=["p1", "p2", "p3", "p4"])
    return df_points, df_cells

def compute_cell_centroids(points, cells):
    n, _ = cells.shape
    x_mid = np.zeros(n)
    y_mid = np.zeros(n)
    z_mid = np.zeros(n)
    x = points[:,0]
    y = points[:,1]
    z = points[:,2]
    for i in range(n):
        x_mid[i] = np.average(x[cells[i]])
        y_mid[i] = np.average(y[cells[i]])
        z_mid[i] = np.average(z[cells[i]])
    df_mid = pd.DataFrame({'x': x_mid, 'y': y_mid, 'z': z_mid})
    return df_mid

def read_xdmf_as_pandas(file_name):
    with ms.xdmf.TimeSeriesReader(file_name) as reader:
        points, cells = reader.read_points_cells()
        mesh = ms.Mesh(points=points, cells=cells)
        x = mesh.points[:,0]
        y = mesh.points[:,1]
        z = mesh.points[:,2]
        df_points = pd.DataFrame({'x': x, 'y': y, 'z': z})
        p1 = mesh.cells["tetra"][:,0]
        p2 = mesh.cells["tetra"][:,1]
        p3 = mesh.cells["tetra"][:,2]
        p4 = mesh.cells["tetra"][:,3]
        df_cells = pd.DataFrame({'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4})
    return df_points, df_cells

def read_scalar_from_cells(file_name):
    with ms.xdmf.TimeSeriesReader(file_name) as reader:
        points, cells = reader.read_points_cells()
        n = cells["tetra"].data.shape[0]
        m = reader.num_steps
        A = np.zeros((n, m))
        time_list = []
        for k in range(reader.num_steps):
            time, _, cell_data = reader.read_data(k)
            time_list.append(time)
            field_name = list(cell_data["tetra"].keys())[0]
            A[:,k] = cell_data["tetra"][field_name].flatten()
        df_scalar = pd.DataFrame(A, columns=time_list)
    return df_scalar

def read_scalar_from_points(file_name, mapping):
    with ms.xdmf.TimeSeriesReader(file_name) as reader:
        points, cells = reader.read_points_cells()
        n = points.shape[0]
        m = reader.num_steps
        A = np.zeros((n, m))
        time_list = []
        for k in range(reader.num_steps):
            time, point_data, _ = reader.read_data(k)
            time_list.append(time)
            field_name = list(point_data.keys())[0]
            A[:,k] = point_data[field_name][:,0]
        df_scalar = pd.DataFrame(A[mapping], columns=time_list)
    return df_scalar

def read_vector_from_points(file_name, point_mapping):
    with ms.xdmf.TimeSeriesReader(file_name) as reader:
        points, cells = reader.read_points_cells()
        n = points.shape[0]
        m = reader.num_steps
        Ax = np.zeros((n, m))
        Ay = np.zeros((n, m))
        Az = np.zeros((n, m))
        time_list = []
        for k in range(reader.num_steps):
            time, point_data, _ = reader.read_data(k)
            time_list.append(time)
            field_name = list(point_data.keys())[0]
            Ax[:,k] = point_data[field_name][:,0]
            Ay[:,k] = point_data[field_name][:,1]
            Az[:,k] = point_data[field_name][:,2]
        df_ux = pd.DataFrame(Ax[point_mapping], columns=time_list)
        df_uy = pd.DataFrame(Ay[point_mapping], columns=time_list)
        df_uz = pd.DataFrame(Az[point_mapping], columns=time_list)
    return df_ux, df_uy, df_uz

# def read_vector_from_points(file_name):
#     with ms.xdmf.TimeSeriesReader(file_name) as reader:
#         points, cells = reader.read_points_cells()
#         n = points.shape[0]
#         m = reader.num_steps
#         Ax = np.zeros((n, m))
#         Ay = np.zeros((n, m))
#         Az = np.zeros((n, m))
#         time_list = []
#         for k in range(reader.num_steps):
#             time, point_data, _ = reader.read_data(k)
#             time_list.append(time)
#             field_name = list(point_data.keys())[0]
#             Ax[:,k] = point_data[field_name][:,0]
#             Ay[:,k] = point_data[field_name][:,1]
#             Az[:,k] = point_data[field_name][:,2]
#         df_ux = pd.DataFrame(Ax, columns=time_list)
#         df_uy = pd.DataFrame(Ay, columns=time_list)
#         df_uz = pd.DataFrame(Az, columns=time_list)
#     return df_ux, df_uy, df_uz


def read_tensor_from_cells(file_name):
    with ms.xdmf.TimeSeriesReader(file_name) as reader:
        points, cells = reader.read_points_cells()
        n = cells["tetra"].data.shape[0]
        m = reader.num_steps
        sxx = np.zeros((n, m))
        syy = np.zeros((n, m))
        szz = np.zeros((n, m))
        sxy = np.zeros((n, m))
        sxz = np.zeros((n, m))
        syz = np.zeros((n, m))
        time_list = []
        for k in range(reader.num_steps):
            time, _, cell_data = reader.read_data(k)
            time_list.append(time)
            # field_name = list(cell_data.keys())[0]
            field_name = list(cell_data["tetra"].keys())[0]
            sxx[:,k] = cell_data["tetra"][field_name][:,0]
            syy[:,k] = cell_data["tetra"][field_name][:,4]
            szz[:,k] = cell_data["tetra"][field_name][:,8]
            sxy[:,k] = cell_data["tetra"][field_name][:,1]
            sxz[:,k] = cell_data["tetra"][field_name][:,2]
            syz[:,k] = cell_data["tetra"][field_name][:,5]
        df_sxx = pd.DataFrame(sxx, columns=time_list)
        df_syy = pd.DataFrame(syy, columns=time_list)
        df_szz = pd.DataFrame(szz, columns=time_list)
        df_sxy = pd.DataFrame(sxy, columns=time_list)
        df_sxz = pd.DataFrame(sxz, columns=time_list)
        df_syz = pd.DataFrame(syz, columns=time_list)
    return df_sxx, df_syy, df_szz, df_sxy, df_sxz, df_syz




