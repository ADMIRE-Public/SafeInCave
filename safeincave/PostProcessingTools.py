"""
Useful to read xdmf files and post-process results.
"""
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
import pandas as pd
import numpy as np
import os

DataFrameType = pd.core.frame.DataFrame

def find_point_mapping(original_points: np.ndarray, new_points: np.ndarray) -> np.ndarray:
    """
    Find an index mapping from `new_points` to `original_points`.

    Each point in `new_points` is matched to its unique equal (within a
    tolerance) point in `original_points`. The result `idx` satisfies
    `original_points[idx[i]] ≈ new_points[i]`.

    Parameters
    ----------
    original_points : numpy.ndarray, shape (N, 3)
        Reference coordinate array.
    new_points : numpy.ndarray, shape (N, 3)
        Coordinate array to be mapped into the reference indexing.

    Returns
    -------
    numpy.ndarray, shape (N,)
        Integer indices such that `original_points[indices]` reorders
        the reference points to align with `new_points`.

    Raises
    ------
    ValueError
        If a point in `new_points` is not found **exactly once** in
        `original_points` (within the tolerance).

    Notes
    -----
    A fixed absolute tolerance of ``1e-10`` is used per coordinate.
    """
    tol = 1e-10
    point_mapping = np.empty(len(new_points), dtype=int)
    for i, node in enumerate(new_points):
        match = np.where(np.all(np.abs(original_points - node) < tol, axis=1))[0]
        if len(match) == 1:
            point_mapping[i] = match[0]
        else:
            raise ValueError(f"Node {i} not found uniquely in the new mesh.")
    return point_mapping

def find_mapping(msh_points: DataFrameType, xdmf_file: str) -> np.ndarray:
    """
    Build a point index mapping between an XDMF mesh and points from a MSH file.

    Parameters
    ----------
    msh_points : pandas.DataFrame
        DataFrame with columns ``['x', 'y', 'z']`` describing the MSH point
        coordinates (row order defines the desired ordering).
    xdmf_file : str
        Path to an XDMF file. The first mesh (points/cells) is read.

    Returns
    -------
    numpy.ndarray, shape (N,)
        Mapping indices such that ``xdmf_points[indices]`` reorders the XDMF
        points to match the order in ``msh_points``.

    Notes
    -----
    Assumes a tetrahedral mesh (cell type ``'tetra'``) exists in the XDMF.
    """
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


def read_msh_as_pandas(file_name: str) -> tuple[DataFrameType, DataFrameType]:
    """
    Read a Gmsh ``.msh`` file into pandas DataFrames.

    Parameters
    ----------
    file_name : str
        Path to the MSH file.

    Returns
    -------
    (df_points, df_cells) : tuple of pandas.DataFrame
        - ``df_points`` with columns ``['x', 'y', 'z']``.
        - ``df_cells`` with columns ``['p1', 'p2', 'p3', 'p4']`` (tetra connectivity).

    Notes
    -----
    Expects tetrahedral cells under key ``'tetra'``.
    """
    msh = ms.read(file_name)
    df_points = pd.DataFrame(msh.points, columns=["x", "y", "z"])
    df_cells = pd.DataFrame(msh.cells["tetra"], columns=["p1", "p2", "p3", "p4"])
    return df_points, df_cells

def compute_cell_centroids(points: np.ndarray, cells: np.ndarray) -> DataFrameType:
    """
    Compute cell centroids as the arithmetic mean of vertex coordinates.

    Parameters
    ----------
    points : numpy.ndarray, shape (N, 3)
        Point coordinates.
    cells : numpy.ndarray, shape (M, 4)
        Tetrahedral cell connectivity (indices into ``points``).

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ``['x', 'y', 'z']`` containing centroids
        for each cell (row-wise).
    """
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

def read_xdmf_as_pandas(file_name: str) -> tuple[DataFrameType, DataFrameType]:
    """
    Read an XDMF mesh into pandas DataFrames (points and tetra cells).

    Parameters
    ----------
    file_name : str
        Path to the XDMF file. The first mesh (points/cells) is read.

    Returns
    -------
    (df_points, df_cells) : tuple of pandas.DataFrame
        - ``df_points`` with columns ``['x', 'y', 'z']``.
        - ``df_cells`` with columns ``['p1', 'p2', 'p3', 'p4']``.

    Notes
    -----
    Assumes tetrahedral cells under key ``'tetra'``.
    """
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

def read_scalar_from_cells(file_name: str) -> DataFrameType:
    """
    Read a cell-wise scalar field time series from an XDMF file.

    Parameters
    ----------
    file_name : str
        Path to the XDMF file containing cell data over multiple time steps.

    Returns
    -------
    pandas.DataFrame
        DataFrame of shape ``(n_cells, n_steps)`` where columns are time
        values (as read from the file) and rows correspond to cells.

    Notes
    -----
    - Assumes tetrahedral cells under key ``'tetra'``.
    - Uses the **first** available cell field at each time step.
    """
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

def read_scalar_from_points(file_name: str, mapping: np.ndarray) -> DataFrameType:
    """
    Read a nodal scalar field time series from an XDMF file and reorder rows.

    Parameters
    ----------
    file_name : str
        Path to the XDMF file containing point data over multiple time steps.
    mapping : numpy.ndarray, shape (n_points,)
        Index mapping to reorder the XDMF point order into a target order
        (e.g., MSH order). The result is ``A[mapping]``.

    Returns
    -------
    pandas.DataFrame
        DataFrame of shape ``(n_points, n_steps)`` with columns as time values.

    Notes
    -----
    Uses the **first** available point field at each time step and the first
    component (column 0) of that field.
    """
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

def read_vector_from_points(file_name: str, point_mapping: np.ndarray) -> tuple[DataFrameType, DataFrameType, DataFrameType]:
    """
    Read a nodal 3D vector field time series from an XDMF file and reorder rows.

    Parameters
    ----------
    file_name : str
        Path to the XDMF file containing point data over multiple time steps.
    point_mapping : numpy.ndarray, shape (n_points,)
        Index mapping to reorder the XDMF point order into a target order.

    Returns
    -------
    (df_ux, df_uy, df_uz) : tuple of pandas.DataFrame
        One DataFrame per component, each of shape ``(n_points, n_steps)``,
        with columns as time values.

    Notes
    -----
    Uses the **first** available point field at each time step and splits
    its three components into separate DataFrames.
    """
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


def read_tensor_from_cells(file_name: str) -> tuple[
        DataFrameType, DataFrameType, DataFrameType, DataFrameType, DataFrameType, DataFrameType
    ]:
    """
    Read a symmetric 3×3 tensor field (per cell) from an XDMF file over time.

    The tensor is split into six component DataFrames using the ordering:
    ``sxx, syy, szz, sxy, sxz, syz``.

    Parameters
    ----------
    file_name : str
        Path to the XDMF file containing cell data over multiple time steps.

    Returns
    -------
    tuple of pandas.DataFrame
        ``(df_sxx, df_syy, df_szz, df_sxy, df_sxz, df_syz)``, each of shape
        ``(n_cells, n_steps)`` with columns as time values.

    Notes
    -----
    - Assumes tetrahedral cells under key ``'tetra'``.
    - The component extraction uses the flattened 3×3 storage in row-major
      order: indices ``[0,4,8,1,2,5]`` map to
      ``(xx, yy, zz, xy, xz, yz)`` respectively.
    """
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




