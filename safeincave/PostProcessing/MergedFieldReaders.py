# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Field-selecting counterparts of the four readers in
``safeincave.PostProcessing.Readers`` (``read_cell_tensor``, ``read_cell_scalar``,
``read_node_scalar``, ``read_node_vector``), which each assume exactly one
field is present and grab whichever key comes first -- fine for the ordinary
one-field-per-file output convention, but wrong for a merged multi-field XDMF
file (e.g. ``{output_folder}/solution/solution.xdmf``), where the caller must
say which field it wants.

These deliberately reuse the ``Readers.py`` function names; they are only
ever imported through this module's own path, so there is no collision.
"""

from __future__ import annotations

import meshio as ms
import numpy as np

from .Readers import compute_cell_centroids


def read_cell_tensor(
    xdmf_field_path: str, field_name: str | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a time series of cell-centered 3x3 tensor fields from an XDMF file.

    Parameters
    ----------
    xdmf_field_path : str
        Path to the XDMF file containing cell data (``cells['tetra']``).
    field_name : str, optional
        Which key to read from ``cell_data['tetra']``. Default: use whichever
        single field is present -- appropriate for per-field files, which
        always contain exactly one field. Pass explicitly when reading a
        merged multi-field file.

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
    If `field_name` is None, assumes a single tensor field is present under
    ``cell_data['tetra']`` at each time step (its key is used blindly); if
    given, that exact key is required to be present, else `KeyError`.
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
        key = field_name if field_name is not None else list(cell_data["tetra"].keys())[0]
        if key not in cell_data["tetra"]:
            raise KeyError(
                f"Field '{key}' not found in {xdmf_field_path}; "
                f"available: {list(cell_data['tetra'].keys())}"
            )
        tensor_field[k, :, :] = cell_data["tetra"][key].reshape((n_cells, 3, 3))

    return centroids, time_list, tensor_field


def read_cell_scalar(
    xdmf_field_path: str, field_name: str | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a time series of cell-centered scalar fields from an XDMF file.

    Parameters
    ----------
    xdmf_field_path : str
        Path to the XDMF file containing cell data (``cells['tetra']``).
    field_name : str, optional
        Which key to read from ``cell_data['tetra']``. Default: use whichever
        single field is present -- appropriate for per-field files, which
        always contain exactly one field. Pass explicitly when reading a
        merged multi-field file.

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
    If `field_name` is None, assumes a single scalar field is present under
    ``cell_data['tetra']`` at each time step (its key is used blindly); if
    given, that exact key is required to be present, else `KeyError`.
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
        key = field_name if field_name is not None else list(cell_data["tetra"].keys())[0]
        if key not in cell_data["tetra"]:
            raise KeyError(
                f"Field '{key}' not found in {xdmf_field_path}; "
                f"available: {list(cell_data['tetra'].keys())}"
            )
        scalar_field[k] = cell_data["tetra"][key].flatten()

    return centroids, time_list, scalar_field


def read_node_scalar(
    xdmf_field_path: str, field_name: str | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a time series of node-based scalar fields from an XDMF file.

    Parameters
    ----------
    xdmf_field_path : str
        Path to the XDMF file containing point data.
    field_name : str, optional
        Which key to read from ``point_data``. Default: use whichever single
        field is present. Pass explicitly when reading a merged multi-field
        file.

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
    If `field_name` is None, assumes a single scalar field exists in
    ``point_data`` at each time step (its key is used blindly) and flattens
    it to 1D per step; if given, that exact key is required, else `KeyError`.
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
        key = field_name if field_name is not None else list(point_data.keys())[0]
        if key not in point_data:
            raise KeyError(
                f"Field '{key}' not found in {xdmf_field_path}; "
                f"available: {list(point_data.keys())}"
            )
        scalar_field[k] = point_data[key].flatten()

    return points, time_list, scalar_field


def read_node_vector(
    xdmf_field_path: str, field_name: str | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a time series of node-based 3D vector fields from an XDMF file.

    Parameters
    ----------
    xdmf_field_path : str
        Path to the XDMF file containing point data.
    field_name : str, optional
        Which key to read from ``point_data``. Default: use whichever single
        field is present. Pass explicitly when reading a merged multi-field
        file.

    Returns
    -------
    points : (n_nodes, 3) ndarray of float
        Node coordinates (x, y, z).
    time_list : (n_steps,) ndarray of float
        Time values for each time step.
    vector_field : (n_steps, n_nodes, 3) ndarray of float
        Vector values at nodes for each time step.

    Notes
    -----
    If `field_name` is None, assumes a single vector field exists in
    ``point_data`` at each time step with shape ``(n_nodes, 3)`` (its key is
    used blindly); if given, that exact key is required, else `KeyError`.
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

        # Add vector
        key = field_name if field_name is not None else list(point_data.keys())[0]
        if key not in point_data:
            raise KeyError(
                f"Field '{key}' not found in {xdmf_field_path}; "
                f"available: {list(point_data.keys())}"
            )
        vector_field[k, :, :] = point_data[key]

    return points, time_list, vector_field
