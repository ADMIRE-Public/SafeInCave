# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Vector/tensor XDMF readers that :mod:`safeincave.Output.DataExtract`
needs beyond what ``safeincave.PostProcessing.Readers`` provides
(``read_cell_tensor``, ``read_cell_scalar``, ``read_node_scalar``,
``read_node_vector`` -- covering neither a cell-centered vector field like
``principal_stresses`` nor a node-centered tensor field like a smoothed
``sig``).

The sibling module ``MergedFieldReaders.py`` complements these by adding
``field_name=`` selection on top of the four ``Readers.py`` functions, for
reading a merged multi-field file.
"""

from __future__ import annotations

import meshio as ms
import numpy as np

from .Readers import compute_cell_centroids


def read_cell_vector(
    xdmf_field_path: str, field_name: str | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a time series of cell-centered 3-component vector fields from an
    XDMF file (e.g. DG0 principal-value fields: ``principal_stresses``,
    ``principal_strains``).

    Parameters
    ----------
    xdmf_field_path : str
        Path to the XDMF file containing cell data (``cells['tetra']``).
    field_name : str, optional
        Which key to read from ``cell_data['tetra']``. Default: use whichever
        single field is present (today's behavior) -- appropriate for
        per-field files, which always contain exactly one field. Pass
        explicitly when reading a merged multi-field file.

    Returns
    -------
    centroids : (n_cells, 3) ndarray of float
        Centroid coordinates of the tetrahedral cells.
    time_list : (n_steps,) ndarray of float
        Time values for each time step.
    vector_field : (n_steps, n_cells, 3) ndarray of float
        Vector values per time step and cell.

    Notes
    -----
    If `field_name` is None, assumes a single vector field is present under
    ``cell_data['tetra']`` at each time step (its key is used blindly); if
    given, that exact key is required to be present, else `KeyError`.
    """
    reader = ms.xdmf.TimeSeriesReader(xdmf_field_path)

    points, cells = reader.read_points_cells()
    n_cells = cells["tetra"].shape[0]
    n_steps = reader.num_steps

    centroids = compute_cell_centroids(cells["tetra"], points)
    vector_field = np.zeros((n_steps, n_cells, 3))
    time_list = np.zeros(n_steps)

    for k in range(reader.num_steps):
        # Read data
        time, point_data, cell_data = reader.read_data(k)

        # Add time
        time_list[k] = time

        # Add vector
        key = field_name if field_name is not None else list(cell_data["tetra"].keys())[0]
        if key not in cell_data["tetra"]:
            raise KeyError(
                f"Field '{key}' not found in {xdmf_field_path}; "
                f"available: {list(cell_data['tetra'].keys())}"
            )
        vector_field[k, :, :] = cell_data["tetra"][key].reshape((n_cells, 3))

    return centroids, time_list, vector_field


def read_node_tensor(
    xdmf_field_path: str, field_name: str | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a time series of node-based 3x3 tensor fields from an XDMF file.

    Some output fields (e.g. ``SaveFieldsNewton`` with ``smooth_output=True``)
    are written node-centered rather than cell-centered, even though the
    underlying quantity is otherwise identical to a ``read_cell_tensor``
    field -- this is the node-data counterpart.

    Parameters
    ----------
    xdmf_field_path : str
        Path to the XDMF file containing point data.
    field_name : str, optional
        Which key to read from ``point_data``. Default: use whichever single
        field is present (today's behavior). Pass explicitly when reading a
        merged multi-field file.

    Returns
    -------
    points : (n_nodes, 3) ndarray of float
        Node coordinates (x, y, z).
    time_list : (n_steps,) ndarray of float
        Time values for each time step.
    tensor_field : (n_steps, n_nodes, 3, 3) ndarray of float
        Tensor values at nodes for each time step.

    Notes
    -----
    If `field_name` is None, assumes a single tensor field exists in
    ``point_data`` at each time step (its key is used blindly); if given,
    that exact key is required, else `KeyError`.
    """
    reader = ms.xdmf.TimeSeriesReader(xdmf_field_path)

    points, cells = reader.read_points_cells()
    n_nodes = points.shape[0]
    n_steps = reader.num_steps

    tensor_field = np.zeros((n_steps, n_nodes, 3, 3))
    time_list = np.zeros(n_steps)

    for k in range(reader.num_steps):
        # Read data
        time, point_data, cell_data = reader.read_data(k)

        # Add time
        time_list[k] = time

        # Add tensor
        key = field_name if field_name is not None else list(point_data.keys())[0]
        if key not in point_data:
            raise KeyError(
                f"Field '{key}' not found in {xdmf_field_path}; "
                f"available: {list(point_data.keys())}"
            )
        tensor_field[k, :, :] = point_data[key].reshape((n_nodes, 3, 3))

    return points, time_list, tensor_field
