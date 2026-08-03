# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
XDMF reader backend.

Reader backends are duck-typed modules providing ``read_cell_tensor``,
``read_cell_scalar``, ``read_node_scalar`` and ``read_node_vector``, each taking
a file path and returning ``(coordinates, time_list, field)`` as documented in
:mod:`safeincave.PostProcessing.Readers`, which dispatches to them by file
extension. Every reader assumes the file holds a single field, as written by
:class:`~safeincave.Output.SaveFields`.
"""

import meshio as ms
import numpy as np

from .Readers import compute_cell_centroids


def read_cell_tensor(field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reader = ms.xdmf.TimeSeriesReader(field_path)
    points, cells = reader.read_points_cells()
    n_cells = cells["tetra"].shape[0]
    n_steps = reader.num_steps

    centroids = compute_cell_centroids(cells["tetra"], points)
    tensor_field = np.zeros((n_steps, n_cells, 3, 3))
    time_list = np.zeros(n_steps)

    for k in range(reader.num_steps):
        time, point_data, cell_data = reader.read_data(k)
        time_list[k] = time
        field_name = list(cell_data["tetra"].keys())[0]
        tensor_field[k, :, :] = cell_data["tetra"][field_name].reshape((n_cells, 3, 3))

    return centroids, time_list, tensor_field


def read_cell_scalar(field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reader = ms.xdmf.TimeSeriesReader(field_path)

    points, cells = reader.read_points_cells()
    n_cells = cells["tetra"].shape[0]
    n_steps = reader.num_steps

    centroids = compute_cell_centroids(cells["tetra"], points)
    scalar_field = np.zeros((n_steps, n_cells))
    time_list = np.zeros(n_steps)

    for k in range(reader.num_steps):
        time, point_data, cell_data = reader.read_data(k)
        time_list[k] = time
        field_name = list(cell_data["tetra"].keys())[0]
        scalar_field[k] = cell_data["tetra"][field_name].flatten()

    return centroids, time_list, scalar_field


def read_node_scalar(field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reader = ms.xdmf.TimeSeriesReader(field_path)

    points, cells = reader.read_points_cells()
    n_nodes = points.shape[0]
    n_steps = reader.num_steps

    scalar_field = np.zeros((n_steps, n_nodes))
    time_list = np.zeros(n_steps)

    for k in range(reader.num_steps):
        time, point_data, cell_data = reader.read_data(k)
        time_list[k] = time
        field_name = list(point_data.keys())[0]
        scalar_field[k] = point_data[field_name].flatten()

    return points, time_list, scalar_field


def read_node_vector(field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reader = ms.xdmf.TimeSeriesReader(field_path)

    points, cells = reader.read_points_cells()
    n_nodes = points.shape[0]
    n_steps = reader.num_steps

    vector_field = np.zeros((n_steps, n_nodes, 3))
    time_list = np.zeros(n_steps)

    for k in range(reader.num_steps):
        time, point_data, cell_data = reader.read_data(k)
        time_list[k] = time
        field_name = list(point_data.keys())[0]
        vector_field[k, :, :] = point_data[field_name]

    return points, time_list, vector_field
