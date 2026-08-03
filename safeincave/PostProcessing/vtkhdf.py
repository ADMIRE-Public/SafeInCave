# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
VTKHDF reader backend.

Counterpart of :mod:`safeincave.Output.vtkhdf`. Reads transient
``UnstructuredGrid`` files with :mod:`h5py` and returns the same arrays as the
XDMF backend, so post-processing scripts only differ by the file extension they
point at.

See :mod:`safeincave.PostProcessing.xdmf` for the reader backend protocol.
"""

import h5py
import numpy as np

from .Readers import compute_cell_centroids


def _read_mesh(root) -> tuple[np.ndarray, np.ndarray]:
    points = np.asarray(root["Points"][:], dtype=float)
    connectivity = np.asarray(root["Connectivity"][:]).reshape(-1, 4)
    return points, connectivity


def _read_times(root) -> np.ndarray:
    """Time values of the series; a file without a Steps group is a single step."""
    if "Steps" in root and "Values" in root["Steps"]:
        return np.asarray(root["Steps"]["Values"][:], dtype=float)
    return np.zeros(1)


def _read_series(root, data_group: str, offsets_group: str, n_entities: int,
                 n_steps: int) -> np.ndarray:
    """
    Read the single field of ``data_group`` as one block per time step.

    Blocks are located with the ``Steps/*Offsets`` entry of the field, falling
    back to contiguous blocks for files written without offsets.
    """
    if data_group not in root or len(root[data_group]) == 0:
        raise ValueError(f"No {data_group} found in VTKHDF file.")

    group = root[data_group]
    field_name = list(group.keys())[0]
    data = group[field_name]

    offsets_path = f"Steps/{offsets_group}/{field_name}"
    if offsets_path in root:
        offsets = np.asarray(root[offsets_path][:], dtype=int)
    else:
        offsets = np.arange(n_steps, dtype=int) * n_entities

    return np.stack(
        [np.asarray(data[offset:offset + n_entities], dtype=float) for offset in offsets]
    )


def read_cell_tensor(field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(field_path, "r") as file:
        root = file["VTKHDF"]
        points, connectivity = _read_mesh(root)
        time_list = _read_times(root)
        n_cells = connectivity.shape[0]
        tensor_field = _read_series(
            root, "CellData", "CellDataOffsets", n_cells, time_list.size
        )
    centroids = compute_cell_centroids(connectivity, points)
    return centroids, time_list, tensor_field.reshape(time_list.size, n_cells, 3, 3)


def read_cell_scalar(field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(field_path, "r") as file:
        root = file["VTKHDF"]
        points, connectivity = _read_mesh(root)
        time_list = _read_times(root)
        n_cells = connectivity.shape[0]
        scalar_field = _read_series(
            root, "CellData", "CellDataOffsets", n_cells, time_list.size
        )
    centroids = compute_cell_centroids(connectivity, points)
    return centroids, time_list, scalar_field.reshape(time_list.size, n_cells)


def read_node_scalar(field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(field_path, "r") as file:
        root = file["VTKHDF"]
        points, _ = _read_mesh(root)
        time_list = _read_times(root)
        n_nodes = points.shape[0]
        scalar_field = _read_series(
            root, "PointData", "PointDataOffsets", n_nodes, time_list.size
        )
    return points, time_list, scalar_field.reshape(time_list.size, n_nodes)


def read_node_vector(field_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(field_path, "r") as file:
        root = file["VTKHDF"]
        points, _ = _read_mesh(root)
        time_list = _read_times(root)
        n_nodes = points.shape[0]
        vector_field = _read_series(
            root, "PointData", "PointDataOffsets", n_nodes, time_list.size
        )
    return points, time_list, vector_field.reshape(time_list.size, n_nodes, 3)
