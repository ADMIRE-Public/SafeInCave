# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Bulk post-run reader for a merged multi-field XDMF file
(``{output_folder}/solution/solution.xdmf``, written when a ``results``
output has ``merged_solutions: true``).

A merged file is **not** one shared per-timestep Grid with several
Attributes (what ``MergedFieldReaders.py``'s ``field_name=`` selection and
``meshio.xdmf.TimeSeriesReader`` both assume) -- ``SaveFields`` actually
writes it as N sibling ``GridType="Collection"`` Grids, one per field, each
holding its own 51-Uniform-Grid time series (verified directly against a
real ``solution.xdmf``: ``Domain/Grid[@Name='mesh']`` plus one
``Domain/Grid[@Name=field][@GridType='Collection']`` per field).
``meshio``'s reader can only ever surface one such collection, so reading a
real merged file through it silently returns just one field -- this is why
``DataExtract._read_field`` already treats a merged-file read as
best-effort, falling back to the per-field file on any failure. This
module reads the XDMF's actual Grid/DataItem structure directly (stdlib
``xml.etree.ElementTree``) and pulls each field's arrays straight out of
the shared HDF5 file via ``h5py``, so every field is read correctly, in one
open of the mesh's underlying HDF5 file.

This module only exposes the plain time-series data (:meth:`Dataset.variable`,
:meth:`Variable.mean`). Cylindrical-region selection and derived quantities
(:meth:`Dataset.cylindrical_geometry`, :meth:`Variable.calculate`) live in the
``SafeInCave_extensions`` overlay of this same file -- prototyped there first,
per this project's convention, before graduating to core.
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from typing import Dict, Optional

import h5py
import numpy as np

from .Readers import compute_cell_centroids

#: Friendly `variable()` names onto a merged file's actual field keys (which
#: are the case's YAML `fields:` spellings themselves -- see `codegen.py`'s
#: `add_output_field(field_name, label)` call, `label` being the original
#: YAML spelling). Only the handful of names worth a shorthand; anything
#: else is looked up by its exact on-disk name.
_VARIABLE_SHORTHANDS = {
    "strain": "total_strain_tensors",
    "stress": "stress_tensors",
    "displacement": "displacements",
}

#: XDMF `AttributeType` -> trailing array shape per point.
_SHAPE_BY_ATTR_TYPE = {"Tensor": (3, 3), "Vector": (3,), "Scalar": ()}


def _resolve_hdf_reference(xdmf_dir: str, text: str):
    """`'solution.h5:/Function/sig/0'` -> `(abs path to solution.h5, '/Function/sig/0')`."""
    filename, h5path = text.strip().split(":", 1)
    return os.path.join(xdmf_dir, filename), h5path


def _parse_merged_xdmf(xdmf_path: str):
    """
    Parse a merged XDMF's actual Grid structure: one `mesh` Grid (topology +
    geometry) plus one `GridType="Collection"` Grid per field, each holding
    one `Grid` per timestep with its own `Time`/`Attribute`/`DataItem`.

    Returns
    -------
    (topology_ref, geometry_ref) : each an (hdf5_path, internal_path) pair.
    fields : list of dict
        `{"name", "center" ('Node'/'Cell'), "attr_type" ('Tensor'/'Vector'/
        'Scalar'), "steps": [(time, hdf5_path, internal_path), ...]}`.
    """
    xdmf_dir = os.path.dirname(os.path.abspath(xdmf_path))
    domain = ET.parse(xdmf_path).getroot().find("Domain")
    grids = domain.findall("Grid")
    mesh_grid = next(g for g in grids if g.attrib.get("Name") == "mesh")
    topology_ref = _resolve_hdf_reference(xdmf_dir, mesh_grid.find(".//Topology/DataItem").text)
    geometry_ref = _resolve_hdf_reference(xdmf_dir, mesh_grid.find(".//Geometry/DataItem").text)

    fields = []
    for grid in grids:
        if grid is mesh_grid:
            continue
        steps = grid.findall("Grid")
        attr0 = steps[0].find("Attribute")
        entries = [
            (
                float(step.find("Time").attrib["Value"]),
                *_resolve_hdf_reference(xdmf_dir, step.find("Attribute/DataItem").text),
            )
            for step in steps
        ]
        fields.append({
            "name": grid.attrib["Name"],
            "center": attr0.attrib.get("Center", "Node"),
            "attr_type": attr0.attrib.get("AttributeType", "Scalar"),
            "steps": entries,
        })
    return (topology_ref, geometry_ref), fields


def _read_merged_fields(xdmf_path: str):
    """
    Read every field in a merged multi-field XDMF file, each in one pass
    over its own timesteps.

    Returns
    -------
    centroids : (n_cells, 3) ndarray
    points : (n_nodes, 3) ndarray
    time_list : (n_steps,) ndarray
    cell_fields : dict[str, ndarray]
        `{field_name: (n_steps, n_cells, ...)}` for every `Center="Cell"` field.
    node_fields : dict[str, ndarray]
        `{field_name: (n_steps, n_nodes, ...)}` for every `Center="Node"` field.
    """
    (topology_ref, geometry_ref), fields = _parse_merged_xdmf(xdmf_path)

    open_files: Dict[str, h5py.File] = {}

    def _dataset(ref):
        path, h5path = ref
        f = open_files.get(path)
        if f is None:
            f = open_files[path] = h5py.File(path, "r")
        return f[h5path][()]

    try:
        topology = _dataset(topology_ref)
        points = _dataset(geometry_ref)
        n_nodes = points.shape[0]
        centroids = compute_cell_centroids(topology, points)
        n_cells = topology.shape[0]

        time_list: Optional[np.ndarray] = None
        cell_fields: Dict[str, np.ndarray] = {}
        node_fields: Dict[str, np.ndarray] = {}
        for field in fields:
            suffix = _SHAPE_BY_ATTR_TYPE[field["attr_type"]]
            n_points = n_cells if field["center"] == "Cell" else n_nodes
            n_steps = len(field["steps"])
            arr = np.zeros((n_steps, n_points) + suffix)
            times = np.zeros(n_steps)
            for k, (t, path, h5path) in enumerate(field["steps"]):
                times[k] = t
                arr[k] = _dataset((path, h5path)).reshape((n_points,) + suffix)
            if time_list is None:
                time_list = times
            (cell_fields if field["center"] == "Cell" else node_fields)[field["name"]] = arr
    finally:
        for f in open_files.values():
            f.close()

    return centroids, points, time_list, cell_fields, node_fields


class TimeVariable:
    """The dataset's time axis, returned by `variable("time")`. Has no
    spatial dimension, so `mean()`/`calculate()` don't apply to it."""

    name = "time"

    def __init__(self, time: np.ndarray):
        self.values = time


class Variable:
    """
    One field's time series, at either the full dataset's cells/nodes or a
    region's subset of them (see the `SafeInCave_extensions` overlay's
    `RegionView`).

    Attributes
    ----------
    name : str
        The on-disk field name (e.g. `'stress_tensors'`).
    location : {'cell', 'node'}
    coords : (n_points, 3) ndarray
    time : (n_steps,) ndarray
    values : ndarray
        `(n_steps, n_points, 3, 3)` for a tensor field, `(n_steps, n_points, 3)`
        for a vector field, `(n_steps, n_points)` for a scalar field.
    """

    def __init__(self, name: str, location: str, coords: np.ndarray, time: np.ndarray,
                 values: np.ndarray, center: Optional[tuple] = None):
        self.name = name
        self.location = location
        self.coords = coords
        self.time = time
        self.values = values
        # Set only when sourced from a RegionView (cylindrical_geometry());
        # needed by the extensions overlay's calculate("radial"/"tangential").
        self._center = center

    def mean(self) -> np.ndarray:
        """Spatial mean over this variable's points, keeping the time series
        (one value per timestep)."""
        return self.values.mean(axis=1)


class Dataset:
    """Every field read from one merged XDMF file, see `read()`."""

    def __init__(self, centroids: np.ndarray, points: np.ndarray, time: np.ndarray,
                 cell_fields: Dict[str, np.ndarray], node_fields: Dict[str, np.ndarray]):
        self.centroids = centroids
        self.points = points
        self.time = time
        self.cell_fields = cell_fields
        self.node_fields = node_fields

    def variable(self, name: str):
        """Look up one field by its on-disk name, or by a friendly shorthand
        (`_VARIABLE_SHORTHANDS`), or `"time"` for the time axis."""
        if name in ("time", "t"):
            return TimeVariable(self.time)
        key = _VARIABLE_SHORTHANDS.get(name, name)
        if key in self.cell_fields:
            return Variable(key, "cell", self.centroids, self.time, self.cell_fields[key])
        if key in self.node_fields:
            return Variable(key, "node", self.points, self.time, self.node_fields[key])
        available = sorted(set(self.cell_fields) | set(self.node_fields))
        raise KeyError(f"Unknown variable {name!r}. Available: {available}")


def read(xdmf_path: str) -> Dataset:
    """
    Read every field saved in a merged multi-field XDMF file
    (`{output_folder}/solution/solution.xdmf`) into a `Dataset`.

    Parameters
    ----------
    xdmf_path : str
        Path to the merged `solution.xdmf` file (requires
        `merged_solutions: true` on the `results` output that wrote it).

    Returns
    -------
    Dataset
    """
    centroids, points, time_list, cell_fields, node_fields = _read_merged_fields(xdmf_path)
    return Dataset(centroids, points, time_list, cell_fields, node_fields)
