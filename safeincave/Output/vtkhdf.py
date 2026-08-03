# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
VTKHDF writer backend.

Writes a transient ``UnstructuredGrid`` following the VTKHDF specification
(https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/), the
format Kitware maintains as the successor of the XDMF/HDF5 pair. Files are
written directly with :mod:`h5py`, so no dependency on ``dolfinx.io.vtkhdf``
(DOLFINx >= 0.10) is needed.

The mesh is static: geometry and topology are stored once, and every time step
appends a block of values to the field datasets. ``Steps/*Offsets`` point at the
block belonging to each step.

See :mod:`safeincave.Output.xdmf` for the writer backend protocol.

Limitations
-----------
- Serial runs only: distributed meshes are not gathered. Use ``"xdmf"`` for
  MPI-parallel runs.
- Tetrahedral meshes only (VTK cell type 10), matching the rest of SafeInCave.
- Degree 0 (cell data) and degree 1 (point data) functions only.
"""

from __future__ import annotations
import numpy as np
import h5py
import os

#: VTK cell type id for a linear tetrahedron.
_VTK_TETRA = 10

#: VTKHDF specification version this writer targets.
_VTKHDF_VERSION = (2, 4)

#: Target size of an HDF5 chunk, in doubles. Chunks holding one time step of a
#: field are ideal for both appending and reading, but must stay bounded for
#: large meshes, otherwise HDF5 reads and caches far more than a step at a time.
_CHUNK_DOUBLES = 128 * 1024


def _append(dset, values: np.ndarray) -> int:
    """Append rows to a resizable dataset and return the row offset used."""
    offset = dset.shape[0]
    dset.resize(offset + values.shape[0], axis=0)
    dset[offset:] = values
    return offset


def _chunk_shape(n_rows: int, row_shape: tuple) -> tuple:
    """Chunk covering one time step of a field, capped at :data:`_CHUNK_DOUBLES`."""
    row_size = int(np.prod(row_shape, dtype=int)) if row_shape else 1
    return (max(1, min(n_rows, _CHUNK_DOUBLES // row_size)),) + row_shape


class VtkhdfWriter:
    """
    Write a time series to a VTKHDF file.

    Parameters
    ----------
    comm : mpi4py.MPI.Comm
        Communicator of the mesh. Must be of size 1.
    filepath : str
        Path of the ``.vtkhdf`` file to create (overwritten if it exists).
    mesh : dolfinx.mesh.Mesh
        Tetrahedral mesh, written once at construction time.

    Raises
    ------
    RuntimeError
        If ``comm`` has more than one rank.
    ValueError
        If the mesh is not tetrahedral.
    """

    extension = "vtkhdf"

    def __init__(self, comm, filepath: str, mesh):
        if comm is not None and comm.size > 1:
            raise RuntimeError(
                "VTKHDF output supports serial runs only. "
                "Use output_format='xdmf' for MPI-parallel runs."
            )

        points = np.ascontiguousarray(mesh.geometry.x, dtype=np.float64)
        connectivity = np.asarray(mesh.geometry.dofmap)
        if connectivity.ndim != 2 or connectivity.shape[1] != 4:
            raise ValueError(
                "VTKHDF output supports tetrahedral meshes only, but the mesh has "
                f"{connectivity.shape[1] if connectivity.ndim == 2 else '?'} nodes per cell."
            )

        self.n_points = points.shape[0]
        self.n_cells = connectivity.shape[0]
        self._times: list[float] = []
        self._perm_cache: dict[int, tuple] = {}

        os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
        self._file = h5py.File(filepath, "w")
        root = self._file.create_group("VTKHDF")
        root.attrs.create("Version", np.array(_VTKHDF_VERSION, dtype=np.int64))
        # Must be a fixed-length ASCII string; VTK rejects h5py's default UTF-8 type.
        root.attrs.create("Type", np.bytes_("UnstructuredGrid"))

        root.create_dataset("NumberOfPoints", data=np.array([self.n_points], np.int64))
        root.create_dataset("NumberOfCells", data=np.array([self.n_cells], np.int64))
        root.create_dataset(
            "NumberOfConnectivityIds", data=np.array([4 * self.n_cells], np.int64)
        )
        root.create_dataset("Points", data=points)
        root.create_dataset(
            "Connectivity", data=connectivity.reshape(-1).astype(np.int64)
        )
        root.create_dataset(
            "Offsets", data=(np.arange(self.n_cells + 1, dtype=np.int64) * 4)
        )
        root.create_dataset("Types", data=np.full(self.n_cells, _VTK_TETRA, np.uint8))

        root.create_group("PointData")
        root.create_group("CellData")

        steps = root.create_group("Steps")
        steps.attrs.create("NSteps", np.int64(0))
        steps.create_dataset(
            "Values", (0,), maxshape=(None,), chunks=(1024,), dtype=np.float64
        )
        for name in ("PartOffsets", "NumberOfParts", "PointOffsets"):
            steps.create_dataset(
                name, (0,), maxshape=(None,), chunks=(1024,), dtype=np.int64
            )
        for name in ("CellOffsets", "ConnectivityIdOffsets"):
            steps.create_dataset(
                name, (0, 1), maxshape=(None, 1), chunks=(1024, 1), dtype=np.int64
            )
        steps.create_group("PointDataOffsets")
        steps.create_group("CellDataOffsets")

        self._root = root
        self._steps = steps

    def write_function(self, field, t: float) -> None:
        """
        Append ``field`` at time ``t``, under the name ``field.name``.

        A new time step is opened whenever ``t`` differs from the previous call,
        so several fields written at the same time (merged solution file) share
        one step.

        Raises
        ------
        ValueError
            If the function is neither degree 0 nor degree 1.
        """
        if self._file is None:
            raise ValueError("Cannot write to a closed VTKHDF file.")

        if not self._times or t != self._times[-1]:
            self._open_step(t)
        step = len(self._times) - 1

        degree = field.function_space.ufl_element().degree
        if degree == 0:
            values = self._cell_values(field)
            group, offsets_group, n_entities = (
                self._root["CellData"],
                self._steps["CellDataOffsets"],
                self.n_cells,
            )
        elif degree == 1:
            values = self._point_values(field)
            group, offsets_group, n_entities = (
                self._root["PointData"],
                self._steps["PointDataOffsets"],
                self.n_points,
            )
        else:
            raise ValueError(
                f"VTKHDF output supports degree 0 and degree 1 functions, got degree {degree} "
                f"for field '{field.name}'. Interpolate onto a P1 space before writing."
            )

        name = str(field.name).replace("/", "_")
        if name not in group:
            self._create_field(
                group, offsets_group, name, values.shape[1:], step, n_entities
            )

        dset, offsets = group[name], offsets_group[name]
        if offsets.shape[0] > step:
            # Same field written twice at the same time: overwrite its block.
            start = int(offsets[step])
            dset[start:start + n_entities] = values
        else:
            offset = _append(dset, values)
            _append(offsets, np.array([offset], np.int64))

        self._file.flush()

    def close(self) -> None:
        if self._file is not None:
            self._file.close()
            self._file = None

    # -- internals ---------------------------------------------------------

    def _open_step(self, t: float) -> None:
        """Record a new time step. The mesh is static, so its offsets stay zero."""
        self._times.append(t)
        n = len(self._times)
        _append(self._steps["Values"], np.array([t], np.float64))
        _append(self._steps["PartOffsets"], np.zeros(1, np.int64))
        _append(self._steps["NumberOfParts"], np.ones(1, np.int64))
        _append(self._steps["PointOffsets"], np.zeros(1, np.int64))
        _append(self._steps["CellOffsets"], np.zeros((1, 1), np.int64))
        _append(self._steps["ConnectivityIdOffsets"], np.zeros((1, 1), np.int64))
        self._steps.attrs.modify("NSteps", np.int64(n))

    def _create_field(self, group, offsets_group, name: str, shape, step: int,
                      n_entities: int) -> None:
        """
        Create the datasets for a field appearing for the first time.

        Every field needs one offset per time step. A field registered after the
        first step gets its earlier steps pointed at block 0, so the file stays
        readable instead of ragged.
        """
        group.create_dataset(
            name,
            (0,) + shape,
            maxshape=(None,) + shape,
            chunks=_chunk_shape(n_entities, shape),
            dtype=np.float64,
        )
        offsets_group.create_dataset(
            name, (0,), maxshape=(None,), chunks=(1024,), dtype=np.int64
        )
        if step > 0:
            _append(offsets_group[name], np.zeros(step, np.int64))

    def _cell_values(self, field) -> np.ndarray:
        values = np.asarray(field.x.array, dtype=np.float64)
        return self._shape_components(values, self.n_cells)

    def _point_values(self, field) -> np.ndarray:
        """
        Reorder degree-1 dof values to match the geometry node ordering.

        Point data must follow the ``Points`` dataset, but the dof numbering of a
        P1 space is independent of the geometry node numbering.
        """
        space = field.function_space
        cached = self._perm_cache.get(id(space))
        if cached is None:
            geometry_nodes = np.asarray(space.mesh.geometry.dofmap)
            dofs = np.asarray(space.dofmap.list)
            perm = np.empty(self.n_points, dtype=np.int64)
            perm[geometry_nodes] = dofs
            # Keep the space alive so its id() cannot be reused by another object.
            self._perm_cache[id(space)] = (space, perm)
        else:
            perm = cached[1]

        n_dofs = space.dofmap.index_map.size_local + space.dofmap.index_map.num_ghosts
        values = np.asarray(field.x.array, dtype=np.float64).reshape(n_dofs, -1)[perm]
        return self._shape_components(values.reshape(-1), self.n_points)

    @staticmethod
    def _shape_components(values: np.ndarray, n_entities: int) -> np.ndarray:
        """Shape a flat dof array as (n_entities,) for scalars, (n_entities, k) otherwise."""
        n_components = values.size // n_entities
        if n_components == 1:
            return values.reshape(n_entities)
        return values.reshape(n_entities, n_components)
