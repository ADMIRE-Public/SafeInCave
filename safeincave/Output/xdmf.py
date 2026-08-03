# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
XDMF writer backend.

Writer backends are duck-typed. A backend class must provide:

- ``extension``: class attribute with the file extension (no dot).
- ``__init__(comm, filepath, mesh)``: create parent directories, open the file
  for writing and write the mesh once.
- ``write_function(field, t)``: append a :class:`dolfinx.fem.Function` at time
  ``t``. The field is written under ``field.name``.
- ``close()``: close the file. Must be safe to call more than once.

This surface matches :class:`dolfinx.io.XDMFFile`, which subclasses of
:class:`~safeincave.Output.SaveFields` (in SafeInCave_extensions) call directly.
"""

from __future__ import annotations
import dolfinx as do
import os


class XdmfWriter:
    """
    Write a time series to an XDMF/HDF5 file pair via :mod:`dolfinx`.

    Parameters
    ----------
    comm : mpi4py.MPI.Comm
        Communicator of the mesh. Parallel I/O is handled by DOLFINx.
    filepath : str
        Path of the ``.xdmf`` file to create (overwritten if it exists).
    mesh : dolfinx.mesh.Mesh
        Mesh written once at construction time.
    """

    extension = "xdmf"

    def __init__(self, comm, filepath: str, mesh):
        os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
        self._file = do.io.XDMFFile(comm, filepath, "w")
        self._file.write_mesh(mesh)

    def write_function(self, field, t: float) -> None:
        self._file.write_function(field, t)

    def close(self) -> None:
        if self._file is not None:
            self._file.close()
            self._file = None
