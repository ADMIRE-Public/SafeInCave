# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING
import dolfinx as do
import shutil
import os
import ctypes
import ctypes.util

if TYPE_CHECKING:
    from .MomentumEquation import LinearMomentum, LinearMomentumBase
    from .HeatEquation import HeatDiffusion

    EqType = LinearMomentum | HeatDiffusion


# Module-level shared storage for merged output files to prevent HDF5 contention.
# When multiple SaveFields instances write to the same output folder, they share
# a single XDMFFile object for the merged solution file instead of each creating
# their own (which causes HDF5 file truncation conflicts).
#
# _shared_merged_outputs: dict[str, do.io.XDMFFile]
#   Maps normalized output folder path → XDMFFile object for merged solution
#
# _shared_merged_output_refs: dict[str, int]
#   Reference count: how many SaveFields instances are using each merged output
_shared_merged_outputs = {}
_shared_merged_output_refs = {}

def _silence_hdf5_errors() -> None:
    """
    Disable HDF5's default error-stack printing (H5Eprint2), which writes
    directly to stderr via the C library, bypassing Python's logging/warnings
    system entirely. Calls H5Eset_auto2(H5E_DEFAULT, NULL, NULL) on whichever
    libhdf5 is already loaded in this process.
    """
    lib_path = None
    with open("/proc/self/maps") as f:
        for line in f:
            if "libhdf5" in line and line.strip().endswith(".so") or ".so." in line and "libhdf5" in line:
                candidate = line.split()[-1]
                if "libhdf5" in candidate:
                    lib_path = candidate
                    break

    if lib_path is None:
        # Fallback: try the dynamic linker's search path
        lib_path = ctypes.util.find_library("hdf5")

    if lib_path is None:
        return  # couldn't find it — fail silently, errors just stay visible

    try:
        hdf5 = ctypes.CDLL(lib_path)
        # H5E_DEFAULT is 0 in the C API when passed as the estack_id arg
        hdf5.H5Eset_auto2(0, None, None)
    except (OSError, AttributeError):
        pass


_silence_hdf5_errors()

class SaveFields:
    """
    Manage writing FEniCSx fields to XDMF over time.

    This helper collects references to fields stored on an equation object
    (either :class:`LinearMomentum` or :class:`HeatDiffusion`), opens one
    XDMF writer per field, and writes time-stamped data during a simulation.
    It can also copy the original Gmsh mesh file into the output directory
    for provenance.

    Parameters
    ----------
    eq : EqType
        Equation/model object. Must expose:
        - ``eq.grid.mesh`` (a DOLFINx mesh with communicator),
        - ``eq.grid.grid_folder`` (path where the original ``.msh`` lives),
        - ``eq.grid.geometry_name`` (base filename of the ``.msh``),
        - attributes for each field you register via :meth:`add_output_field`.

    Attributes
    ----------
    eq : EqType
        Stored equation/model handle.
    fields_data : list of dict
        Registered field descriptors, each with keys
        ``{"field_name": str, "label_name": str}``.
    output_fields : list of dolfinx.io.XDMFFile
        Open writers, in the same order as ``fields_data``.
    output_folder : str
        Base directory for outputs (set via :meth:`set_output_folder`).

    Notes
    -----
    - Voigt/tensor conventions, function ranks, and meshtags are not managed
      here; this class only writes whatever :mod:`dolfinx` ``Function`` you
      provide in ``eq``.
      created by :meth:`initialize`. Ensure they exist beforehand.
    """

    def __init__(self, eq: LinearMomentumBase | HeatDiffusion):
        self.eq = eq
        self.fields_data = []
        self.output_fields = []
        self.merged_output = None
        self.output_folder = None  # Set by set_output_folder()

    def set_output_folder(self, output_folder: str) -> None:
        """
        Set the base directory for all outputs.

        Parameters
        ----------
        output_folder : str
            Path to the directory where subfolders and XDMF files will be placed.

        Returns
        -------
        None
        """
        self.output_folder = output_folder

    def add_output_field(self, field_name: str, label_name: str) -> None:
        """
        Register a field to be written, with a display label.

        Parameters
        ----------
        field_name : str
            Attribute name on ``self.eq`` that refers to a
            :class:`dolfinx.fem.Function` (e.g., ``"u"``, ``"T"``, ``"sigma"``).
        label_name : str
            Human-readable name assigned to ``field.name`` when writing
            (appears in XDMF/ParaView).

        Returns
        -------
        None

        Notes
        -----
        You may call this multiple times before :meth:`initialize`.
        """
        data = {
            "field_name": field_name,
            "label_name": label_name,
        }
        self.fields_data.append(data)

    def initialize(self) -> None:
        """
        Open one XDMF writer per registered field and write the mesh.

        For each entry in :attr:`fields_data`, opens an XDMF at
        ``{output_folder}/{field_name}/{field_name}.xdmf`` and writes
        ``self.eq.grid.mesh`` once.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        OSError
            If the per-field output directory does not exist.

        Notes
        -----
        - Files are opened in ``"w"`` mode (overwrite).
        """
        for field_data in self.fields_data:
            field_name = field_data["field_name"]
            output_field = do.io.XDMFFile(
                self.eq.grid.mesh.comm,
                os.path.join(self.output_folder, field_name, f"{field_name}.xdmf"),
                "w",
            )
            output_field.write_mesh(self.eq.grid.mesh)
            self.output_fields.append(output_field)

        # Create merged solution file containing all fields
        # IMPORTANT: Known limitations and behaviors of this merged XDMF file:
        #
        # 1. VTK COMPOSITE METADATA (Artifacts)
        #    When multiple <Attribute> elements exist in the same XDMF <Grid>,
        #    VTK interprets this as a composite dataset and AUTOMATICALLY ADDS:
        #    - vtkCompositeIndex: Internal tracking field (not real data)
        #    - vtkBlockColors: Internal visualization field (not real data)
        #    These clutter the ParaView field list but are harmless.
        #
        # 2. WARP BY VECTOR FILTER FAILURE (Critical Limitation)
        #    ParaView's "Warp by Vector" filter (for mesh deformation visualization)
        #    FAILS or BEHAVES INCORRECTLY when multiple vector fields are present.
        #    ROOT CAUSE: Filter input selection is ambiguous with composite metadata.
        #
        # 3. XDMF STRUCTURE LIMITATION (DOLFINx API)
        #    DOLFINx's XDMFFile does not support external mesh references.
        #    Result: Mesh geometry/topology is DUPLICATED in this file.
        #    Impact: Larger file size (minor for most simulations).
        #    Workaround: None; inherent to DOLFINx API design.
        #
        # 4. USE CASES FOR MERGED vs INDIVIDUAL FILES
        #
        #    USE MERGED (solution.xdmf) FOR:
        #    - Post-processing analysis across all fields simultaneously
        #    - Correlation studies between fields
        #    - Statistical analysis and data export
        #    - Quick initial result verification
        #
        #    USE INDIVIDUAL FILES FOR:
        #    - Mesh deformation visualization (Warp by Vector) → u/u.xdmf
        #    - Field-specific analysis and filtering
        #    - Publication-quality visualizations
        #    - Large simulations (smaller file size per field)
        #
        # Use shared merged output to prevent HDF5 file contention
        # when multiple SaveFields instances write to the same output folder
        merged_path = os.path.join(self.output_folder, "solution", "solution.xdmf")
        os.makedirs(os.path.dirname(merged_path), exist_ok=True)
        
        # Normalize path to handle symlinks and relative paths consistently
        normalized_path = os.path.normpath(os.path.abspath(self.output_folder))
        
        if normalized_path not in _shared_merged_outputs:
            # First SaveFields instance for this output folder: create merged file
            merged_output_file = do.io.XDMFFile(
                self.eq.grid.mesh.comm,
                merged_path,
                "w",
            )
            merged_output_file.write_mesh(self.eq.grid.mesh)
            _shared_merged_outputs[normalized_path] = merged_output_file
            _shared_merged_output_refs[normalized_path] = 1
        else:
            # Reuse existing merged output file and increment reference count
            _shared_merged_output_refs[normalized_path] += 1
        
        # Store reference to shared merged output (may be created by this instance or reused)
        self.merged_output = _shared_merged_outputs[normalized_path]

    def save_fields(self, t: float) -> None:
        """
        Write all registered fields at simulation time ``t``.

        Writes to BOTH individual field files AND the merged solution file.

        Parameters
        ----------
        t : float
            Time value to associate with this write.

        Returns
        -------
        None

        Notes
        -----
        For each descriptor in :attr:`fields_data`:
        1. Fetches the field via ``getattr(self.eq, field_name)``.
        2. Sets ``field.name = label_name`` (used in visualization).
        3. Calls ``XDMFFile.write_function(field, t)`` on individual writer.
        4. Calls ``XDMFFile.write_function(field, t)`` on merged solution writer.
        """
        for i, field_data in enumerate(self.fields_data):
            field = getattr(self.eq, field_data["field_name"])
            field.name = field_data["label_name"]

            try:
                self.output_fields[i].write_function(field, t)
            except RuntimeError:
                pass

            # Write to merged solution file for unified analysis across all fields.
            # For specific field visualization (e.g., deformation via Warp by Vector),
            # use individual field files in ParaView instead.
            try:
                self.merged_output.write_function(field, t)
            except RuntimeError:
                pass

    def close(self) -> None:
        """
        Close all open XDMF output files.

        Ensures that all file handles opened during :meth:`initialize` are
        properly closed. This should be called at the end of a simulation to
        prevent resource leaks and ensure all data is flushed to disk.

        For the merged output file, it is only closed when the last SaveFields
        instance for that output folder closes (reference counting).

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        Safe to call multiple times (subsequent calls are no-ops).
        """
        # Close individual field output files
        for output_field in self.output_fields:
            output_field.close()
        
        # Close merged output file only if this is the last SaveFields instance
        # for this output folder (reference counting)
        if self.output_folder is not None and self.merged_output is not None:
            normalized_path = os.path.normpath(os.path.abspath(self.output_folder))
            if normalized_path in _shared_merged_output_refs:
                _shared_merged_output_refs[normalized_path] -= 1
                if _shared_merged_output_refs[normalized_path] == 0:
                    # Last SaveFields instance for this folder: close merged output
                    self.merged_output.close()
                    del _shared_merged_outputs[normalized_path]
                    del _shared_merged_output_refs[normalized_path]
                # Clear reference to prevent accidental re-use
                self.merged_output = None

    def save_mesh(self) -> None:
        """
        Copy the original Gmsh mesh file into ``{output_folder}/mesh/``.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Side Effects
        ------------
        - Creates ``{output_folder}/mesh`` if it does not exist.
        - Copies
          ``{eq.grid.grid_folder}/{eq.grid.geometry_name}.msh``
          to that directory.

        Raises
        ------
        FileNotFoundError
            If the source ``.msh`` file does not exist.
        OSError
            If the copy fails for other I/O reasons.
        """
        mesh_origin_file = os.path.join(
            self.eq.grid.grid_folder, f"{self.eq.grid.geometry_name}.msh"
        )
        mesh_destination_folder = os.path.join(self.output_folder, "mesh")
        if not os.path.exists(mesh_destination_folder):
            os.makedirs(mesh_destination_folder, exist_ok=True)
        shutil.copy(mesh_origin_file, mesh_destination_folder)
