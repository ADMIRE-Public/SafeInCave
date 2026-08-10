# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING
import shutil
import os

from .xdmf import XdmfWriter
from .vtkhdf import VtkhdfWriter

if TYPE_CHECKING:
    from ..Equations.Momentum import LinearMomentum, LinearMomentumBase
    from ..Equations.Heat import HeatDiffusion

    EqType = LinearMomentum | HeatDiffusion


# Available output formats. Each backend writes one file per registered field and
# exposes the same interface (see safeincave.Output.xdmf for the protocol), so a
# new format only needs a module here and an entry in this registry.
WRITER_BACKENDS = {
    XdmfWriter.extension: XdmfWriter,
    VtkhdfWriter.extension: VtkhdfWriter,
}


# Module-level shared storage for merged output files to prevent HDF5 contention.
# When multiple SaveFields instances write to the same output folder, they share
# a single writer object for the merged solution file instead of each creating
# their own (which causes HDF5 file truncation conflicts).
#
# All three are keyed by (normalized output folder path, output format), so that
# instances writing different formats into one folder do not share a writer.
#
# _shared_merged_outputs: dict[key, writer]
#   Maps key → writer backend object for the merged solution
#
# _shared_merged_output_refs: dict[key, int]
#   Reference count: how many SaveFields instances are using each merged output
#
# _shared_merged_output_field_names: dict[key, bool]
#   Tracks whether the duplicate-label mapping has already been built for a key.
#   Used to automatically rename duplicate field names by appending suffixes.
_shared_merged_outputs = {}
_shared_merged_output_refs = {}
_shared_merged_output_field_names = {}

class SaveFields:
    """
    Manage writing FEniCSx fields to disk over time.

    This helper collects references to fields stored on an equation object
    (either :class:`LinearMomentum` or :class:`HeatDiffusion`), opens one
    writer per field, and writes time-stamped data during a simulation.
    It can also copy the original Gmsh mesh file into the output directory
    for provenance.

    The file format is selected with ``output_format``; the writing code itself
    is format agnostic and lives in one module per format
    (:mod:`safeincave.Output.xdmf`, :mod:`safeincave.Output.vtkhdf`).

    Parameters
    ----------
    eq : EqType
        Equation/model object. Must expose:
        - ``eq.grid.mesh`` (a DOLFINx mesh with communicator),
        - ``eq.grid.grid_folder`` (path where the original ``.msh`` lives),
        - ``eq.grid.geometry_name`` (base filename of the ``.msh``),
        - attributes for each field you register via :meth:`add_output_field`.
    output_format : str, optional
        Output file format, one of the keys of :data:`WRITER_BACKENDS`
        (``"xdmf"``, the default, or ``"vtkhdf"``).

    Attributes
    ----------
    eq : EqType
        Stored equation/model handle.
    fields_data : list of dict
        Registered field descriptors, each with keys
        ``{"field_name": str, "label_name": str}``.
    output_fields : list
        Open writer backend objects, in the same order as ``fields_data``.
    merged_solutions : bool
        Whether to create and write a merged solution file. Set by Simulator.
    output_folder : str
        Base directory for outputs (set via :meth:`set_output_folder`).
    output_format : str
        Selected output file format.

    Raises
    ------
    ValueError
        If ``output_format`` is not a known format.

    Notes
    -----
    - The ``merged_solutions`` flag is typically set by the Simulator class
      so that all SaveFields instances share a consistent configuration.
    - Individual field files are always created regardless of merged_solutions setting.
    - Voigt/tensor conventions, function ranks, and meshtags are not managed
      here; this class only writes whatever :mod:`dolfinx` ``Function`` you
      provide in ``eq``.
      created by :meth:`initialize`. Ensure they exist beforehand.
    """

    def __init__(self, eq: LinearMomentumBase | HeatDiffusion, output_format: str = "xdmf"):
        if output_format not in WRITER_BACKENDS:
            raise ValueError(
                f"Unknown output format '{output_format}'. "
                f"Available formats: {', '.join(sorted(WRITER_BACKENDS))}."
            )
        self.eq = eq
        self.output_format = output_format
        self._writer_cls = WRITER_BACKENDS[output_format]
        self.fields_data = []
        self.output_fields = []
        self.merged_output = None
        self.merged_solutions = False  # Set by Simulator
        self.smooth_output = False  # Set by Simulator
        self.output_folder = None  # Set by set_output_folder()
        self._merged_field_name_map = {}  # Maps original label → unique name for merged output

    def set_output_folder(self, output_folder: str) -> None:
        """
        Set the base directory for all outputs.

        Parameters
        ----------
        output_folder : str
            Path to the directory where subfolders and output files will be placed.

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
            (appears in ParaView).

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
        Open one writer per registered field and write the mesh.

        For each entry in :attr:`fields_data`, opens a file at
        ``{output_folder}/{field_name}/{field_name}.{ext}`` (where ``ext`` is the
        extension of the selected :attr:`output_format`) and writes
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
            If the per-field output directory cannot be created.

        Notes
        -----
        - Files are opened in ``"w"`` mode (overwrite).
        """
        extension = self._writer_cls.extension
        for field_data in self.fields_data:
            field_name = field_data["field_name"]
            output_field = self._writer_cls(
                self.eq.grid.mesh.comm,
                os.path.join(self.output_folder, field_name, f"{field_name}.{extension}"),
                self.eq.grid.mesh,
            )
            self.output_fields.append(output_field)

        # Conditionally create merged solution file containing all fields
        # This behavior is controlled by the merged_solutions parameter
        if self.merged_solutions:
            # WHEN TO USE THE MERGED FILE vs THE INDIVIDUAL FILES
            #
            #    USE MERGED (solution.*) FOR:
            #    - Post-processing analysis across all fields simultaneously
            #    - Correlation studies between fields
            #    - Statistical analysis and data export
            #    - Quick initial result verification
            #
            #    USE INDIVIDUAL FILES FOR:
            #    - Mesh deformation visualization (Warp by Vector) → u/u.*
            #    - Field-specific analysis and filtering
            #    - Publication-quality visualizations
            #    - Large simulations (smaller file size per field)
            #
            # KNOWN XDMF-SPECIFIC LIMITATIONS OF THE MERGED FILE
            #
            # 1. VTK COMPOSITE METADATA (Artifacts)
            #    When multiple <Attribute> elements exist in the same XDMF <Grid>,
            #    VTK interprets this as a composite dataset and AUTOMATICALLY ADDS
            #    vtkCompositeIndex and vtkBlockColors. These clutter the ParaView
            #    field list but are harmless.
            #
            # 2. WARP BY VECTOR FILTER FAILURE
            #    ParaView's "Warp by Vector" filter (for mesh deformation visualization)
            #    FAILS or BEHAVES INCORRECTLY when multiple vector fields are present,
            #    because filter input selection is ambiguous with composite metadata.
            #
            # 3. DUPLICATED MESH
            #    DOLFINx's XDMFFile does not support external mesh references, so the
            #    mesh geometry/topology is duplicated in this file.
            #
            # The VTKHDF backend has none of these three limitations.
            #
            # Use shared merged output to prevent HDF5 file contention
            # when multiple SaveFields instances write to the same output folder
            merged_path = os.path.join(
                self.output_folder, "solution", f"solution.{extension}"
            )

            # Normalize path to handle symlinks and relative paths consistently.
            # The format is part of the key so that instances writing different
            # formats into one folder get one merged file each.
            merged_key = (
                os.path.normpath(os.path.abspath(self.output_folder)),
                self.output_format,
            )

            if merged_key not in _shared_merged_outputs:
                # First SaveFields instance for this output folder: create merged file
                _shared_merged_outputs[merged_key] = self._writer_cls(
                    self.eq.grid.mesh.comm,
                    merged_path,
                    self.eq.grid.mesh,
                )
                _shared_merged_output_refs[merged_key] = 1
            else:
                # Reuse existing merged output file and increment reference count
                _shared_merged_output_refs[merged_key] += 1

            # Store reference to shared merged output (may be created by this instance or reused)
            self.merged_output = _shared_merged_outputs[merged_key]

            # Build field name mapping to handle duplicate field names in merged output.
            # Detect duplicates and rename them with node/element suffixes to avoid HDF5 errors.
            if merged_key not in _shared_merged_output_field_names:
                # Extract all field labels and detect field types (node/element)
                name_to_indices_and_types = {}
                for i, field_data in enumerate(self.fields_data):
                    label_name = field_data["label_name"]
                    field_name = field_data["field_name"]
                    field = getattr(self.eq, field_name)
                    
                    # Detect if field is element-based or node-based
                    # Element-based: DG(0) or P(0) - degree 0
                    # Node-based: P(1), CG(1), etc. - degree >= 1
                    try:
                        ufl_elem = field.function_space.ufl_element()
                        degree = ufl_elem.degree
                        field_type = "(element)" if degree == 0 else "(node)"
                    except Exception:
                        field_type = None  # Unable to determine
                    
                    if label_name not in name_to_indices_and_types:
                        name_to_indices_and_types[label_name] = []
                    name_to_indices_and_types[label_name].append((i, field_type))
                
                # Build mapping from field index to unique name
                for label_name, occurrences in name_to_indices_and_types.items():
                    if len(occurrences) == 1:
                        # Single field with this label, use as-is
                        field_index, _ = occurrences[0]
                        self._merged_field_name_map[field_index] = label_name
                    else:
                        # Multiple fields with same label
                        field_types = [ft for _, ft in occurrences]
                        
                        if len(set(field_types)) > 1:
                            # Different field types, use field type as suffix
                            for field_index, field_type in occurrences:
                                if field_type:
                                    self._merged_field_name_map[field_index] = f"{label_name} {field_type}"
                                else:
                                    self._merged_field_name_map[field_index] = label_name
                        else:
                            # Same field type, use numeric suffix as fallback
                            for idx, (field_index, _) in enumerate(occurrences):
                                if idx == 0:
                                    self._merged_field_name_map[field_index] = label_name
                                else:
                                    self._merged_field_name_map[field_index] = f"{label_name}_{idx}"
                
                _shared_merged_output_field_names[merged_key] = True  # Mark as processed

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
        3. Calls ``write_function(field, t)`` on the individual writer.
        4. Calls ``write_function(field, t)`` on the merged solution writer.
           If duplicate field names exist, automatically appends suffixes to unique-ify them.
        """
        for i, field_data in enumerate(self.fields_data):
            field = getattr(self.eq, field_data["field_name"])
            original_name = field_data["label_name"]

            if self.smooth_output:
                try:
                    ufl_elem = field.function_space.ufl_element()
                    if ufl_elem.degree == 0:
                        field = self.eq.grid.smooth_dg0_to_cg1(field)
                except Exception:
                    pass

            field.name = original_name

            self.output_fields[i].write_function(field, t)

            # Write to merged solution file only if enabled.
            # For specific field visualization (e.g., deformation via Warp by Vector),
            # use individual field files in ParaView instead.
            if self.merged_output is not None:
                # Use pre-computed unique name mapping to handle duplicates
                # Map by field index to preserve duplicates (labels used as keys don't work)
                unique_name = self._merged_field_name_map.get(i, original_name)
                field.name = unique_name
                self.merged_output.write_function(field, t)

    def close(self) -> None:
        """
        Close all open output files.

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
        # for this output folder (reference counting) and if merged solutions are enabled
        if self.output_folder is not None and self.merged_output is not None and self.merged_solutions:
            merged_key = (
                os.path.normpath(os.path.abspath(self.output_folder)),
                self.output_format,
            )
            if merged_key in _shared_merged_output_refs:
                _shared_merged_output_refs[merged_key] -= 1
                if _shared_merged_output_refs[merged_key] == 0:
                    # Last SaveFields instance for this folder: close merged output
                    self.merged_output.close()
                    del _shared_merged_outputs[merged_key]
                    del _shared_merged_output_refs[merged_key]
                    del _shared_merged_output_field_names[merged_key]
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
        
        # Only copy if mesh file exists (may not exist for programmatically generated meshes)
        if os.path.exists(mesh_origin_file):
            shutil.copy(mesh_origin_file, mesh_destination_folder)
