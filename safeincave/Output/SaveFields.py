# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING
import dolfinx as do
import shutil
import os

if TYPE_CHECKING:
    from ..Equations.Momentum import LinearMomentum, LinearMomentumBase
    from ..Equations.Heat import HeatDiffusion

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
#
# _shared_merged_output_field_names: dict[str, set[str]]
#   Tracks field names (labels) written to each merged output to detect duplicates.
#   Maps normalized output folder path → set of field names written.
#   Used to automatically rename duplicate field names by appending suffixes.
_shared_merged_outputs = {}
_shared_merged_output_refs = {}
_shared_merged_output_field_names = {}

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
    merged_solutions : bool
        Whether to create and write a merged solution file. Set by Simulator.
    output_folder : str
        Base directory for outputs (set via :meth:`set_output_folder`).

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

    def __init__(self, eq: LinearMomentumBase | HeatDiffusion):
        self.eq = eq
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

        # Conditionally create merged solution file containing all fields
        # This behavior is controlled by the merged_solutions parameter
        if self.merged_solutions:
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
            
            # Build field name mapping to handle duplicate field names in merged output.
            # Detect duplicates and rename them with node/element suffixes to avoid HDF5 errors.
            if normalized_path not in _shared_merged_output_field_names:
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
                
                _shared_merged_output_field_names[normalized_path] = True  # Mark as processed

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
        # for this output folder (reference counting) and if merged solutions are enabled
        if self.output_folder is not None and self.merged_output is not None and self.merged_solutions:
            normalized_path = os.path.normpath(os.path.abspath(self.output_folder))
            if normalized_path in _shared_merged_output_refs:
                _shared_merged_output_refs[normalized_path] -= 1
                if _shared_merged_output_refs[normalized_path] == 0:
                    # Last SaveFields instance for this folder: close merged output
                    self.merged_output.close()
                    del _shared_merged_outputs[normalized_path]
                    del _shared_merged_output_refs[normalized_path]
                    del _shared_merged_output_field_names[normalized_path]
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
