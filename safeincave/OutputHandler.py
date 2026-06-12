# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import TYPE_CHECKING
import dolfinx as do
import shutil
import os

if TYPE_CHECKING:
    from .MomentumEquation import LinearMomentum, LinearMomentumBase
    from .HeatEquation import HeatDiffusion

    EqType = LinearMomentum | HeatDiffusion


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
        #    These clutter the ParaView field list but are harmless. To hide them:
        #    Right-click field name in ParaView → Hide
        #
        # 2. WRAP BY VECTOR FILTER FAILURE (Critical Limitation)
        #    ParaView's "Wrap by Vector" filter (for mesh deformation visualization)
        #    FAILS or BEHAVES INCORRECTLY when multiple vector fields are present.
        #    ROOT CAUSE: Filter input selection is ambiguous with composite metadata.
        #    WORKAROUND: Use individual field files instead:
        #    - output_folder/u/u.xdmf for displacement deformation
        #    Do NOT use this merged file for deformation visualization.
        #
        # 3. XDMF STRUCTURE LIMITATION (DOLFINx API)
        #    DOLFINx's XDMFFile does not support external mesh references.
        #    Result: Mesh geometry/topology is DUPLICATED in this file.
        #    Impact: Larger file size (minor for most simulations).
        #    Workaround: None; inherent to DOLFINx API design.
        #
        # 4. USE CASES FOR MERGED vs INDIVIDUAL FILES
        #    USE MERGED (solution.xdmf) FOR:
        #    - Post-processing analysis across all fields simultaneously
        #    - Correlation studies between fields
        #    - Statistical analysis and data export
        #    - Quick initial result verification
        #
        #    USE INDIVIDUAL FILES FOR:
        #    - Mesh deformation visualization (Wrap by Vector) → u/u.xdmf
        #    - Field-specific analysis and filtering
        #    - Publication-quality visualizations
        #    - Large simulations (smaller file size per field)
        #
        merged_path = os.path.join(self.output_folder, "solution", "solution.xdmf")
        os.makedirs(os.path.dirname(merged_path), exist_ok=True)
        self.merged_output = do.io.XDMFFile(
            self.eq.grid.mesh.comm,
            merged_path,
            "w",
        )
        self.merged_output.write_mesh(self.eq.grid.mesh)

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
            self.output_fields[i].write_function(field, t)
            # Write to merged solution file for unified analysis across all fields.
            # For specific field visualization (e.g., deformation via Wrap by Vector),
            # use individual field files in ParaView instead.
            self.merged_output.write_function(field, t)

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
