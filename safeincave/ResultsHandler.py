"""
Useful to read vtk files and post-process results.
"""
# Copyright 2024 The safeincave community.
#
# This file is part of safeincave.
#
# Licensed under the GNU GENERAL PUBLIC LICENSE, Version 3 (the "License"); you may not
# use this file except in compliance with the License.  You may obtain a copy
# of the License at
#
#     https://spdx.org/licenses/GPL-3.0-or-later.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations under
# the License.

import xml.etree.ElementTree as ET
import meshio
import pandas as pd
import numpy as np
import os

def convert_vtk_to_pandas(pvd_path, pvd_file):
    """
    This function reads vtk files containing the time dependent *displacement* solution
    and convert it to pandas dataframes

    Parameters
    ----------
    pvd_path : str
        Path to the .pvd file, which must be at the same directory as the .vtk files.
    pvd_file : str
        Name of the .pvd file (usually *displacement.pvd*).

    Returns
    -------
    df_coord : pandas.core.frame.DataFrame
        Spatial coordinates (*x*, *y*, *z*) of the grid nodes.
    df_ux : pandas.core.frame.DataFrame
        Component *x* of the displacement vector solution.
    df_uy : pandas.core.frame.DataFrame
        Component *y* of the displacement vector solution.
    df_uz : pandas.core.frame.DataFrame
        Component *z* of the displacement vector solution.
    """
    tree = ET.parse(os.path.join(pvd_path, pvd_file))
    root = tree.getroot()
    vtu_files = [os.path.join(pvd_path, child.get("file")) for child in root.findall(".//DataSet")]
    time_steps = [float(child.get("timestep")) for child in root.findall(".//DataSet")]

    mesh = meshio.read(os.path.join(vtu_files[0]))
    x = mesh.points[:,0]
    y = mesh.points[:,1]
    z = mesh.points[:,2]
    df_coord = pd.DataFrame({'x': x, 'y': y, 'z': z})
    n_points = len(x)

    ux_dict, uy_dict, uz_dict = {}, {}, {}
    for time_step, vtu_file in zip(time_steps, vtu_files):
        mesh = meshio.read(vtu_file)
        u = mesh.point_data["Displacement"]
        ux_dict[time_step] = u[:,0]
        uy_dict[time_step] = u[:,1]
        uz_dict[time_step] = u[:,2]
    df_ux = pd.DataFrame(ux_dict)
    df_uy = pd.DataFrame(uy_dict)
    df_uz = pd.DataFrame(uz_dict)

    return df_coord, df_ux, df_uy, df_uz

def read_vector_from_points(pvd_path, pvd_file):
    """
    This function reads vtk files containing the time dependent *displacement* solution
    at grid nodes and convert it to pandas dataframes.

    Parameters
    ----------
    pvd_path : str
        Path to the .pvd file, which must be at the same directory as the .vtk files.
    pvd_file : str
        Name of the .pvd file (usually *displacement.pvd*).

    Returns
    -------
    df_coord : pandas.core.frame.DataFrame
        Spatial coordinates (*x*, *y*, *z*) of the grid nodes.
    df_ux : pandas.core.frame.DataFrame
        Component *x* of the displacement vector solution.
    df_uy : pandas.core.frame.DataFrame
        Component *y* of the displacement vector solution.
    df_uz : pandas.core.frame.DataFrame
        Component *z* of the displacement vector solution.
    """
    tree = ET.parse(os.path.join(pvd_path, pvd_file))
    root = tree.getroot()
    vtu_files = [os.path.join(pvd_path, child.get("file")) for child in root.findall(".//DataSet")]
    time_steps = [float(child.get("timestep")) for child in root.findall(".//DataSet")]

    mesh = meshio.read(os.path.join(vtu_files[0]))
    x = mesh.points[:,0]
    y = mesh.points[:,1]
    z = mesh.points[:,2]
    df_coord = pd.DataFrame({'x': x, 'y': y, 'z': z})
    n_points = len(x)
    
    # Get field name
    for key in mesh.point_data.keys():
        field_name = key

    ux_dict, uy_dict, uz_dict = {}, {}, {}
    for time_step, vtu_file in zip(time_steps, vtu_files):
        mesh = meshio.read(vtu_file)
        u = mesh.point_data[field_name]
        ux_dict[time_step] = u[:,0]
        uy_dict[time_step] = u[:,1]
        uz_dict[time_step] = u[:,2]
    df_ux = pd.DataFrame(ux_dict)
    df_uy = pd.DataFrame(uy_dict)
    df_uz = pd.DataFrame(uz_dict)

    return df_coord, df_ux, df_uy, df_uz

def read_scalar_from_cells(pvd_path, pvd_file):
    """
    This function reads vtk files containing the time dependent solution of a scalar
    function defined on elements and convert it to pandas dataframes.

    Parameters
    ----------
    pvd_path : str
        Path to the .pvd file, which must be at the same directory as the .vtk files.
    pvd_file : str
        Name of the .pvd file (usually *displacement.pvd*).

    Returns
    -------
    df_coord : pandas.core.frame.DataFrame
        Spatial coordinates (*x*, *y*, *z*) of the centroids of all grid elements.
    df_scalar : pandas.core.frame.DataFrame
        Scalar values at each grid element.
    """
    tree = ET.parse(os.path.join(pvd_path, pvd_file))
    root = tree.getroot()
    vtu_files = [os.path.join(pvd_path, child.get("file")) for child in root.findall(".//DataSet")]
    time_steps = [float(child.get("timestep")) for child in root.findall(".//DataSet")]

    mesh = meshio.read(os.path.join(vtu_files[0]))

    # Get field name
    for key in mesh.cell_data["tetra"].keys():
        field_name = key

    # Vertex coordinates
    x = mesh.points[:,0]
    y = mesh.points[:,1]
    z = mesh.points[:,2]

    # Cell conectivities
    cells = mesh.cells["tetra"]
    n_elems, _ = cells.shape

    # Compute centroid coordinates of cells
    x_cells = np.zeros(n_elems)
    y_cells = np.zeros(n_elems)
    z_cells = np.zeros(n_elems)
    for i, cell in enumerate(cells):
        x_cells[i] = sum(x[cell])/4
        y_cells[i] = sum(y[cell])/4
        z_cells[i] = sum(z[cell])/4

    df_coord = pd.DataFrame({'x': x_cells, 'y': y_cells, 'z': z_cells})

    field_dict = {}
    for time_step, vtu_file in zip(time_steps, vtu_files):
        mesh = meshio.read(vtu_file)
        field_dict[time_step] = mesh.cell_data["tetra"][field_name]
    df_scalar = pd.DataFrame(field_dict)

    return df_coord, df_scalar

def read_tensor_from_cells_old(pvd_path, pvd_file):
    """
    This function reads vtk files containing the time dependent solution of a 
    rank-2 tensor function **A** defined on elements and convert it to pandas dataframes.

    Parameters
    ----------
    pvd_path : str
        Path to the .pvd file, which must be at the same directory as the .vtk files.
    pvd_file : str
        Name of the .pvd file (usually *displacement.pvd*).

    Returns
    -------
    df_coord : pandas.core.frame.DataFrame
        Spatial coordinates (*x*, *y*, *z*) of the centroids of all grid elements.
    df_sx : pandas.core.frame.DataFrame
        Values of component :math:`A_{xx}`.
    df_sy : pandas.core.frame.DataFrame
        Values of component :math:`A_{yy}`.
    df_sz : pandas.core.frame.DataFrame
        Values of component :math:`A_{zz}`.
    df_sxy : pandas.core.frame.DataFrame
        Values of component :math:`A_{xy}`.
    df_sxz : pandas.core.frame.DataFrame
        Values of component :math:`A_{xz}`.
    df_syz : pandas.core.frame.DataFrame
        Values of component :math:`A_{yz}`.
    """
    tree = ET.parse(os.path.join(pvd_path, pvd_file))
    root = tree.getroot()
    vtu_files = [os.path.join(pvd_path, child.get("file")) for child in root.findall(".//DataSet")]
    time_steps = [float(child.get("timestep")) for child in root.findall(".//DataSet")]

    mesh = meshio.read(os.path.join(vtu_files[0]))

    # Get field name
    for key in mesh.cell_data["tetra"].keys():
        field_name = key

    # Vertex coordinates
    x = mesh.points[:,0]
    y = mesh.points[:,1]
    z = mesh.points[:,2]

    # Cell conectivities
    cells = mesh.cells["tetra"]
    n_elems, _ = cells.shape

    # Compute centroid coordinates of cells
    x_cells = np.zeros(n_elems)
    y_cells = np.zeros(n_elems)
    z_cells = np.zeros(n_elems)
    for i, cell in enumerate(cells):
        x_cells[i] = sum(x[cell])/4
        y_cells[i] = sum(y[cell])/4
        z_cells[i] = sum(z[cell])/4

    df_coord = pd.DataFrame({'x': x_cells, 'y': y_cells, 'z': z_cells})

    s_x, s_y, s_z, s_xy, s_xz, s_yz = {}, {}, {}, {}, {}, {}
    for time_step, vtu_file in zip(time_steps, vtu_files):
        mesh = meshio.read(vtu_file)
        stress = mesh.cell_data["tetra"][field_name]
        s_x[time_step] = stress[:,0]
        s_y[time_step] = stress[:,4]
        s_z[time_step] = stress[:,8]
        s_xy[time_step] = stress[:,1]
        s_xz[time_step] = stress[:,2]
        s_yz[time_step] = stress[:,5]
    df_sx = pd.DataFrame(s_x)
    df_sy = pd.DataFrame(s_y)
    df_sz = pd.DataFrame(s_z)
    df_sxy = pd.DataFrame(s_xy)
    df_sxz = pd.DataFrame(s_xz)
    df_syz = pd.DataFrame(s_yz)

    return df_coord, df_sx, df_sy, df_sz, df_sxy, df_sxz, df_syz

def read_tensor_from_cells(pvd_path, pvd_file):
    """
    This function reads vtk files containing the time dependent solution of a 
    rank-2 tensor function **A** defined on elements and convert it to pandas dataframes.

    Parameters
    ----------
    pvd_path : str
        Path to the .pvd file, which must be at the same directory as the .vtk files.
    pvd_file : str
        Name of the .pvd file (usually *displacement.pvd*).

    Returns
    -------
    df_coord : pandas.core.frame.DataFrame
        Spatial coordinates (*x*, *y*, *z*) of the centroids of all grid elements.
    df_sx : pandas.core.frame.DataFrame
        Values of component :math:`A_{xx}`.
    df_sy : pandas.core.frame.DataFrame
        Values of component :math:`A_{yy}`.
    df_sz : pandas.core.frame.DataFrame
        Values of component :math:`A_{zz}`.
    df_sxy : pandas.core.frame.DataFrame
        Values of component :math:`A_{xy}`.
    df_sxz : pandas.core.frame.DataFrame
        Values of component :math:`A_{xz}`.
    df_syz : pandas.core.frame.DataFrame
        Values of component :math:`A_{yz}`.
    """
    tree = ET.parse(os.path.join(pvd_path, pvd_file))
    root = tree.getroot()
    vtu_files = [os.path.join(pvd_path, child.get("file")) for child in root.findall(".//DataSet")]
    time_steps = [float(child.get("timestep")) for child in root.findall(".//DataSet")]

    mesh = meshio.read(os.path.join(vtu_files[0]))

    # Get field name
    for key in mesh.cell_data.keys():
        field_name = key

    # Vertex coordinates
    x = mesh.points[:,0]
    y = mesh.points[:,1]
    z = mesh.points[:,2]

    # Element connectivities
    connectivity = mesh.cells[0].data

    # Number of elements
    n_elems= len(connectivity)

    # Compute centroid coordinates of cells
    x_cells = np.zeros(n_elems)
    y_cells = np.zeros(n_elems)
    z_cells = np.zeros(n_elems)
    for i, cell in enumerate(connectivity):
        x_cells[i] = sum(x[cell])/4
        y_cells[i] = sum(y[cell])/4
        z_cells[i] = sum(z[cell])/4
        i += 1

    df_coord = pd.DataFrame({'x': x_cells, 'y': y_cells, 'z': z_cells})

    s_x, s_y, s_z, s_xy, s_xz, s_yz = {}, {}, {}, {}, {}, {}
    for time_step, vtu_file in zip(time_steps, vtu_files):
        mesh = meshio.read(vtu_file)
        stress = mesh.cell_data[field_name][0]
        s_x[time_step] = stress[:,0]
        s_y[time_step] = stress[:,4]
        s_z[time_step] = stress[:,8]
        s_xy[time_step] = stress[:,1]
        s_xz[time_step] = stress[:,2]
        s_yz[time_step] = stress[:,5]
    df_sx = pd.DataFrame(s_x)
    df_sy = pd.DataFrame(s_y)
    df_sz = pd.DataFrame(s_z)
    df_sxy = pd.DataFrame(s_xy)
    df_sxz = pd.DataFrame(s_xz)
    df_syz = pd.DataFrame(s_yz)

    return df_coord, df_sx, df_sy, df_sz, df_sxy, df_sxz, df_syz


if __name__ == '__main__':
    # Example of usage
    pvd_path = os.path.join("..", "apps", "sugar_cube", "output", "fem", "theta_0.5", "vtk")
    pvd_file = "displacement.pvd"

    df_coord, df_ux, df_uy, df_uz = convert_vtk_to_pandas(pvd_path, pvd_file)
    print(df_uz)
