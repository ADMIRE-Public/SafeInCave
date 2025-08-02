import os
import sys
sys.path.append(os.path.join("..", "..", "..", "safeincave"))
import numpy as np
import pandas as pd
import meshio
import time
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Button, Slider
from PostProcessingTools import read_xdmf_as_pandas, read_msh_as_pandas, read_vector_from_points, read_scalar_from_points, find_mapping, compute_cell_centroids
import json

def read_json(file_name):
	with open(file_name, "r") as j_file:
		data = json.load(j_file)
	return data

hour = 60*60
day = 24*hour
MPa = 1e6

def apply_grey_theme(fig, axes, transparent=True, grid_color="0.92", back_color='0.85'):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.grid(True, color=grid_color)
			ax.set_axisbelow(True)
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black')
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			ax.tick_params(axis='x', colors='black', which='both')
			ax.tick_params(axis='y', colors='black', which='both')
			ax.yaxis.label.set_color('black')
			ax.xaxis.label.set_color('black')
			ax.set_facecolor(back_color)


def rotate_z(coord, angle):
	rotation_matrix_z = np.array([
	    [np.cos(angle), -np.sin(angle), 0],
	    [np.sin(angle), np.cos(angle), 0],
	    [0, 0, 1]
	])
	return np.dot(coord, rotation_matrix_z.T)

def get_selected_points(point_base, point_coords):
	selected_idx = []
	thetas = np.linspace(0, np.pi/2, 80)
	for theta in thetas:
		rotating_point = rotate_z(point_base, theta)
		x_p, y_p, z_p = rotating_point
		d = np.sqrt(  (point_coords[:,0] - x_p)**2
		            + (point_coords[:,1] - y_p)**2
		            + (point_coords[:,2] - z_p)**2 )
		point_idx = d.argmin()
		if abs(point_coords[point_idx,2] - z_p) < 0.01:
			if point_idx not in selected_idx:
				selected_idx.append(point_idx)
	return selected_idx

def get_p_q(point_base, point_coords, df_p, df_q):
	# Get selected points
	selected_idxs = get_selected_points(point_base, point_coords)

	# Calculate average mean stress along selected points
	p_selected = df_p.iloc[selected_idxs].values
	p_avg = np.average(p_selected, axis=0)

	# Calculate average von Mises stress along selected points
	q_selected = df_q.iloc[selected_idxs].values
	q_avg = np.average(q_selected, axis=0)

	return p_avg, q_avg


def main():
	output_path = os.path.join("output", "case_2", "equilibrium")

	# Read mesh
	msh_file_name = os.path.join(output_path, "mesh", "geom.msh")
	points_msh, cells_msh = read_msh_as_pandas(msh_file_name)

	# Stress data
	xdmf_file_name = os.path.join(output_path, "p_nodes", "p_nodes.xdmf")
	points_xdmf, cells_xdmf = read_xdmf_as_pandas(xdmf_file_name)
	df_p = read_scalar_from_points(os.path.join(output_path, "p_nodes", "p_nodes.xdmf"))
	df_q = read_scalar_from_points(os.path.join(output_path, "q_nodes", "q_nodes.xdmf"))
	point_coords = points_xdmf.values

	# Define point
	depth = 800 + 430
	point_base = np.array([74.63, 0.0, 267.4 - depth])

	# Get selected points
	selected_idxs = get_selected_points(point_base, point_coords)

	# Get average values
	p_avg, q_avg = get_p_q(point_base, point_coords, df_p, df_q)

	print(q_avg)




if __name__ == '__main__':
	main()
