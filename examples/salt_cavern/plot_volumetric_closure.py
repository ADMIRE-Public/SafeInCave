import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
from ResultsHandler import convert_vtk_to_pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
import json
import meshio
import pandas as pd

minute = 60
hour = 60*minute
day = 24*hour
MPa = 1e6

def read_json(file_name):
	with open(file_name, "r") as j_file:
		data = json.load(j_file)
	return data

def apply_grey_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.grid(True, color='0.92')
			ax.set_axisbelow(True)
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black')
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			ax.tick_params(axis='x', colors='black', which='both')
			ax.tick_params(axis='y', colors='black', which='both')
			ax.yaxis.label.set_color('black')
			ax.xaxis.label.set_color('black')
			ax.set_facecolor("0.85")

def apply_dark_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#dede00")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.grid(True, color='0.87')
			ax.set_axisbelow(True)
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black')
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			ax.tick_params(axis='x', colors='black', which='both')
			ax.tick_params(axis='y', colors='black', which='both')
			ax.yaxis.label.set_color('black')
			ax.xaxis.label.set_color('black')
			ax.set_facecolor("0.78")

def trapezoidal_volume(x, y):
	""" This function calculates the volume of a solid of revolution (around y=0 axis) based on the trapezoidal rule. """
	volume = 0.0
	area = 0.0
	n = len(x)
	for i in range(1, n):
		R = 0.5*(y[i] + y[i-1])
		A = np.pi*R**2
		d = x[i] - x[i-1]
		area += R*d
		volume += A*d
	return volume

def compute_closure(x0, y0, z0, u, v, w, times):
	vol_0 = trapezoidal_volume(z0.values, x0.values)
	# vol_0 = trapezoidal_volume(z0.values + w[times[0]].values, x0.values + u[times[0]].values)
	volumes = []
	for t in times[1:]:
		z = z0.values + w[t].values
		x = x0.values + u[t].values
		vol = trapezoidal_volume(z, x)
		volumes.append(100*abs(vol_0 - vol)/vol_0)
	return volumes

def get_wall_indices(mesh_file):
	mesh = meshio.read(mesh_file)
	wall_ind = np.unique(mesh.cells["line"].flatten())
	# wall_coord = mesh.points[wall_ind]
	return wall_ind

def reorder_data(df_coord, u, v, w, wall_ind):
	# Initial cavern shape
	x0 = df_coord.iloc[wall_ind]["x"]
	y0 = df_coord.iloc[wall_ind]["y"]
	z0 = df_coord.iloc[wall_ind]["z"]
	# Reorder all coordinates according to coordinate z
	sorted_z0_ind = z0.sort_values().index
	x0 = x0[sorted_z0_ind]
	y0 = y0[sorted_z0_ind]
	z0 = z0[sorted_z0_ind]
	# Reorder all displacements according to coordinate z
	u = u.iloc[wall_ind].loc[sorted_z0_ind]
	v = v.iloc[wall_ind].loc[sorted_z0_ind]
	w = w.iloc[wall_ind].loc[sorted_z0_ind]
	return x0, y0, z0, u, v, w

def plot_cavern_shape(ax, x0, y0, z0, u, v, w, t_final, case_color="blue", show_0=False):
	# Final cavern shape
	factor = 50
	xf = x0 + factor*u[t_final]
	yf = y0 + factor*v[t_final]
	zf = z0 + factor*w[t_final]
	if show_0:
		ax.plot(x0, z0, "-", color="black", linewidth=2.0, label="Initial shape")
		ax.plot(-x0, z0, "-", color="black", linewidth=2.0)
	ax.plot(xf, zf, "-", color=case_color, linewidth=2.0, label=f"Final shape")
	ax.plot(-xf, zf, "-", color=case_color, linewidth=2.0)
	# ax.axis("equal")
	ax.set_xlabel("x (m)", size=12, fontname="serif")
	ax.set_ylabel("z (m)", size=12, fontname="serif")
	ax.grid(True)
	ax.legend(loc=1, shadow=True, fancybox=True)

def plot_results(ax_shape, ax_closure, results_folder, mesh_folder, case_name="name", line_style="-", case_color="purple"):
	# Read displacement results
	pvd_path = os.path.join(results_folder, "displacement")
	pvd_file = "displacement.pvd"
	df_coord, u, v, w = convert_vtk_to_pandas(pvd_path, pvd_file)

	# Get indices of wall profile
	wall_ind = get_wall_indices(os.path.join(mesh_folder, "geom.msh"))

	# Get reordered data over cavern wall
	x0, y0, z0, u, v, w = reorder_data(df_coord, u, v, w, wall_ind)

	# Get times
	times = u.columns.values
	t_final = times[-1]

	# Compute cavern volumes over time
	volumes = compute_closure(x0, y0, z0, u, v, w, times)

	# Plot cavern shape
	plot_cavern_shape(ax_shape, x0, y0, z0, u, v, w, t_final, case_color=case_color, show_0=True)
	ax_shape.set_title(case_name, size=14, fontname="serif")

	# Plot cavern closure
	ax_closure.plot(times[2:]/hour, volumes[1:] - 0*volumes[1], line_style, color=case_color, label=case_name, linewidth="2.0")
	ax_closure.set_xlabel("Time (h)", size=12, fontname="serif")
	ax_closure.set_ylabel("Cavern closure (%)", size=12, fontname="serif")
	ax_closure.legend(loc=0, shadow=True, fancybox=True)
	ax_closure.grid(True)

def main():
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3))
	fig.subplots_adjust(top=0.985, bottom=0.145, left=0.060, right=0.990, hspace=0.35, wspace=0.260)

	results_folder = os.path.join("output", "case_0", "operation", "vtk")
	mesh_folder = os.path.join("..", "..", "grids", "cavern_regular")
	plot_results(ax1, ax2, results_folder, mesh_folder, case_name="With viscoplasticity", line_style=".-", case_color="#377eb8")

	ax1.axis("equal")

	apply_grey_theme(fig, [ax1, ax2], transparent=True)
	# apply_dark_theme(fig, [ax_inset], transparent=True)

	plt.show()

def main_2():
	fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9, 3))
	fig.subplots_adjust(top=0.985, bottom=0.145, left=0.070, right=0.990, hspace=0.35, wspace=0.310)

	results_folder = os.path.join("output", "case_0", "operation", "vtk")
	mesh_folder = os.path.join("..", "..", "grids", "cavern_regular")
	plot_results(ax1, ax3, results_folder, mesh_folder, case_name="With viscoplasticity", line_style=".-", case_color="#377eb8")

	results_folder = os.path.join("output", "case_1", "operation", "vtk")
	mesh_folder = os.path.join("..", "..", "grids", "cavern_regular")
	plot_results(ax2, ax3	, results_folder, mesh_folder, case_name="No viscoplasticity", line_style=".-", case_color="#ff7f00")

	# ax1.set_xlim(-63, 186)
	# ax1.set_ylim(193, 482)
	ax1.axis("equal")
	ax2.axis("equal")

	apply_grey_theme(fig, [ax1, ax2, ax3], transparent=True)
	# apply_dark_theme(fig, [ax_inset], transparent=True)

	plt.show()

if __name__ == '__main__':
	# main()
	main_2()