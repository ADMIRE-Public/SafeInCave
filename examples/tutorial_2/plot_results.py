import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
import numpy as np
import pandas as pd
import meshio
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Button, Slider
from ResultsHandler import read_vector_from_points, read_scalar_from_cells
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

def calculate_convergence_data(displacement_data, mesh):
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

	df_coord, u, v, w = displacement_data

	# Get indices of wall profile
	wall_ind = np.unique(mesh.cells["line"].flatten())

	# Get reordered data over cavern wall
	x0, y0, z0, u, v, w = reorder_data(df_coord, u, v, w, wall_ind)

	# Get times
	times = u.columns.values
	t_final = times[-1]
	t_initial = times[1]

	# Compute cavern volumes over time
	vol_0 = trapezoidal_volume(z0.values, x0.values)
	volumes = []
	for t in times[1:]:
		z = z0.values + w[t].values
		x = x0.values + u[t].values
		vol = trapezoidal_volume(z, x)
		volumes.append(100*abs(vol_0 - vol)/vol_0)

	# Plot cavern shape
	expansion_factor = 50
	xi = x0 + 0*expansion_factor*u[t_initial]
	yi = y0 + 0*expansion_factor*v[t_initial]
	zi = z0 + 0*expansion_factor*w[t_initial]
	xf = x0 + expansion_factor*u[t_final]
	yf = y0 + expansion_factor*v[t_final]
	zf = z0 + expansion_factor*w[t_final]

	return xi, zi, xf, zf, times, volumes

def plot_subsidence(ax, displacement_data, index=0):
	df_coord, u, v, w = displacement_data
	point_A = df_coord[(df_coord["z"] == 660) & (df_coord["x"] == 0) & (df_coord["y"] == 0)].index[0]
	w_A = w.iloc[point_A].values[0:] - w.iloc[point_A].values[0]
	t = w.iloc[point_A].index.values[0:]
	max_val = len(t)
	i = min(max_val, max(0, index))
	ax.plot(t/day, w_A*100, "-", color="#377eb8", linewidth="2.0")
	ax.scatter(t[i]/day, w_A[i]*100, c="white", edgecolors="black", zorder=10000)
	ax.set_xlabel("Time (days)", size=12, fontname="serif")
	ax.set_ylabel("Subsidence (cm)", size=12, fontname="serif")
	ax.grid(True)

def plot_convergence(ax, times, volumes, index=0):
	max_val = len(times[2:])
	i = min(max_val, max(1, index))
	ax.plot(times[2:]/day, volumes[1:], "-", color="#377eb8", linewidth="2.0")
	ax.scatter(times[i+1]/day, volumes[i], c="white", edgecolors="black", zorder=10000)
	ax.set_xlabel("Time (days)", size=12, fontname="serif")
	ax.set_ylabel("Convergence (%)", size=12, fontname="serif")

def plot_cavern_shape(ax, xi, zi, xf, zf):
	ax.plot(xi, zi, "-", color="black", linewidth=2.0, label="Initial shape")
	ax.plot(xf, zf, "-", color="#377eb8", linewidth=2.0, label=f"Final shape")
	ax.set_xlabel("x (m)", size=12, fontname="serif")
	ax.set_ylabel("z (m)", size=12, fontname="serif")
	ax.legend(loc=1, shadow=True, fancybox=True)
	ax.axis("equal")
	ax.legend(bbox_to_anchor=(1.09, 1.11), ncol=2)

def plot_dilatancy_boundary(ax):
	dilation_points = np.array([
								[-5.039370078740159, 0.08637236084453548],
								[-0.3149606299212593, 2.5047984644913655],
								[4.251968503937007, 4.232245681381968],
								[10.078740157480315, 6.132437619961621],
								[16.850393700787404, 7.946257197696745],
								[23.307086614173233, 9.500959692898277],
								[29.133858267716544, 11.142034548944345],
								[33.385826771653555, 12.1785028790787],
								[39.21259842519686, 13.474088291746646],
								[45.669291338582696, 14.856046065259122],
								[50.708661417322844, 16.06525911708254],
								[56.2204724409449, 17.188099808061427],
								[63.6220472440945, 18.570057581573902],
								[70.07874015748034, 19.606525911708257],
								[75.27559055118112, 20.815738963531675],
								[79.84251968503939, 21.679462571976973] ])
	dilation_points[:,1] *= 1.7
	# dilation_points[:,0] *= np.sqrt(3)
	ax.plot(dilation_points[:,0], dilation_points[:,1], "-", color="black")


def plot_paths(ax, stress_data, point_color, index=0, label_start=False):
	cells_coord, sigma_v, q = stress_data
	sigma_v = -sigma_v

	# Find closest cell to point
	x_p, y_p, z_p = point_color[0]
	d = np.sqrt(  (cells_coord["x"].values - x_p)**2
	            + (cells_coord["y"].values - y_p)**2
	            + (cells_coord["z"].values - z_p)**2 )
	cell_p = d.argmin()
	x_cell_p = cells_coord["x"].values[cell_p]
	y_cell_p = cells_coord["y"].values[cell_p]
	z_cell_p = cells_coord["z"].values[cell_p]

	plot_dilatancy_boundary(ax)
	ax.plot(sigma_v[cell_p], q[cell_p], "-", color=point_color[1])
	ax.scatter(sigma_v[cell_p][index], q[cell_p][index], c="white", edgecolors="black", zorder=10000)
	ax.set_xlabel("Mean stress (MPa)", size=10, fontname="serif")
	ax.set_ylabel("Von Mises stress (MPa)", size=10, fontname="serif")
	if label_start:
		ax.plot(sigma_v[cell_p][0], q[cell_p][0], "o", color="blue", label="Start")
		ax.plot(sigma_v[cell_p][-1], q[cell_p][-1], "^", color="red", label="Finish")
		ax.legend(bbox_to_anchor=(0.92, 1.25), ncol=2)
	else:
		ax.plot(sigma_v[cell_p][0], q[cell_p][0], "o", color="blue")
		ax.plot(sigma_v[cell_p][-1], q[cell_p][-1], "^", color="red")

def plot_probe_points(ax, points):
	for point, color in points:
		ax.scatter(point[0], point[2], marker="o", alpha=1.0, edgecolors="k", color=color, zorder=10000)

def plot_gas_pressure(ax, time_steps, time, pressure, index=0):
	max_val = len(time[2:])
	i = min(max_val, max(1, index))
	t_interp = time_steps[max(0,index-1)]
	p_interp = np.interp(t_interp, time, pressure)
	ax.scatter(t_interp/day, p_interp/MPa, c="white", edgecolors="black")
	# ax.scatter(time[i+1]/day, pressure[i]/MPa, c="white", edgecolors="black")
	ax.plot(time/day, pressure/MPa, "-", color="black")
	ax.set_xlabel("Time (days)", size=10, fontname="serif")
	ax.set_ylabel("Gas pressure (MPa)", size=10, fontname="serif")

def get_relevant_points():
	x1, z1 = 0, 430
	x2, z2 = 0, 205.1
	x3, z3 = 74.63, 267.4
	x4, z4 = 57.62, 301.3
	x5, z5 = 45, 345
	x6, z6 = 42.8, 393.4

	point_1 = (x1, 0.0, z1)
	point_2 = (x2, 0.0, z2)
	point_3 = (x3, 0.0, z3)
	point_4 = (x4, 0.0, z4)
	point_5 = (x5, 0.0, z5)
	point_6 = (x6, 0.0, z6)

	points = [
		(point_1, "deepskyblue"),
		(point_2, "tomato"),
		(point_3, "orange"),
		(point_4, "steelblue"),
		(point_5, "purple"),
		(point_6, "magenta"),
	]
	return points



def plot_results_panel(results_folder):
	# Define paths
	output_path = os.path.join("output", results_folder, "vtk")

	# Read displacement results
	df_coord_nodes, df_ux, df_uy, df_uz = read_vector_from_points(os.path.join(output_path, "displacement"), "displacement.pvd")
	displacement_data = df_coord_nodes, df_ux, df_uy, df_uz

	# Read simulation time steps
	time_steps = df_ux.columns.values
	n_time = len(time_steps)

	# Read input file
	input_file = read_json(os.path.join("output", results_folder, "input_file.json"))
	grid_path = input_file["grid"]["path"]

	# Get indices of wall profile
	mesh = meshio.read(os.path.join(grid_path, "geom.msh"))
	wall_ind = np.unique(mesh.cells["line"].flatten())

	# # Read stress results
	df_coord_cells, df_sigma_v = read_scalar_from_cells(os.path.join(output_path, "p_smooth"), "p_smooth.pvd")
	df_coord_cells, df_q = read_scalar_from_cells(os.path.join(output_path, "q_smooth"), "q_smooth.pvd")
	stress_data = (df_coord_cells, df_sigma_v.values, df_q.values)

	# Read gas pressure
	gas_time = np.array(input_file["time_settings"]["time_list"])
	gas_pressure = np.array(input_file["boundary_conditions"]["Cavern"]["values"])

	# Calculate cavern convergence results
	xi, zi, xf, zf, times, volumes = calculate_convergence_data(displacement_data, mesh)

	# Plot pressure schedule
	fig = plt.figure(figsize=(14, 8))
	fig.subplots_adjust(top=0.930, bottom=0.120, left=0.050, right=0.967, hspace=0.44, wspace=0.64)

	gs = GridSpec(14, 19, figure=fig)
	ax0 = fig.add_subplot(gs[0:9,0:4])

	ax00 = fig.add_subplot(gs[0:4,5:9])
	ax01 = fig.add_subplot(gs[0:4,10:14])
	ax02 = fig.add_subplot(gs[0:4,15:])

	ax10 = fig.add_subplot(gs[5:9,5:9])
	ax11 = fig.add_subplot(gs[5:9,10:14])
	ax12 = fig.add_subplot(gs[5:9,15:])

	ax30 = fig.add_subplot(gs[10:,0:5])
	ax31 = fig.add_subplot(gs[10:,6:12])
	ax32 = fig.add_subplot(gs[10:,13:19])

	apply_grey_theme(fig, [ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30, ax31, ax32], transparent=True)

	plot_gas_pressure(ax30, time_steps, gas_time, gas_pressure)
	plot_subsidence(ax31, displacement_data)
	plot_cavern_shape(ax0, xi, zi, xf, zf)
	plot_convergence(ax32, times, volumes)

	points = get_relevant_points()

	plot_paths(ax00, stress_data, points[0])
	plot_paths(ax10, stress_data, points[5])
	plot_paths(ax02, stress_data, points[2])
	plot_paths(ax01, stress_data, points[1], label_start=True)
	plot_paths(ax12, stress_data, points[3])
	plot_paths(ax11, stress_data, points[4])

	plot_probe_points(ax0, points)


	# The function to be called anytime a slider's value changes
	def update_plot(val):
		global index
		index = max(0, int(val))

		xmin,xmax = ax30.get_xlim()
		ymin,ymax = ax30.get_ylim()
		ax30.cla()
		plot_gas_pressure(ax30, time_steps, gas_time, gas_pressure, index)
		t_interp = time_steps[max(0,index-1)]
		p_interp = np.interp(t_interp, gas_time, gas_pressure)
		ax30.scatter(t_interp/day, p_interp/MPa, c="white", edgecolors="black")
		ax30.set_xlim(xmin,xmax)
		ax30.set_ylim(ymin,ymax)

		stress_path_list = [(ax00, 0, False), (ax10, 5, False), (ax02, 2, False), (ax01, 1, True), (ax12, 3, False), (ax11, 4, False)]
		for ax, i, label in stress_path_list:
			xmin,xmax = ax.get_xlim()
			ymin,ymax = ax.get_ylim()
			ax.cla()
			plot_paths(ax, stress_data, points[i], index, label_start=label)
			ax.set_xlim(xmin,xmax)
			ax.set_ylim(ymin,ymax)

		xmin,xmax = ax32.get_xlim()
		ymin,ymax = ax32.get_ylim()
		ax32.cla()
		plot_convergence(ax32, times, volumes, index)
		ax32.set_xlim(xmin,xmax)
		ax32.set_ylim(ymin,ymax)

		xmin,xmax = ax31.get_xlim()
		ymin,ymax = ax31.get_ylim()
		ax31.cla()
		plot_subsidence(ax31, displacement_data, index)
		ax31.set_xlim(xmin,xmax)
		ax31.set_ylim(ymin,ymax)

		apply_grey_theme(fig, [ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30, ax31, ax32], transparent=True)

		fig.canvas.draw_idle()

	# Make a horizontal slider to control the frequency.
	axtime = fig.add_axes([0.09, 0.02, 0.75, 0.03])
	time_slider = Slider(
	    ax=axtime,
	    label='Time step',
	    valmin=0,
	    valmax=n_time-1,
	    valinit=0,
	)

	# register the update function with each slider
	time_slider.on_changed(update_plot)

	plt.show()

def main():
	results_folder = os.path.join("case_0", "equilibrium")
	plot_results_panel(results_folder)

if __name__ == '__main__':
	main()
