import safeincave as sf
import safeincave.PostProcessingTools as post
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Button, Slider
import meshio

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
	# wall_ind = np.unique(mesh.cells_dict["line"].flatten())
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
	ax.legend(loc=0, ncol=1, prop={"size": 8})

def plot_dilatancy_boundary(ax):
	dilation_points = np.array([
				[-5.289256198347102, 0.06089309878213811],
				[-3.3057851239669382, 1.8876860622462814],
				[-0.11019283746556141, 3.897158322056839],
				[3.4159779614325068, 5.7848443843031205],
				[6.611570247933887, 7.4289580514208495],
				[10.24793388429752, 8.951285520974295],
				[13.553719008264462, 10.351826792963472],
				[17.190082644628095, 11.81326116373478],
				[20.495867768595044, 13.152909336941818],
				[23.80165289256198, 14.431664411366718],
				[27.21763085399449, 15.588633288227339],
				[30.523415977961424, 16.745602165087963],
				[34.0495867768595, 17.96346414073072],
				[37.465564738292, 19.120433017591342],
				[40.881542699724505, 20.15561569688769],
				[44.29752066115701, 21.190798376184034],
				[47.82369146005509, 22.34776725304466],
				[51.23966942148759, 23.322056833558864],
				[54.87603305785122, 24.41813261163735],
				[58.402203856749296, 25.39242219215156],
				[61.8181818181818, 26.305818673883632],
				[65.45454545454544, 27.341001353179976],
				[68.9807162534435, 28.254397834912048],
				[72.50688705234158, 29.289580514208392],
				[76.03305785123966, 30.202976995940464],
				[79.33884297520659, 30.994587280108256] ])
	I1 = dilation_points[:,0]
	J2_sqrt = dilation_points[:,1]
	p_points = I1/3
	q_points = J2_sqrt*np.sqrt(3)
	ax.plot(p_points, q_points, "-", color="black")


def plot_paths(ax, stress_data, point_color, index=0):
	cells_coord, sigma_v, q = stress_data
	# sigma_v = -sigma_v/MPa
	# sigma_v /= -MPa
	# q /= MPa
	# print(sigma_v)
	# print(q)

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

	ax.plot(sigma_v[cell_p][0], q[cell_p][0], "o", color="blue", label="Start")
	ax.plot(sigma_v[cell_p][-1], q[cell_p][-1], "^", color="red", label="Finish")
	ax.legend(loc=2, ncol=2, prop={"size": 8})
	# if label_start:
	# 	ax.plot(sigma_v[cell_p][0], q[cell_p][0], "o", color="blue", label="Start")
	# 	ax.plot(sigma_v[cell_p][-1], q[cell_p][-1], "^", color="red", label="Finish")
	# 	ax.legend(loc=2, ncol=2, prop={"size": 8})
	# 	# ax.legend(bbox_to_anchor=(0.92, 1.25), ncol=2)
	# else:
	# 	ax.plot(sigma_v[cell_p][0], q[cell_p][0], "o", color="blue")
	# 	ax.plot(sigma_v[cell_p][-1], q[cell_p][-1], "^", color="red")

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

def get_simulation_times(log_file):
	lines = log_file.split("\n")
	times = []
	for line in lines:
		if "seconds)" in line:
			i1 = line.find("(") + 1
			i2 = line.find(" seconds)")
			time = float(line[i1:i2])
			times.append(time)
	return sum(times)



def plot_results_panel(results_folder, stage="operation"):
	# Define paths
	output_path = os.path.join("output", results_folder, stage)

	# # Read log file
	# with open(os.path.join("output", results_folder, "log.txt"), "r") as file:
	#     log_file = file.read()

	# # Read CPU time
	# cpu_time = get_simulation_times(log_file)
	# cpu_gmtime = time.strftime("%H:%M:%S", time.gmtime(cpu_time))

	# Read mesh
	msh_file_name = os.path.join(output_path, "mesh", "geom.msh")
	points_msh, cells_msh = post.read_msh_as_pandas(msh_file_name)

	# Build mapping
	xdmf_file_name = os.path.join(output_path, "u", "u.xdmf")
	mapping = post.find_mapping(points_msh, xdmf_file_name)

	# Displacement data
	df_ux, df_uy, df_uz = post.read_vector_from_points(xdmf_file_name, mapping)
	displacement_data = points_msh, df_ux, df_uy, df_uz

	# Stress data
	points_xdmf, cells_xdmf = post.read_xdmf_as_pandas(xdmf_file_name)
	mid_cells = post.compute_cell_centroids(points_xdmf.values, cells_xdmf.values)
	df_p = post.read_scalar_from_cells(os.path.join(output_path, "p_elems", "p_elems.xdmf"))
	df_q = post.read_scalar_from_cells(os.path.join(output_path, "q_elems", "q_elems.xdmf"))
	stress_data = (mid_cells, -df_p.values/MPa, df_q.values/MPa)

	# Read simulation time steps
	time_steps = df_ux.columns.values
	n_time = len(time_steps)

	# # Read input file
	# input_file = read_json(os.path.join("output", results_folder, stage, "input_file.json"))
	# grid_path = input_file["grid"]["path"]

	# Get indices of wall profile
	# mesh = meshio.read(os.path.join(grid_path, "geom.msh"))
	mesh = meshio.read(msh_file_name)

	# # Read gas pressure
	# gas_time = np.array(input_file["time_settings"]["time_list"])
	# gas_pressure = np.array(input_file["boundary_conditions"]["Cavern"]["values"])

	# Calculate cavern convergence results
	xi, zi, xf, zf, times, volumes = calculate_convergence_data(displacement_data, mesh)

	# Plot pressure schedule
	fig = plt.figure(figsize=(16, 9))
	fig.subplots_adjust(top=0.975, bottom=0.120, left=0.060, right=0.986, hspace=0.44, wspace=0.64)

	gs = GridSpec(18, 19, figure=fig)
	ax_logo = fig.add_subplot(gs[0:2,0:4])
	ax_info_1 = fig.add_subplot(gs[0:2,5:9])
	ax_info_2 = fig.add_subplot(gs[0:2,10:14])
	ax_info_3 = fig.add_subplot(gs[0:2,15:])

	ax0 = fig.add_subplot(gs[3:12,0:4])

	ax00 = fig.add_subplot(gs[3:7,5:9])
	ax01 = fig.add_subplot(gs[3:7,10:14])
	ax02 = fig.add_subplot(gs[3:7,15:])

	ax10 = fig.add_subplot(gs[8:12,5:9])
	ax11 = fig.add_subplot(gs[8:12,10:14])
	ax12 = fig.add_subplot(gs[8:12,15:])

	ax30 = fig.add_subplot(gs[14:,0:5])
	ax31 = fig.add_subplot(gs[14:,6:12])
	ax32 = fig.add_subplot(gs[14:,13:19])

	apply_grey_theme(fig, [ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30, ax31, ax32], transparent=True)

	# plot_gas_pressure(ax30, time_steps, gas_time, gas_pressure)
	plot_subsidence(ax31, displacement_data)
	plot_cavern_shape(ax0, xi, zi, xf, zf)
	plot_convergence(ax32, times, volumes)

	img = plt.imread(os.path.join("..", "..", "..", "assets", "logo_2.png"))
	ax_logo.imshow(img)
	ax_logo.text(910, 295, "Version 2.0.0")
	ax_logo.axis('off')

	# Plot grid info
	# n_elems = mesh.cells[0].data.shape[0]
	n_elems = mesh.cells["tetra"].shape[0]
	n_nodes = len(mesh.points)

	# print()
	# print(mesh.__dict__.keys())
	# print()
	# print(mesh.field_data)
	# print()
	# print(mesh.field_data.keys())

	# print()

	region_names = ""
	for field_name in mesh.field_data.keys():
		dimension = mesh.field_data[field_name][1]
		if dimension == 3:
			region_names += field_name + ", "
	region_names = region_names[:-2]

	# region_names = ""
	# for region_name in input_file["grid"]["regions"].keys():
	# 	region_names += region_name + ", "
	# region_names = region_names[:-2]

	ax_info_1.text(0, 0.8, "Mesh info:", size=12, fontname="serif")
	# ax_info_1.text(0, 0.5, f"- Location: {grid_path}", size=10, fontname="serif")
	ax_info_1.text(0, 0.2, f"- Number of elems: {n_elems}", size=10, fontname="serif")
	ax_info_1.text(0, -0.1, f"- Number of nodes: {n_nodes}", size=10, fontname="serif")
	ax_info_1.text(0, -0.4, f"- Regions: {region_names}", size=10, fontname="serif")
	ax_info_1.axis('off')

	# # Plot constitutive model info
	# elements = []
	# for elem_type in ["elastic", "viscoelastic", "inelastic"]:
	# 	for elem_name in input_file["constitutive_model"][elem_type].keys():
	# 		if input_file["constitutive_model"][elem_type][elem_name]["active"] == True:
	# 			elements.append(elem_name)

	# dh = 0.3
	# h = 0.8
	# i = 1
	# ax_info_2.text(0, h, "Constitutive model: ", size=12, fontname="serif")
	# for element in elements:
	# 	ax_info_2.text(0, h-i*dh, f"- {element}", size=10, fontname="serif")
	# 	i += 1
	ax_info_2.axis('off')

	# ax_info_3.text(0, 0.8, "Simulation info:", size=12, fontname="serif")
	# ax_info_3.text(0, 0.5, f"- CPU time: {cpu_gmtime}", size=10, fontname="serif")
	ax_info_3.axis('off')

	points = get_relevant_points()

	plot_paths(ax00, stress_data, points[0])
	plot_paths(ax10, stress_data, points[5])
	plot_paths(ax02, stress_data, points[2])
	plot_paths(ax01, stress_data, points[1])
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
		# plot_gas_pressure(ax30, time_steps, gas_time, gas_pressure, index)
		# t_interp = time_steps[max(0,index-1)]
		# p_interp = np.interp(t_interp, gas_time, gas_pressure)
		# ax30.scatter(t_interp/day, p_interp/MPa, c="white", edgecolors="black")
		# ax30.set_xlim(xmin,xmax)
		# ax30.set_ylim(ymin,ymax)

		stress_path_list = [(ax00, 0), (ax10, 5), (ax02, 2), (ax01, 1), (ax12, 3), (ax11, 4)]
		for ax, i in stress_path_list:
			xmin,xmax = ax.get_xlim()
			ymin,ymax = ax.get_ylim()
			ax.cla()
			plot_paths(ax, stress_data, points[i], index)
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
	plot_results_panel("case_2", "operation")

if __name__ == '__main__':
	main()
