import safeincave as sf
import safeincave.PostProcessingTools as post
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Button, Slider
import meshio

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


def create_layout():
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
	return fig, ax_logo, ax_info_1, ax_info_2, ax_info_3, ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30, ax31, ax32


class WallProfileData():
	def __init__(self, output_folder, scale=1):
		# Read displacements from xdmf
		points, self.time_list, u_field = post.read_node_vector(os.path.join(output_folder, "u", "u.xdmf"))

		# Read msh node coordinates from msh file
		reader_msh = meshio.read(os.path.join(output_folder, "mesh", "geom.msh"))
		points_msh = reader_msh.points

		# Get wall indices from msh file and remap it to xdmf data
		wall_idx = np.unique(reader_msh.cells["line"].flatten())
		mapping = post.build_mapping(points_msh, points)
		for i in range(len(wall_idx)):
			wall_idx[i] = mapping[wall_idx[i]]

		# Get coordinates from wall nodes
		self.wall_points = points[wall_idx]

		# Extract displacements from wall nodes
		self.wall_u = u_field[:,wall_idx]

		# Reorder wall nodes and displacements according to z coordinate
		sorted_idx = np.argsort(self.wall_points[:,2])
		self.wall_points = self.wall_points[sorted_idx]
		self.wall_u = self.wall_u[:,sorted_idx]

		# Calculate initial cavern shape
		self.scale = scale

		self.compute_volumes()

	def get_wall_coordinates(self, time_step: int):
		return self.wall_points + self.scale*self.wall_u[time_step,:]

	def compute_volumes(self):
		wall_t0 = self.get_wall_coordinates(time_step=0)
		vol_0 = self.__trapezoidal_volume(wall_t0[:,2], wall_t0[:,0])
		self.volumes = []
		for t_step in range(len(self.time_list)):
			wall_t = self.get_wall_coordinates(time_step=t_step)
			vol = self.__trapezoidal_volume(wall_t[:,2], wall_t[:,0])
			self.volumes.append(100*abs(vol_0 - vol)/vol_0)

	def __trapezoidal_volume(self, x, y):
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


def plot_cavern_shape(ax, wall: WallProfileData, t_step: int):
	wall_t0 = wall.get_wall_coordinates(time_step=0)
	wall_tf = wall.get_wall_coordinates(time_step=t_step)

	ax.plot(wall_t0[:,0], wall_t0[:,2], "-", color="black", label="Initial shape")
	ax.plot(wall_tf[:,0], wall_tf[:,2], "-", color="steelblue", label="Final shape")
	ax.set_xlabel("x (m)", size=12, fontname="serif")
	ax.set_ylabel("z (m)", size=12, fontname="serif")
	ax.legend(loc=1, shadow=True, fancybox=True, prop={"size": 8})
	ax.axis("equal")


def plot_convergence(ax, wall: WallProfileData, t_step: int):
	ax.plot(wall.time_list/day, wall.volumes, "-", color="#377eb8", linewidth="2.0")
	ax.scatter(wall.time_list[t_step]/day, wall.volumes[t_step], c="white", edgecolors="black", zorder=10000)
	ax.set_xlabel("Time (days)", size=12, fontname="serif")
	ax.set_ylabel("Convergence (%)", size=12, fontname="serif")


def plot_logo(ax):
	img = plt.imread(os.path.join("..", "..", "..", "assets", "logo_2.png"))
	ax.imshow(img)
	ax.text(910, 295, "Version 2.0.0")
	ax.axis('off')


def show_metadata(ax_mesh, ax_model, ax_solver, output_folder):
	# reader_msh = meshio.read(os.path.join(output_folder, "mesh", "geom.msh"))
	with open(os.path.join(output_folder, "log.txt")) as f:
		log = f.readlines()

	dh = 0.3
	h = 0.8

	for i, line in enumerate(log):
		if "| Mesh info:" in line:
			line_split = log[i+4].split("|")
			n_elems = int(line_split[1])
			n_nodes = int(line_split[2])
			location = line_split[3].strip(" ")
			ax_mesh.text(0, h-0*dh, "Mesh info:", size=12, fontname="serif")
			ax_mesh.text(0, h-1*dh, f"- Location: {location}", size=10, fontname="serif")
			ax_mesh.text(0, h-2*dh, f"- Number of elems: {n_elems}", size=10, fontname="serif")
			ax_mesh.text(0, h-3*dh, f"- Number of nodes: {n_nodes}", size=10, fontname="serif")
			ax_mesh.axis('off')
		elif "| Solver info:" in line:
			line_split = log[i+4].split("|")
			ksp_solver = line_split[1].strip(" ")
			ksp_pc = line_split[2].strip(" ")
			tol = float(line_split[3])
			max_ite = int(line_split[4])
			cpu_time = log[-1].strip(" ").strip("Total time: ")[:8]
			ax_solver.text(0, h-0*dh, "Simulation Info: ", size=12, fontname="serif")
			ax_solver.text(0, h-1*dh, f"- Solver: {ksp_solver}, PC: {ksp_pc}", size=10, fontname="serif")
			# ax_solver.text(0, h-2*dh, f"- PC: {ksp_pc}", size=10, fontname="serif")
			ax_solver.text(0, h-2*dh, f"- tol: {tol}, max_ite: {max_ite}", size=10, fontname="serif")
			ax_solver.text(0, h-3*dh, f"- CPU Time: {cpu_time}", size=10, fontname="serif")
			ax_solver.axis('off')
		elif "| Constitutive model:" in line:
			# line_split = log[i+4].split("|")
			elastic_elems = log[i+4].split("|")[2].strip(" ")
			nonelastic_elems = log[i+5].split("|")[2].strip(" ")
			thermoelastic_elems = log[i+6].split("|")[2].strip(" ")
			ax_model.text(0, h-0*dh, "Constitutive model: ", size=12, fontname="serif")
			ax_model.text(0, h-1*dh, f"- Elastic: {elastic_elems}", size=10, fontname="serif")
			ax_model.text(0, h-2*dh, f"- Non-elastic: {nonelastic_elems}", size=10, fontname="serif")
			ax_model.text(0, h-3*dh, f"- Thermoelastic: {thermoelastic_elems}", size=10, fontname="serif")
			ax_model.axis('off')


class SubsidenceData():
	def __init__(self, output_folder):
		points, self.time_list, u_field = post.read_node_vector(os.path.join(output_folder, "u", "u.xdmf"))
		top_idx = post.find_closest_point([0,0,660], points)
		self.subsidence = u_field[:,top_idx,2] - u_field[0,top_idx,2]


def plot_subsidence(ax, sub: SubsidenceData, t_step: int):
	ax.plot(sub.time_list/day, 100*sub.subsidence, "-")
	ax.scatter(sub.time_list[t_step]/day, 100*sub.subsidence[t_step], c="white", edgecolors="black", zorder=10000)
	ax.set_xlabel("Time (days)", size=12, fontname="serif")
	ax.set_ylabel("Subsidence (cm)", size=12, fontname="serif")


class StressPath():
	def __init__(self, output_folder, probes):
		self.probes = probes
		self.points, self.time_list, self.p_elems = post.read_cell_scalar(os.path.join(output_folder, "p_elems", "p_elems.xdmf"))
		self.points, self.time_list, self.q_elems = post.read_cell_scalar(os.path.join(output_folder, "q_elems", "q_elems.xdmf"))

		self.p_probes = np.zeros((probes.shape[0], self.time_list.size))
		self.q_probes = np.zeros((probes.shape[0], self.time_list.size))
		for i, probe in enumerate(probes):
			idx = post.find_closest_point(probe, self.points)
			self.p_probes[i,:] = -self.p_elems[:,idx]/MPa
			self.q_probes[i,:] = self.q_elems[:,idx]/MPa

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


COLORS = ["deepskyblue", "tomato", "orange", "steelblue", "purple", "magenta"]

def plot_stress_paths(ax, stress_path: StressPath, i: int, t_step: int):
	ax.plot(stress_path.p_probes[i], stress_path.q_probes[i], "-", color=COLORS[i])
	ax.scatter(stress_path.p_probes[i,t_step], stress_path.q_probes[i,t_step], c="white", edgecolors="black", zorder=10000)
	plot_dilatancy_boundary(ax)
	ax.set_xlabel("Mean stress (MPa)", size=10, fontname="serif")
	ax.set_ylabel("Von Mises stress (MPa)", size=10, fontname="serif")

def plot_probes(ax, probes):
	for i, probe in enumerate(probes):
		ax.scatter(probe[0], probe[2], c=COLORS[i], edgecolors="black", zorder=10000)


def main():
	fig, ax_logo, ax_info_1, ax_info_2, ax_info_3, ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30, ax31, ax32 = create_layout()

	# Specify folder to load the results from
	output_folder = os.path.join("output", "case_0", "operation")
	# output_folder = os.path.join("output", "case_gui", "operation")

	# Probe points around the cavern wall
	probes = np.array([
		[0, 0, 430],
		[42.8, 0, 393.4],
		[45, 0, 345],
		[57.62, 0, 301.3],
		[74.63, 0, 267.4],
		[0, 0, 205.1],
	])

	wall_profile = WallProfileData(output_folder, scale=1)
	subsidence_data = SubsidenceData(output_folder)
	stress_path = StressPath(output_folder, probes)

	n_steps = stress_path.time_list.size

	plot_logo(ax_logo)
	show_metadata(ax_info_1, ax_info_2, ax_info_3, output_folder)

	plot_cavern_shape(ax0, wall_profile, t_step=-1)
	plot_convergence(ax31, wall_profile, t_step=0)
	plot_subsidence(ax32, subsidence_data, t_step=0)
	plot_stress_paths(ax00, stress_path, i=0, t_step=0)
	plot_stress_paths(ax01, stress_path, i=1, t_step=0)
	plot_stress_paths(ax02, stress_path, i=2, t_step=0)
	plot_stress_paths(ax10, stress_path, i=3, t_step=0)
	plot_stress_paths(ax11, stress_path, i=4, t_step=0)
	plot_stress_paths(ax12, stress_path, i=5, t_step=0)
	plot_probes(ax0, probes)

	def update_plot(val):
		t_step = int(val)

		xmin,xmax = ax31.get_xlim()
		ymin,ymax = ax31.get_ylim()
		ax31.cla()
		plot_convergence(ax31, wall_profile, t_step=t_step)
		ax31.set_xlim(xmin,xmax)
		ax31.set_ylim(ymin,ymax)

		xmin,xmax = ax32.get_xlim()
		ymin,ymax = ax32.get_ylim()
		ax32.cla()
		plot_subsidence(ax32, subsidence_data, t_step=t_step)
		ax32.set_xlim(xmin,xmax)
		ax32.set_ylim(ymin,ymax)

		for i, ax in enumerate([ax00, ax01, ax02, ax10, ax11, ax12]):
			xmin,xmax = ax.get_xlim()
			ymin,ymax = ax.get_ylim()
			ax.cla()
			plot_stress_paths(ax, stress_path, i=i, t_step=t_step)
			ax.set_xlim(xmin,xmax)
			ax.set_ylim(ymin,ymax)




		apply_grey_theme(fig, [ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30, ax31, ax32], transparent=True)

	# Make a horizontal slider to control the frequency.
	axtime = fig.add_axes([0.09, 0.02, 0.75, 0.03])
	time_slider = Slider(
	    ax=axtime,
	    label='Time step',
	    valmin=0,
	    valmax=n_steps-1,
	    valinit=0,
	)

	# register the update function with each slider
	time_slider.on_changed(update_plot)

	plt.show()


if __name__ == '__main__':
	main()