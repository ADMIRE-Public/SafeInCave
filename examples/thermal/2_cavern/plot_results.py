import safeincave as sf
import safeincave.PostProcessingTools as post
import matplotlib.pyplot as plt
import numpy as np
import os

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

def reorder_data(df_coord, T, line_idx):
	# Initial cavern shape
	x0 = df_coord["x"]
	# Reorder all coordinates according to coordinate x
	sorted_x0_ind = x0.sort_values().index
	x0 = x0[sorted_x0_ind]
	# Reorder all displacements according to coordinate x
	T = T.loc[sorted_x0_ind]
	return x0, T


def main():
	# Choose results folder
	results_folder = os.path.join("output", "case_0")

	# Read temperature field
	points, time_list, T_field = post.read_node_scalar(os.path.join(results_folder, "T", "T.xdmf"))

	print(points.shape)
	print(time_list.shape)
	print(T_field.shape)

	# Extract point at top of the cavern z = 430
	top_idx = post.find_closest_point(np.array([0, 0, 430]), points)
	bottom_idx = post.find_closest_point(np.array([0, 0, 205.19]), points)

	print(top_idx, points[top_idx])
	print(bottom_idx, points[bottom_idx])

	t = time_list/86400 # Convert time to days
	T_top = T_field[:, top_idx]
	T_bottom = T_field[:, bottom_idx]

	# Plot pressure schedule
	fig, ax = plt.subplots(1, 1, figsize=(5, 3.5))
	fig.subplots_adjust(top=0.90, bottom=0.15, left=0.143, right=0.980, hspace=0.35, wspace=0.293)

	# # Plot figures
	ax.plot(t, T_top, "-", color="#377eb8", label="T_top (z=430 m)")
	ax.plot(t, T_bottom, "-", color="#ff7f00", label="T_bottom (z=205.19 m)")
	ax.set_xlabel("Time (days)", fontname="serif", fontsize=12)
	ax.set_ylabel("Temperature (K)", fontname="serif", fontsize=12)
	ax.legend(loc=0, shadow=True, fancybox=True, prop={"size":8})

	apply_grey_theme(fig, [ax], transparent=True, grid_color="0.92", back_color='0.85')

	plt.show()


if __name__ == '__main__':
	main()