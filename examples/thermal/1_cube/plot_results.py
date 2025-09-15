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

# def find_points(points, fun):



def main():
	results_folder = os.path.join("output", "case_0")

	points, t, T = post.read_node_scalar(os.path.join(results_folder, "T", "T.xdmf"))
	t /= (60*60*24)

	# Extract indices along line (y,z)=(0,1)
	line_idx = np.where((points[:,1] == 0.0) & (points[:,2] == 1.0))[0]

	# Extract points along line
	line_points = points[line_idx]

	# Extract x coordinates along line
	x0 = line_points[:,0]

	# Extract T along line
	T0 = T[:,line_idx]

	# Sort indices according to increasing x
	sorted_idx = np.argsort(x0)

	# Sort x and T
	x_sorted = x0[sorted_idx]
	T_sorted = T0[:,sorted_idx]

	# Plot pressure schedule
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))
	fig.subplots_adjust(top=0.90, bottom=0.15, left=0.11, right=0.980, hspace=0.35, wspace=0.293)

	# Plot figures
	ax1.plot(x_sorted, T_sorted[0,:], ".-", color="#377eb8", label=f"t={round(t[0],2)} day(s)")
	ax1.plot(x_sorted, T_sorted[3,:], ".-", color="#ff7f00", label=f"t={round(t[3],2)} day(s)")
	ax1.plot(x_sorted, T_sorted[15,:], ".-", color="#4daf4a", label=f"t={round(t[15],2)} day(s)")
	ax1.plot(x_sorted, T_sorted[30,:], ".-", color="#f781bf", label=f"t={round(t[30],2)} day(s)")
	ax1.set_xlabel("x (m)", fontname="serif", fontsize=12)
	ax1.set_ylabel("Temperature (K)", fontname="serif", fontsize=12)
	ax1.legend(loc=0, shadow=True, fancybox=True, prop={"size":8})

	ax2.plot(t, T_sorted[:,5])
	ax2.set_xlabel("Time (days)", fontname="serif", fontsize=12)
	ax2.set_ylabel("Temperature (K)", fontname="serif", fontsize=12)

	
	apply_grey_theme(fig, [ax1, ax2], transparent=True, grid_color="0.92", back_color='0.85')

	plt.show()


if __name__ == '__main__':
	main()