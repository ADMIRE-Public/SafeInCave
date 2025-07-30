import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
import numpy as np
import pandas as pd
from Utils import read_json, MPa
import matplotlib.pyplot as plt
import meshio
from PostProcessingTools import (read_vector_from_points,
								find_mapping,
								read_scalar_from_cells,
								read_xdmf_as_pandas,
								compute_cell_centroids,
								read_msh_as_pandas)

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

def plot_displacements(ax, results_folder):
	# Read mesh
	msh_file_name = os.path.join(results_folder, "mesh", "geom.msh")
	points_msh, cells_msh = read_msh_as_pandas(msh_file_name)

	# Read displacements
	xdmf_file_name = os.path.join(results_folder, "u", "u.xdmf")
	mapping = find_mapping(points_msh, cells_msh, xdmf_file_name)
	u, v, w = read_vector_from_points(xdmf_file_name, mapping)

	point_A = points_msh[(points_msh["z"] == 1) & (points_msh["x"] == 0) & (points_msh["y"] == 0)].index[0]
	point_B = points_msh[(points_msh["z"] == 1) & (points_msh["x"] == 0) & (points_msh["y"] == 1)].index[0]
	point_C = points_msh[(points_msh["z"] == 1) & (points_msh["x"] == 1) & (points_msh["y"] == 1)].index[0]
	point_D = points_msh[(points_msh["z"] == 1) & (points_msh["x"] == 1) & (points_msh["y"] == 0)].index[0]
	print(point_A, point_B, point_C, point_D)
	print("Point A: ", points_msh.iloc[point_A].values)
	print("Point B: ", points_msh.iloc[point_B].values)
	print("Point C: ", points_msh.iloc[point_C].values)
	print("Point D: ", points_msh.iloc[point_D].values)

	w_A = w.iloc[point_A].values[1:]
	w_B = w.iloc[point_B].values[1:]
	w_C = w.iloc[point_C].values[1:]
	w_D = w.iloc[point_D].values[1:]

	t = w.iloc[point_A].index.values[1:]/60

	ax.plot(t, w_A*1000, ".-", color="#377eb8", label="Point A")
	ax.plot(t, w_B*1000, ".-", color="#ff7f00", label="Point B")
	ax.plot(t, w_C*1000, ".-", color="#4daf4a", label="Point C")
	ax.plot(t, w_D*1000, ".-", color="#f781bf", label="Point D")
	ax.set_xlabel("Time (minutes)", size=12, fontname="serif")
	ax.set_ylabel("Displacement (mm)", size=12, fontname="serif")
	ax.grid(True)
	ax.legend(loc=0, shadow=True, fancybox=True)


def plot_q(ax, results_folder):
	points_xdmf, cells_xdmf = read_xdmf_as_pandas(os.path.join(results_folder, "q_elems", "q_elems.xdmf"))
	mid_cells = compute_cell_centroids(points_xdmf.values, cells_xdmf.values)

	ind_A = mid_cells.index[mid_cells["y"] < 0.5]
	ind_B = mid_cells.index[mid_cells["y"] > 0.5]

	q_df = read_scalar_from_cells(os.path.join(results_folder, "q_elems", "q_elems.xdmf"))

	q_A = np.average(q_df.iloc[ind_A].values, axis=0)[1:]/MPa
	q_B = np.average(q_df.iloc[ind_B].values, axis=0)[1:]/MPa

	t = q_df.iloc[0].index.values[1:]/60

	ax.plot(t, q_A, ".-", color="steelblue", label=r"$\Omega_A$")
	ax.plot(t, q_B, ".-", color="lightcoral", label=r"$\Omega_B$")
	ax.set_xlabel("Time (minutes)", size=12, fontname="serif")
	ax.set_ylabel("Average Von Mises stress (MPa)", size=12, fontname="serif")
	ax.legend(loc=0, fancybox=True, shadow=True)



def main():
	results_folder = os.path.join("output", "case_0")

	

	# Plot pressure schedule
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))
	fig.subplots_adjust(top=0.970, bottom=0.135, left=0.093, right=0.980, hspace=0.35, wspace=0.225)

	plot_displacements(ax1, results_folder)
	plot_q(ax2, results_folder)
	apply_grey_theme(fig, [ax1, ax2], transparent=True, grid_color="0.92", back_color='0.85')

	plt.show()


if __name__ == '__main__':
	main()