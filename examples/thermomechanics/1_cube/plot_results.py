import os
import sys
sys.path.append(os.path.join("..", "..", "..", "safeincave"))
import numpy as np
import pandas as pd
from Utils import read_json, MPa
import matplotlib.pyplot as plt
import meshio
from PostProcessingTools import (read_vector_from_points,
								find_mapping,
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


def main():
	# Choose results folder
	results_folder = os.path.join("output", "case_0")


	# Plot pressure schedule
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.5))
	fig.subplots_adjust(top=0.92, bottom=0.155, left=0.093, right=0.980, hspace=0.35, wspace=0.295)

	# Read mesh
	msh_file_name = os.path.join(results_folder, "mesh", "geom.msh")
	points_msh, cells_msh = read_msh_as_pandas(msh_file_name)

	# Read displacements
	xdmf_file_name = os.path.join(results_folder, "u", "u.xdmf")
	mapping = find_mapping(points_msh, cells_msh, xdmf_file_name)
	u, v, w = read_vector_from_points(xdmf_file_name, mapping)

	point_1 = points_msh[(points_msh["z"] == 1) & (points_msh["x"] == 1) & (points_msh["y"] == 1)].index[0]
	point_2 = points_msh[(points_msh["z"] == 1) & (points_msh["x"] == 1) & (points_msh["y"] == 0)].index[0]
	print(point_1, point_2)
	print("Point 1: ", points_msh.iloc[point_1].values)
	print("Point 2: ", points_msh.iloc[point_2].values)

	ux_1 = u.iloc[point_1].values*1000
	uy_1 = v.iloc[point_1].values*1000
	uz_1 = w.iloc[point_1].values*1000

	ux_2 = u.iloc[point_2].values*1000
	uy_2 = v.iloc[point_2].values*1000
	uz_2 = w.iloc[point_2].values*1000

	t = w.iloc[point_2].index.values/60/60/24

	ax1.semilogx(t, ux_1, "-", color="#377eb8", label=r"$u_x$")
	ax1.semilogx(t, uy_1, "-", color="#ff7f00", label=r"$u_y$")
	ax1.semilogx(t, uz_1, "-", color="#4daf4a", label=r"$u_z$")
	ax1.set_xlabel("Time (days)", size=12, fontname="serif")
	ax1.set_ylabel("Displacement (mm)", size=12, fontname="serif")
	ax1.set_title("Point 1 (1,1,1)", size=12, fontname="serif")
	ax1.legend(loc=0, shadow=True, fancybox=True)

	ax2.semilogx(t, ux_2, "-", color="#377eb8", label=r"$u_x$")
	ax2.semilogx(t, uy_2, "-", color="#ff7f00", label=r"$u_y$")
	ax2.semilogx(t, uz_2, "-", color="#4daf4a", label=r"$u_z$")
	ax2.set_xlabel("Time (days)", size=12, fontname="serif")
	ax2.set_ylabel("Displacement (mm)", size=12, fontname="serif")
	ax2.set_title("Point 2 (1,0,1)", size=12, fontname="serif")
	ax2.legend(loc=0, shadow=True, fancybox=True)

	apply_grey_theme(fig, [ax1, ax2], transparent=True, grid_color="0.92", back_color='0.85')

	plt.show()


if __name__ == '__main__':
	main()