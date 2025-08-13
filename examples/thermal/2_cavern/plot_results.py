# import os
# import sys
# sys.path.append(os.path.join("..", "..", "..", "safeincave"))
# import numpy as np
# import pandas as pd
# from Utils import read_json, MPa
# import matplotlib.pyplot as plt
# import meshio
# from PostProcessingTools import (read_scalar_from_points,
# 								find_mapping,
# 								read_xdmf_as_pandas,
# 								read_msh_as_pandas)

import safeincave as sf
import safeincave.PostProcessingTools as post
import matplotlib.pyplot as plt
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

	# Plot pressure schedule
	fig, ax = plt.subplots(1, 1, figsize=(5, 3.5))
	fig.subplots_adjust(top=0.90, bottom=0.15, left=0.143, right=0.980, hspace=0.35, wspace=0.293)

	# Read mesh
	msh_file_name = os.path.join(results_folder, "mesh", "geom.msh")
	points_msh, cells_msh = post.read_msh_as_pandas(msh_file_name)

	# Read displacements
	xdmf_file_name = os.path.join(results_folder, "T", "T.xdmf")
	mapping = post.find_mapping(points_msh, xdmf_file_name)
	T = post.read_scalar_from_points(xdmf_file_name, mapping)

	# Extract point at top of the cavern z = 430
	top_idx = points_msh[(points_msh["z"] == 430) & (points_msh["x"] == 0) & (points_msh["y"] == 0)].index

	# Extract point at top of the cavern z = 205.19
	bottom_idx = points_msh[(abs(points_msh["z"] - 205.19)< 0.1) & (points_msh["x"] == 0) & (points_msh["y"] == 0)].index

	# Get temperature and coordinates along the line (y,z)=(0,1)
	T_top = T.loc[top_idx]
	T_bottom = T.loc[bottom_idx]


	# # Reorder according to increasing x coodinate values
	# x, T_line = reorder_data(x, T_line, line_idx)

	# Extract time in days
	t = T_top.iloc[0].index.values/60/60/24

	# Plot figures
	ax.semilogx(t, T_top.values[0], "-", color="#377eb8", label="T_top (z=430 m)")
	ax.semilogx(t, T_bottom.values[0], "-", color="#ff7f00", label="T_bottom (z=205.19 m)")
	ax.set_xlabel("x (m)", fontname="serif", fontsize=12)
	ax.set_ylabel("Temperature (K)", fontname="serif", fontsize=12)
	ax.legend(loc=0, shadow=True, fancybox=True, prop={"size":8})

	# ax2.plot(t, T_line.values[0,:])
	# ax2.set_xlabel("Time (days)", fontname="serif", fontsize=12)
	# ax2.set_ylabel("Temperature (K)", fontname="serif", fontsize=12)
	
	apply_grey_theme(fig, [ax], transparent=True, grid_color="0.92", back_color='0.85')

	plt.show()


if __name__ == '__main__':
	main()