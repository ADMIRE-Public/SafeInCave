import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ResultsHandler import read_vector_from_points, read_mesh_as_pandas

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


# Define results file
file_name = os.path.join("output", "case_0", "operation", "vtk", "displacement", "u.xdmf")

# Read mesh
df_points, df_cells = read_mesh_as_pandas(file_name)

# Read displacements
u, v, w = read_vector_from_points(file_name)

point_A = df_points[(df_points["z"] == 1) & (df_points["x"] == 0) & (df_points["y"] == 0)].index[0]
point_B = df_points[(df_points["z"] == 1) & (df_points["x"] == 0) & (df_points["y"] == 1)].index[0]
point_C = df_points[(df_points["z"] == 1) & (df_points["x"] == 1) & (df_points["y"] == 1)].index[0]
point_D = df_points[(df_points["z"] == 1) & (df_points["x"] == 1) & (df_points["y"] == 0)].index[0]
print(point_A, point_B, point_C, point_D)
print("Point A: ", df_points.iloc[point_A].values)
print("Point B: ", df_points.iloc[point_B].values)
print("Point C: ", df_points.iloc[point_C].values)
print("Point D: ", df_points.iloc[point_D].values)

w_A = w.iloc[point_A].values[1:]
w_B = w.iloc[point_B].values[1:]
w_C = w.iloc[point_C].values[1:]
w_D = w.iloc[point_D].values[1:]

t = w.iloc[point_A].index.values[1:]

# Plot pressure schedule
fig, ax = plt.subplots(1, 1, figsize=(5, 3.5))
fig.subplots_adjust(top=0.970, bottom=0.135, left=0.140, right=0.980, hspace=0.35, wspace=0.225)

ax.plot(t/60, w_A*1000, ".-", color="#377eb8", label="Point A")
ax.plot(t/60, w_B*1000, ".-", color="#ff7f00", label="Point B")
ax.plot(t/60, w_C*1000, ".-", color="#4daf4a", label="Point C")
ax.plot(t/60, w_D*1000, ".-", color="#f781bf", label="Point D")
ax.set_xlabel("Time (minutes)", size=12, fontname="serif")
ax.set_ylabel("Displacement (mm)", size=12, fontname="serif")
ax.grid(True)
ax.legend(loc=0, shadow=True, fancybox=True)

apply_grey_theme(fig, [ax], transparent=True, grid_color="0.92", back_color='0.85')

plt.show()
