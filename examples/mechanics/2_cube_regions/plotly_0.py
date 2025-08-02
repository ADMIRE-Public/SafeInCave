import os
import sys
sys.path.append(os.path.join("..", "..", "..", "safeincave"))
import numpy as np
import pandas as pd
from Utils import read_json, MPa
import meshio
from PostProcessingTools import (read_vector_from_points,
								find_mapping,
								read_scalar_from_cells,
								read_xdmf_as_pandas,
								compute_cell_centroids,
								read_msh_as_pandas)
import plotly.graph_objects as go
from plotly.subplots import make_subplots



def plot_displacements(fig, results_folder):
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

	w_A = w.iloc[point_A].values[1:]*1000
	w_B = w.iloc[point_B].values[1:]*1000
	w_C = w.iloc[point_C].values[1:]*1000
	w_D = w.iloc[point_D].values[1:]*1000

	t = w.iloc[point_A].index.values[1:]/60

	line = go.Scatter(x=t, y=w_A, name="<b>Point A</b>", line=dict(color='#377eb8', width=2), showlegend=True)
	fig.add_trace(line, row=1, col=1)

	line = go.Scatter(x=t, y=w_B, name="<b>Point B</b>", line=dict(color='#ff7f00', width=2), showlegend=True)
	fig.add_trace(line, row=1, col=1)

	line = go.Scatter(x=t, y=w_C, name="<b>Point C</b>", line=dict(color='#4daf4a', width=2), showlegend=True)
	fig.add_trace(line, row=1, col=1)

	line = go.Scatter(x=t, y=w_D, name="<b>Point D</b>", line=dict(color='#f781bf', width=2), showlegend=True)
	fig.add_trace(line, row=1, col=1)

	fig.update_xaxes(title_text="<b>Time (minutes)</b>", row=1, col=1)
	fig.update_yaxes(title_text="<b>Displacement (mm)</b>", row=1, col=1)

	# ax.plot(t, w_A*1000, ".-", color="#377eb8", label="Point A")
	# ax.plot(t, w_B*1000, ".-", color="#ff7f00", label="Point B")
	# ax.plot(t, w_C*1000, ".-", color="#4daf4a", label="Point C")
	# ax.plot(t, w_D*1000, ".-", color="#f781bf", label="Point D")
	# ax.set_xlabel("Time (minutes)", size=12, fontname="serif")
	# ax.set_ylabel("Displacement (mm)", size=12, fontname="serif")
	# ax.grid(True)
	# ax.legend(loc=0, shadow=True, fancybox=True)


def plot_q(fig, results_folder):
	points_xdmf, cells_xdmf = read_xdmf_as_pandas(os.path.join(results_folder, "q_elems", "q_elems.xdmf"))
	mid_cells = compute_cell_centroids(points_xdmf.values, cells_xdmf.values)

	ind_A = mid_cells.index[mid_cells["y"] < 0.5]
	ind_B = mid_cells.index[mid_cells["y"] > 0.5]

	q_df = read_scalar_from_cells(os.path.join(results_folder, "q_elems", "q_elems.xdmf"))

	q_A = np.average(q_df.iloc[ind_A].values, axis=0)[1:]/MPa
	q_B = np.average(q_df.iloc[ind_B].values, axis=0)[1:]/MPa

	t = q_df.iloc[0].index.values[1:]/60

	line = go.Scatter(x=t, y=q_A, name="Omega A", line=dict(color='#377eb8', width=2), showlegend=True)
	fig.add_trace(line, row=1, col=2)

	line = go.Scatter(x=t, y=q_B, name="Omega B", line=dict(color='#ff7f00', width=2), showlegend=True)
	fig.add_trace(line, row=1, col=2)

	fig.update_xaxes(title_text="<b>Time (minutes)</b>", row=1, col=2)
	fig.update_yaxes(title_text="<b>Average Von Mises stress (MPa)</b>", row=1, col=2)

	# ax.plot(t, q_A, ".-", color="steelblue", label=r"$\Omega_A$")
	# ax.plot(t, q_B, ".-", color="lightcoral", label=r"$\Omega_B$")
	# ax.set_xlabel("Time (minutes)", size=12, fontname="serif")
	# ax.set_ylabel("Average Von Mises stress (MPa)", size=12, fontname="serif")
	# ax.legend(loc=0, fancybox=True, shadow=True)



def main():
	results_folder = os.path.join("output", "case_0")

	
	# Create subplot with 1 row and 2 columns
	fig = make_subplots(rows=1, cols=2, subplot_titles=("First Plot - Two Curves", "Second Plot - One Bar Chart"))



	plot_displacements(fig, results_folder)
	plot_q(fig, results_folder)

	fig.write_html("plotly.html")


if __name__ == '__main__':
	main()