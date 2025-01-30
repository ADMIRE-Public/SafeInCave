import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
import numpy as np
import pandas as pd
from Utils import read_json
import matplotlib.pyplot as plt
from ResultsHandler import read_vector_from_points

hour = 60*60
MPa = 1e6

def plot_FEM(ax, output_folder):

	# Read displacement results
	pvd_path = os.path.join(output_folder, "operation", "vtk", "displacement")
	pvd_file = "displacement.pvd"
	df, u, v, w = read_vector_from_points(pvd_path, pvd_file)

	point_P = df[(df["z"]==1) & (df["x"]==1) & (df["y"]==1)].index[0]

	Lz = df["z"].values.max()
	Lx = df["x"].values.max()

	u_P = u.iloc[point_P].values[1:]
	w_P = w.iloc[point_P].values[1:]

	eps_h = -u_P / Lx
	eps_v = -w_P / Lz

	t = w.iloc[point_P].index.values[1:]/hour

	# Read input file
	input_file = read_json("input_file.json")

	# Extract vertical stresses
	time_list = np.array(input_file["time_settings"]["time_list"])/hour
	sigma_v_list = np.array(input_file["boundary_conditions"]["TOP"]["values"])/MPa
	sigma_v = np.interp(t, time_list, sigma_v_list)

	ax.plot(t, eps_v*100, ".", color="steelblue", label=r"$\varepsilon_1$ - FEM")
	ax.plot(t, eps_h*100, ".", color="lightcoral", label=r"$\varepsilon_3$ - FEM")
	ax.set_xlabel("Time (h)", size=12, fontname="serif")
	ax.set_ylabel("Total strain (%)", size=12, fontname="serif")
	ax.legend(loc=0, shadow=True, fancybox=True, prop={"size": 8})
	ax.grid(True, color="0.92")
	ax.set_facecolor("0.85")

def plot_MP(ax1, ax2, ax3, output_folder):
	results = read_json(os.path.join(output_folder, "material_point", "results.json"))

	# Extract results
	time = np.array(results["time"]) / hour
	epsilon_1 = np.array(results["total"]["epsilon_1"]) * 100
	epsilon_3 = np.array(results["total"]["epsilon_3"]) * 100
	sigma_1 = np.array(results["stress"]["sigma_1"]) / MPa
	sigma_3 = np.array(results["stress"]["sigma_3"]) / MPa

	# Plot stresses
	ax1.plot(time, sigma_1, ".-", color="steelblue", label=r"$\sigma_1$")
	ax1.plot(time, sigma_3, ".-", color="lightcoral", label=r"$\sigma_3$")
	ax1.set_xlabel("Time (h)", size=10, fontname="serif")
	ax1.set_ylabel("Stress (MPa)", size=10, fontname="serif")
	ax1.legend(loc=0, shadow=True, fancybox=True, prop={"size": 8})
	ax1.grid(True, color="0.92")
	ax1.set_facecolor("0.85")

	# Plot total strains
	ax2.plot(time, epsilon_1, "-", color="steelblue", label=r"$\varepsilon_1$ - MP")
	ax2.plot(time, epsilon_3, "-", color="lightcoral", label=r"$\varepsilon_3$ - MP")
	ax2.set_xlabel("Time (h)", size=10, fontname="serif")
	ax2.set_ylabel("Total strain (%)", size=10, fontname="serif")
	ax2.legend(loc=0, shadow=True, fancybox=True, prop={"size": 8})
	ax2.grid(True, color="0.92")
	ax2.set_facecolor("0.85")

	# Plot individual element contributions
	for element in results.keys():
		if element != "total" and element != "stress" and element != "time":
			eps_1 = np.array(results[element]["epsilon_1"]) * 100
			eps_3 = np.array(results[element]["epsilon_3"]) * 100
			ax3.plot(time, eps_1, "-", label=element)
	ax3.set_xlabel("Time (h)", size=10, fontname="serif")
	ax3.set_ylabel("Partial strains (%)", size=10, fontname="serif")
	ax3.legend(loc=0, shadow=True, fancybox=True, prop={"size": 8})
	ax3.grid(True, color="0.92")
	ax3.set_facecolor("0.85")

def main():
	results_folder = os.path.join("output", "case_0")

	# Plot loading schedule
	fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 3))
	fig.subplots_adjust(top=0.970, bottom=0.155, left=0.043, right=0.980, hspace=0.35, wspace=0.225)
	fig.patch.set_alpha(0.0)

	plot_MP(ax1, ax2, ax3, results_folder)
	plot_FEM(ax2, results_folder)

	plt.show()


if __name__ == '__main__':
	main()