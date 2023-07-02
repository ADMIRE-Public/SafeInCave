import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import pandas as pd
import numpy as np
import os

hour = 60*60

def apply_grey_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.grid(True, color='0.92')
			ax.set_axisbelow(True)
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black')
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			ax.tick_params(axis='x', colors='black', which='both')
			ax.tick_params(axis='y', colors='black', which='both')
			ax.yaxis.label.set_color('black')
			ax.xaxis.label.set_color('black')
			ax.set_facecolor("0.85")

def plot_model(axis, folder):
	skip = 1
	eps_files = ["eps_tot", "eps_e", "eps_ve", "eps_cr_dis", "eps_cr_sol", "eps_damage"]
	eps_names = [r"$\varepsilon_{tot}$", r"$\varepsilon_{e}$", r"$\varepsilon_{ve}$", r"$\varepsilon_{cr}$", r"$\varepsilon_{ps}$", r"$\varepsilon_{d}$"]
	color_names = ["steelblue", "orange", "lightcoral", "violet", "dodgerblue", "mediumseagreen"]

	for eps_file, eps_name, color_name in zip(eps_files, eps_names, color_names):
		try:
			eps_tot = pd.read_excel(os.path.join(folder, f"{eps_file}.xlsx"))
			time = eps_tot["Time"].values[::skip]/hour
			eps_a = -100*eps_tot["22"].values[::skip]
			eps_r = -100*eps_tot["00"].values[::skip]
			axis[0].plot(time, eps_a, "-", color=color_name, label=eps_name, linewidth=2.0)
			axis[1].plot(time, eps_r, "-", color=color_name, linewidth=2.0)
		except:
			pass

	axis[0].set_xlabel("Time (hours)", size=12, fontname="serif")
	axis[0].set_ylabel("Axial strain (%)", size=12, fontname="serif")
	axis[0].grid(True)
	axis[0].legend(bbox_to_anchor=(2.65, 1.0), shadow=True, fancybox=True, ncol=1)

	axis[1].set_xlabel("Time (hours)", size=12, fontname="serif")
	axis[1].set_ylabel("Radial strain (%)", size=12, fontname="serif")
	axis[1].grid(True)

def main():
	fig, axis = plt.subplots(1, 2, figsize=(8, 3))
	fig.subplots_adjust(top=0.975, bottom=0.155, left=0.085, right=0.875, hspace=0.2, wspace=0.32)

	results_folder = os.path.join("output", "case_0", "fem", "avg")
	plot_model(axis, results_folder)
	apply_grey_theme(fig, axis.flatten(), transparent=True)
	plt.show()

if __name__ == "__main__":
	main()