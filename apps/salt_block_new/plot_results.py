import os
import sys
sys.path.append(os.path.join("..", "..", "libs"))
from ResultsHandler import *
import matplotlib.pyplot as plt
import numpy as np

minute = 60
hour = 60*minute
day = 24*hour
MPa = 1e6

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

def plot_results(ax, pvd_path, color_name="steelblue", label_name=""):
	pvd_file = "displacement.pvd"

	df_coord, df_ux, df_uy, df_uz = convert_vtk_to_pandas(pvd_path, pvd_file)

	# Find cube dimensions
	Lx = df_coord["x"].max()
	Ly = df_coord["y"].max()
	Lz = df_coord["z"].max()

	# Find index of boundaries EAST (x = Lx), NORTH (y = Ly) and TOP (z = Lz)
	ind_east = df_coord["x"] == Lx
	ind_west = df_coord["y"] == Ly
	ind_top = df_coord["z"] == Lz

	eps_rad = -df_ux[ind_east].mean().values / Lx
	eps_axi = -df_uz[ind_top].mean().values / Lz
	time = df_ux.columns.values

	ax.plot(time/hour, 100*eps_axi, "-", color=color_name, label=label_name)
	ax.plot(time/hour, 100*eps_rad, "-", color=color_name)

def main():
	fig1, ax = plt.subplots(1, 1, figsize=(10, 4))
	fig1.subplots_adjust(top=0.945, bottom=0.115, left=0.085, right=0.980, hspace=0.35, wspace=0.225)

	plot_results(ax, os.path.join("output", "case_e_ve_vp_cr", "vtk", "displacement"), color_name="lightcoral", label_name="IMP")

	ax.set_xlabel("Time (hour)", size=12, fontname="serif")
	ax.set_ylabel("Total strain (%)", size=12, fontname="serif")
	ax.grid(True)
	ax.legend(loc=0, shadow=True, fancybox=True)

	apply_grey_theme(fig1, [ax], transparent=True)

	plt.show()

if __name__ == '__main__':
	main()