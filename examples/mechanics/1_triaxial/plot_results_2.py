import os
import sys
sys.path.append(os.path.join("..", "..", "..", "safeincave"))
from PostProcessingTools_2 import read_cell_scalar, read_cell_tensor, find_closest_point
import numpy as np
import matplotlib.pyplot as plt

hour = 60*60
MPa = 1e6

def plot_eps_tot(ax, output_folder):
	centroids, time_list, eps_tot = read_cell_tensor(os.path.join(output_folder, "eps_tot", "eps_tot.xdmf"))

	target_point = [1.0, 1.0, 1.0]
	p_idx = find_closest_point(target_point, centroids)

	time_list /= hour
	eps_3 = -100*eps_tot[:,p_idx,0,0]
	eps_1 = -100*eps_tot[:,p_idx,2,2]

	ax.plot(time_list, eps_1, "-", label=r"$\varepsilon_1$")
	ax.plot(time_list, eps_3, "-", label=r"$\varepsilon_3$")
	ax.set_xlabel("Time (h)", size=12, fontname="serif")
	ax.set_ylabel("Total strain (%)", size=12, fontname="serif")
	ax.legend(loc=0, shadow=True, fancybox=True, prop={"size": 8})
	ax.grid(True, color="0.92")
	ax.set_facecolor("0.85")

def plot_strains(ax, output_folder):
	centroids, time_list, eps_ve = read_cell_tensor(os.path.join(output_folder, "eps_ve", "eps_ve.xdmf"))
	centroids, time_list, eps_cr = read_cell_tensor(os.path.join(output_folder, "eps_cr", "eps_cr.xdmf"))
	centroids, time_list, eps_vp = read_cell_tensor(os.path.join(output_folder, "eps_vp", "eps_vp.xdmf"))
	centroids, time_list, eps_tot = read_cell_tensor(os.path.join(output_folder, "eps_tot", "eps_tot.xdmf"))

	eps_e = eps_tot - eps_ve - eps_cr - eps_vp
	time_list /= hour

	target_point = [1.0, 1.0, 1.0]
	p_idx = find_closest_point(target_point, centroids)

	eps_tot = -100*(eps_tot[:,p_idx,2,2] - eps_tot[:,p_idx,0,0])
	eps_ve = -100*(eps_ve[:,p_idx,2,2] - eps_ve[:,p_idx,0,0])
	eps_cr = -100*(eps_cr[:,p_idx,2,2] - eps_cr[:,p_idx,0,0])
	eps_vp = -100*(eps_vp[:,p_idx,2,2] - eps_vp[:,p_idx,0,0])
	eps_e = -100*(eps_e[:,p_idx,2,2] - eps_e[:,p_idx,0,0])

	ax.plot(time_list, eps_tot, "-", label=r"$\varepsilon_\mathrm{tot}$")
	ax.plot(time_list, eps_ve, "-", label=r"$\varepsilon_\mathrm{ve}$")
	ax.plot(time_list, eps_cr, "-", label=r"$\varepsilon_\mathrm{cr}$")
	ax.plot(time_list, eps_vp, "-", label=r"$\varepsilon_\mathrm{vp}$")
	ax.plot(time_list, eps_e, "-", label=r"$\varepsilon_\mathrm{e}$")
	ax.set_xlabel("Time (h)", size=12, fontname="serif")
	ax.set_ylabel(r"$\varepsilon_1 - \varepsilon_3$ (%)", size=12, fontname="serif")
	ax.legend(loc=0, shadow=True, fancybox=True, prop={"size": 8})
	ax.grid(True, color="0.92")
	ax.set_facecolor("0.85")

def plot_Fvp(ax, output_folder):
	centroids, time_list, Fvp = read_cell_scalar(os.path.join(output_folder, "Fvp", "Fvp.xdmf"))

	time_list /= hour
	target_point = [1.0, 1.0, 1.0]
	p_idx = find_closest_point(target_point, centroids)
	Fvp = Fvp[:,p_idx]

	ax.plot(time_list, Fvp, "-", label=r"$\varepsilon_\mathrm{e}$")
	ax.plot(time_list, len(time_list)*[0], "--", color="tomato")
	ax.set_xlabel("Time (h)", size=12, fontname="serif")
	ax.set_ylabel(r"Yield function, $F_{vp}$", size=12, fontname="serif")
	ax.grid(True, color="0.92")
	ax.set_facecolor("0.85")




def main():
	# Define paths
	output_folder = os.path.join("output", "case_0")

	# Plot loading schedule
	fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 3))
	fig.subplots_adjust(top=0.970, bottom=0.155, left=0.062, right=0.980, hspace=0.35, wspace=0.312)
	fig.patch.set_alpha(0.0)

	plot_eps_tot(ax1, output_folder)
	plot_strains(ax2, output_folder)
	plot_Fvp(ax3, output_folder)

	

	plt.show()


if __name__ == '__main__':
	main()