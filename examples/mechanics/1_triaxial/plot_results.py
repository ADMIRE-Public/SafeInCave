import safeincave as sf
import safeincave.PostProcessingTools as post
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

hour = 60*60
MPa = 1e6

def find_closest_point(target_point: list, points: pd.DataFrame) -> int:
	x_p, y_p, z_p = target_point
	d = np.sqrt(  (points["x"].values - x_p)**2
	            + (points["y"].values - y_p)**2
	            + (points["z"].values - z_p)**2 )
	cell_p = d.argmin()
	return cell_p


def plot_strains(ax, output_folder):
	points_xdmf, cells_xdmf = post.read_xdmf_as_pandas(os.path.join(output_folder, "eps_tot", "eps_tot.xdmf"))
	mid_cells = post.compute_cell_centroids(points_xdmf.values, cells_xdmf.values)

	target_point = [0.5, 0.5, 0.5]
	cell_id = find_closest_point(target_point, mid_cells)

	ve_xx, ve_yy, ve_zz, ve_xy, ve_xz, ve_yz = post.read_tensor_from_cells(os.path.join(output_folder, "eps_ve", "eps_ve.xdmf"))
	cr_xx, cr_yy, cr_zz, cr_xy, cr_xz, cr_yz = post.read_tensor_from_cells(os.path.join(output_folder, "eps_cr", "eps_cr.xdmf"))
	vp_xx, vp_yy, vp_zz, vp_xy, vp_xz, vp_yz = post.read_tensor_from_cells(os.path.join(output_folder, "eps_vp", "eps_vp.xdmf"))
	tot_xx, tot_yy, tot_zz, tot_xy, tot_xz, tot_yz = post.read_tensor_from_cells(os.path.join(output_folder, "eps_tot", "eps_tot.xdmf"))

	eps_ve = 100*(ve_xx - ve_zz).iloc[cell_id].values
	eps_cr = 100*(cr_xx - cr_zz).iloc[cell_id].values
	eps_vp = 100*(vp_xx - vp_zz).iloc[cell_id].values
	eps_tot = 100*(tot_xx - tot_zz).iloc[cell_id].values

	t = tot_xx.iloc[cell_id].index.values/hour

	ax.plot(t, eps_tot, "-", label=r"$\varepsilon_\mathrm{tot}$")
	ax.plot(t, eps_ve, "-", label=r"$\varepsilon_\mathrm{ve}$")
	ax.plot(t, eps_cr, "-", label=r"$\varepsilon_\mathrm{cr}$")
	ax.plot(t, eps_vp, "-", label=r"$\varepsilon_\mathrm{vp}$")
	ax.set_xlabel("Time (h)", size=12, fontname="serif")
	ax.set_ylabel(r"$\varepsilon_1 - \varepsilon_3$ (%)", size=12, fontname="serif")
	ax.legend(loc=0, shadow=True, fancybox=True, prop={"size": 8})
	ax.grid(True, color="0.92")
	ax.set_facecolor("0.85")



def plot_eps_tot(ax, output_folder):
	points_xdmf, cells_xdmf = post.read_xdmf_as_pandas(os.path.join(output_folder, "eps_tot", "eps_tot.xdmf"))
	mid_cells = post.compute_cell_centroids(points_xdmf.values, cells_xdmf.values)

	# target_point = [0.5, 0.5, 0.5]
	target_point = [1.0, 1.0, 1.0]
	cell_id = find_closest_point(target_point, mid_cells)

	tot_xx, tot_yy, tot_zz, tot_xy, tot_xz, tot_yz = post.read_tensor_from_cells(os.path.join(output_folder, "eps_tot", "eps_tot.xdmf"))

	# eps_tot = 100*(tot_xx - tot_zz).iloc[cell_id].values
	eps_1 = -100*tot_zz.iloc[cell_id].values
	eps_3 = -100*tot_xx.iloc[cell_id].values

	t = tot_xx.iloc[cell_id].index.values/hour

	ax.plot(t, eps_1, "-", label=r"$\varepsilon_1$")
	ax.plot(t, eps_3, "-", label=r"$\varepsilon_3$")
	ax.set_xlabel("Time (h)", size=12, fontname="serif")
	ax.set_ylabel("Total strain (%)", size=12, fontname="serif")
	ax.legend(loc=0, shadow=True, fancybox=True, prop={"size": 8})
	ax.grid(True, color="0.92")
	ax.set_facecolor("0.85")






def plot_Fvp(ax, output_folder):
	points_xdmf, cells_xdmf = post.read_xdmf_as_pandas(os.path.join(output_folder, "Fvp", "Fvp.xdmf"))
	mid_cells = post.compute_cell_centroids(points_xdmf.values, cells_xdmf.values)

	# target_point = [0.5, 0.5, 0.5]
	target_point = [1.0, 1.0, 1.0]
	cell_id = find_closest_point(target_point, mid_cells)

	df_Fvp = post.read_scalar_from_cells(os.path.join(output_folder, "Fvp", "Fvp.xdmf"))
	Fvp = df_Fvp.iloc[cell_id].values

	t = df_Fvp.iloc[cell_id].index.values/hour

	ax.plot(t, Fvp, "-")
	ax.plot(t, len(t)*[0], "--", color="tomato")
	ax.set_xlabel("Time (h)", size=12, fontname="serif")
	ax.set_ylabel("Yield function (-)", size=12, fontname="serif")
	ax.grid(True, color="0.92")
	ax.set_facecolor("0.85")


def main():
	results_folder = os.path.join("output", "case_0")
	# results_folder = os.path.join("output", "test")

	# Plot loading schedule
	fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 3))
	fig.subplots_adjust(top=0.970, bottom=0.155, left=0.062, right=0.980, hspace=0.35, wspace=0.312)
	fig.patch.set_alpha(0.0)

	plot_eps_tot(ax1, results_folder)
	plot_strains(ax2, results_folder)
	plot_Fvp(ax3, results_folder)

	plt.show()


if __name__ == '__main__':
	main()