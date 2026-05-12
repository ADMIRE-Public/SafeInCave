import os
import numpy as np
import matplotlib.pyplot as plt
import safeincave as sf
import safeincave.PostProcessingTools as post
from safeincave.Utils import MPa, kPa, read_json
# import sys
# sys.path.append(os.path.join("solgeom"))
from terzaghi_solution import Solution
mm = 1e-3


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
			
def plot_pressure(ax, terza, fem_folder, t_indices=[1, 5, 10]):
	# Read FEM results
    points, time_list, p_field = post.read_node_scalar(os.path.join(fem_folder, "P", "P.xdmf"))
    mask = (points[:,0] == 1.0) & (points[:,1] == 1.0)
    z = points[mask,2]
    p = p_field[:,mask]/kPa
	
    z_exact = terza.getPositionValues(z)
	
    for i, tidx in enumerate(t_indices):
        # Calculate analytical solution
        print(time_list[tidx])
        p_exact = terza.getPressureValuesConstTime(time_list[tidx], ny=z_exact)
        p_exact /= kPa
        if i == 0:
            ax.plot(p[tidx,:], z, ".-", marker="x", markerfacecolor="none", color="steelblue", label="FEM")
            ax.plot(p_exact, z_exact, ".-", marker="o", markerfacecolor="none", color="lightcoral", label="Exact")
        else:
            ax.plot(p[tidx,:], z, ".-", marker="x", markerfacecolor="none", color="steelblue")
            ax.plot(p_exact, z_exact, ".-", marker="o", markerfacecolor="none", color="lightcoral")

    ax.set_xlabel("Pressure (kPa)", size=12, fontname="serif")
    ax.set_ylabel("z (m)", size=12, fontname="serif")
    ax.legend(loc=0, shadow=True, fancybox=True)
    ax.grid(True, color="0.92")
    ax.set_facecolor("0.85")
		

def plot_displacement(ax, terza, fem_folder, t_indices=[1, 5, 10]):
      # Read FEM results
    points, time_list, u_field = post.read_node_vector(os.path.join(fem_folder, "u", "u.xdmf"))
    mask = (points[:,0] == 1.0) & (points[:,1] == 1.0)
    z = points[mask,2]
    u = u_field[:,mask,:]/mm
	
    z_exact = terza.getPositionValues(z)
	
    for i, tidx in enumerate(t_indices):
        # Calculate analytical solution
        u_exact = terza.getDisplacementValuesConstTime(time_list[tidx], ny=z_exact)
        u_exact /= mm
        if i == 0:
            ax.plot(u[tidx,:,2], z, "-", marker="x", markerfacecolor="none", color="steelblue", label="FEM")
            ax.plot(u_exact, z_exact, "-", marker="o", markerfacecolor="none", color="lightcoral", label="Exact")
        else:
            ax.plot(u[tidx,:,2], z, "-", marker="x", markerfacecolor="none", color="steelblue")
            ax.plot(u_exact, z_exact, "-", marker="o", markerfacecolor="none", color="lightcoral")

    ax.set_xlabel("Displacement (mm)", size=12, fontname="serif")
    ax.set_ylabel("z (m)", size=12, fontname="serif")
    ax.legend(loc=0, shadow=True, fancybox=True)
    ax.grid(True, color="0.92")
    ax.set_facecolor("0.85")


	


def main():
    results_folder = os.path.join("output", "case_g_stab_1", "operation")
    # results_folder = os.path.join("output", "case_0")
    print(results_folder)

    props = read_json("props.json")
    solid_json = {"ROCK": props["ROCK"]}
    fluid_json = {"WATER": props["WATER"]}

    # Read analytical solution
    terza = Solution(height=1.0, tao_0=10*kPa, solid=solid_json, fluid=fluid_json, gravity=9.81)
    # terza = Solution(height=1.0, tao_0=10*kPa, solid=solid_json, fluid=fluid_json, gravity=0.0)

    # Plot pressure schedule
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))
    fig.subplots_adjust(top=0.970, bottom=0.135, left=0.076, right=0.980, hspace=0.35, wspace=0.225)
	
    # t_indices = [0, 2, 5, 10, 19, -1]
    t_indices = [0, 2, -1]
    # t_indices = [0, 2]
    # t_indices = [0]
    # t_indices = [1]
    plot_pressure(ax1, terza, results_folder, t_indices)
    plot_displacement(ax2, terza, results_folder, t_indices)


    # ax2.plot(p[t_idx,:], z, ".-", color="steelblue", label="FEM")
    # ax2.plot(p_exact, z_exact, ".-", color="lightcoral", label="Exact")
    # ax2.set_xlabel("Pressure (kPa)", size=10, fontname="serif")
    # ax2.set_ylabel("z (m)", size=10, fontname="serif")
    # ax2.legend(loc=0, shadow=True, fancybox=True)
    # ax2.grid(True, color="0.92")
    # ax2.set_facecolor("0.85")

    plt.show()


if __name__ == '__main__':
	main()