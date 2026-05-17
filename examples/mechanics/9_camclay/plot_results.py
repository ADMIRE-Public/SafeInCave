"""
Plot results from the cube-triaxial Modified Cam-Clay example.

Reads centroid-cell data from ./output/case_0/<formulation>/ and produces:
    - axial vs lateral total strain history
    - p-q stress path with the initial MCC yield ellipse + CSL overlay
    - preconsolidation pressure p_c evolution
    - yield function F evolution
"""

import safeincave.PostProcessingTools as post
import os
import numpy as np
import matplotlib.pyplot as plt


hour = 60 * 60
MPa  = 1e6


# MCC parameters used in main.py (kept in sync manually for plotting overlays)
M_CSL = 0.702
PC0_MPA = 13.0


def _styled(ax):
    ax.grid(True, color="0.92")
    ax.set_axisbelow(True)
    ax.set_facecolor("0.85")


def plot_eps_tot(ax, output_folder):
    centroids, time_list, eps_tot = post.read_cell_tensor(
        os.path.join(output_folder, "eps_tot", "eps_tot.xdmf"))

    target_point = [1.0, 1.0, 1.0]
    p_idx = post.find_closest_point(target_point, centroids)

    time_h = time_list / hour
    # safeincave convention: tension-positive, so compressive strain is negative.
    eps_axial   = -100 * eps_tot[:, p_idx, 2, 2]   # along loading direction (z)
    eps_lateral = -100 * eps_tot[:, p_idx, 0, 0]   # x lateral

    ax.plot(time_h, eps_axial,   "-", label=r"$\varepsilon_1$ (axial)")
    ax.plot(time_h, eps_lateral, "-", label=r"$\varepsilon_3$ (lateral)")
    ax.set_xlabel("Time (h)", size=12, fontname="serif")
    ax.set_ylabel("Total strain (%)", size=12, fontname="serif")
    ax.legend(loc=0, shadow=True, fancybox=True, prop={"size": 8})
    _styled(ax)


def plot_stress_path(ax, output_folder):
    centroids, time_list, p_elems = post.read_cell_scalar(
        os.path.join(output_folder, "p_elems", "p_elems.xdmf"))
    centroids, time_list, q_elems = post.read_cell_scalar(
        os.path.join(output_folder, "q_elems", "q_elems.xdmf"))

    target_point = [1.0, 1.0, 1.0]
    p_idx = post.find_closest_point(target_point, centroids)

    p_path = -p_elems[:, p_idx] / MPa    # compression-positive convention
    q_path =  q_elems[:, p_idx] / MPa

    ax.plot(p_path, q_path, "-", color="#377eb8", linewidth=1.6, label="Stress path")
    ax.scatter(p_path[0],  q_path[0],  c="white", edgecolors="black", zorder=5, label="t=0")
    ax.scatter(p_path[-1], q_path[-1], c="black", zorder=5, label="t=end")

    # Initial MCC ellipse and CSL
    p_ell = np.linspace(0.0, PC0_MPA, 200)
    q_ell = M_CSL * np.sqrt(np.clip(p_ell * (PC0_MPA - p_ell), 0.0, None))
    ax.plot(p_ell, q_ell, "-", color="black", linewidth=1.0,
            label=f"F=0 (p_c={PC0_MPA:.0f} MPa)")
    p_csl = np.linspace(0.0, PC0_MPA, 50)
    ax.plot(p_csl, M_CSL * p_csl, "--", color="dimgray", linewidth=0.8,
            label="CSL")

    ax.set_xlabel("Mean stress p (MPa)",      size=12, fontname="serif")
    ax.set_ylabel("Deviatoric stress q (MPa)", size=12, fontname="serif")
    ax.legend(loc=0, shadow=True, fancybox=True, prop={"size": 7})
    _styled(ax)


def plot_pc(ax, output_folder):
    centroids, time_list, pc = post.read_cell_scalar(
        os.path.join(output_folder, "pc", "pc.xdmf"))

    target_point = [1.0, 1.0, 1.0]
    p_idx = post.find_closest_point(target_point, centroids)

    time_h = time_list / hour
    pc_mpa = pc[:, p_idx] / MPa

    ax.plot(time_h, pc_mpa, "-", color="seagreen", linewidth=1.6)
    ax.axhline(PC0_MPA, color="black", linestyle="--", linewidth=0.7,
               label=f"p_c0 = {PC0_MPA:.0f} MPa")
    ax.set_xlabel("Time (h)", size=12, fontname="serif")
    ax.set_ylabel(r"$p_c$ (MPa)", size=12, fontname="serif")
    ax.legend(loc=0, shadow=True, fancybox=True, prop={"size": 8})
    _styled(ax)


def plot_Fvp(ax, output_folder):
    centroids, time_list, F = post.read_cell_scalar(
        os.path.join(output_folder, "Fvp", "Fvp.xdmf"))

    target_point = [1.0, 1.0, 1.0]
    p_idx = post.find_closest_point(target_point, centroids)

    time_h = time_list / hour
    F_val = F[:, p_idx] / (MPa * MPa)   # F has units of Pa^2; show MPa^2

    ax.plot(time_h, F_val, "-", color="tomato", linewidth=1.6)
    ax.axhline(0.0, color="black", linestyle="--", linewidth=0.7)
    ax.set_xlabel("Time (h)", size=12, fontname="serif")
    ax.set_ylabel(r"Yield function $F$ (MPa$^2$)", size=12, fontname="serif")
    _styled(ax)


def main():
    output_folder = os.path.join("output", "case_0", "P1P1")
    # output_folder = os.path.join("output", "case_0", "P1")
    # output_folder = os.path.join("output", "case_0", "P1P1_Stab")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 7))
    fig.subplots_adjust(top=0.96, bottom=0.09, left=0.075, right=0.985,
                        hspace=0.30, wspace=0.28)
    fig.patch.set_alpha(0.0)

    plot_eps_tot(ax1,     output_folder)
    plot_stress_path(ax2, output_folder)
    plot_pc(ax3,          output_folder)
    plot_Fvp(ax4,         output_folder)

    plt.show()


if __name__ == "__main__":
    main()
