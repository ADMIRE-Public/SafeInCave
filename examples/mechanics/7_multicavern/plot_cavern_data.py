import numpy as np
from safeincave.Utils import read_json, day, MPa
import matplotlib.pyplot as plt
import os


def apply_grey_theme(fig, axes, transparent=True, grid_color="0.92", back_color="0.85"):
    fig.patch.set_facecolor("#212121ff")
    if transparent:
        fig.patch.set_alpha(0.0)
    for ax in axes:
        if ax is not None:
            ax.grid(True, color=grid_color)
            ax.set_axisbelow(True)
            ax.spines["bottom"].set_color("black")
            ax.spines["top"].set_color("black")
            ax.spines["right"].set_color("black")
            ax.spines["left"].set_color("black")
            ax.tick_params(axis="x", colors="black", which="both")
            ax.tick_params(axis="y", colors="black", which="both")
            ax.yaxis.label.set_color("black")
            ax.xaxis.label.set_color("black")
            ax.set_facecolor(back_color)


def main():
    file_name = os.path.join("output", "case_0", "cavern_data.json")
    data = read_json(file_name)

    ton = 1e3
    Mton = ton * 1e6

    t_water = np.array(data["Cavern_quarter"]["Data"]["Time"]) / day
    p_water = np.array(data["Cavern_quarter"]["Data"]["Pressure"]) / MPa
    v_water = np.array(data["Cavern_quarter"]["Data"]["Volume"])
    T_water = np.array(data["Cavern_quarter"]["Data"]["Temperature"]) - 273
    v_loss_water = 100 * (v_water[0] - v_water) / v_water[0]
    m_water = np.array(data["Cavern_quarter"]["Data"]["Mass"]) / Mton

    t_hydrogen = np.array(data["Cavern_half"]["Data"]["Time"]) / day
    p_hydrogen = np.array(data["Cavern_half"]["Data"]["Pressure"]) / MPa
    v_hydrogen = np.array(data["Cavern_half"]["Data"]["Volume"])
    T_hydrogen = np.array(data["Cavern_half"]["Data"]["Temperature"]) - 273
    v_loss_hydrogen = 100 * (v_hydrogen[0] - v_hydrogen) / v_hydrogen[0]
    m_hydrogen = np.array(data["Cavern_half"]["Data"]["Mass"]) / Mton

    t_methane = np.array(data["Cavern_full"]["Data"]["Time"]) / day
    p_methane = np.array(data["Cavern_full"]["Data"]["Pressure"]) / MPa
    v_methane = np.array(data["Cavern_full"]["Data"]["Volume"])
    T_methane = np.array(data["Cavern_full"]["Data"]["Temperature"]) - 273
    v_loss_methane = 100 * (v_methane[0] - v_methane) / v_methane[0]
    m_methane = np.array(data["Cavern_full"]["Data"]["Mass"]) / Mton

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    fig.subplots_adjust(
        top=0.935, bottom=0.155, left=0.057, right=0.980, hspace=0.35, wspace=0.38
    )
    fig.patch.set_alpha(0.0)

    ax1 = ax[0, 0]
    ax2 = ax[0, 1]
    ax3 = ax[1, 0]
    ax4 = ax[1, 1]

    ax1.plot(t_water, p_water, "-", color="lightcoral", label="Water")
    ax1.plot(t_methane, p_methane, "-", color="forestgreen", label="Methane")
    ax1.plot(t_hydrogen, p_hydrogen, "-", color="steelblue", label="Hydrogen")
    ax1.set_xlabel("Time (days)", fontname="serif", size=12)
    ax1.set_ylabel("Pressure (MPa)", fontname="serif", size=12)

    ax2.plot(t_water, v_loss_water, "-", color="lightcoral", label="Water")
    ax2.plot(t_methane, v_loss_methane, "-", color="forestgreen", label="Methane")
    ax2.plot(t_hydrogen, v_loss_hydrogen, "-", color="steelblue", label="Hydrogen")
    ax2.set_xlabel("Time (days)", fontname="serif", size=12)
    ax2.set_ylabel("Volume loss (%)", fontname="serif", size=12)

    ax3.semilogy(t_water, m_water, "-", color="lightcoral", label="Water")
    ax3.semilogy(t_methane, m_methane, "-", color="forestgreen", label="Methane")
    ax3.semilogy(t_hydrogen, m_hydrogen, "-", color="steelblue", label="Hydrogen")
    ax3.set_xlabel("Time (days)", fontname="serif", size=12)
    ax3.set_ylabel("Mass (Mton)", fontname="serif", size=12)

    ax4.plot(t_water, T_water, "-", color="lightcoral", label="Water")
    ax4.plot(t_methane, T_methane, "-", color="forestgreen", label="Methane")
    ax4.plot(t_hydrogen, T_hydrogen, "-", color="steelblue", label="Hydrogen")
    ax4.set_xlabel("Time (days)", fontname="serif", size=12)
    ax4.set_ylabel("Temperature (°C)", fontname="serif", size=12)

    apply_grey_theme(fig, [ax1, ax2, ax3, ax4], transparent=True)

    plt.show()


if __name__ == "__main__":
    main()
