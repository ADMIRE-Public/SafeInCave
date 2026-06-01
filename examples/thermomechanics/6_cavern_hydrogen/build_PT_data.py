import numpy as np
from safeincave.Utils import day, MPa
import safeincave as sf
import matplotlib.pyplot as plt


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


def get_PT_hydrogen():
    data = np.loadtxt("PT_H2.csv", delimiter=",", skiprows=1)
    t_hist = data[:, 0] * day
    P_hist = (data[:, 1] - 5) * MPa
    T_hist = data[:, 2]
    return list(t_hist), list(P_hist), list(T_hist)


def main():
    # Read data
    t_hist, P_hist, T_hist = get_PT_hydrogen()

    # Initialize data
    data = {}
    data["time"] = t_hist
    data["pressure"] = P_hist
    data["temperature"] = T_hist

    sf.Utils.save_json(data, "input_PT.json")

    fig, axis = plt.subplots(1, 2, figsize=(10, 4))
    fig.subplots_adjust(
        top=0.935, bottom=0.155, left=0.064, right=0.980, hspace=0.35, wspace=0.252
    )
    fig.patch.set_alpha(0.0)

    axis[0].plot(
        np.array(data["time"]) / day,
        np.array(data["pressure"]) / MPa,
        "-",
        color="steelblue",
        linewidth=2.0,
    )
    axis[0].set_xlabel("Time (days)", fontname="serif", size=12)
    axis[0].set_ylabel("Pressure (MPa)", fontname="serif", size=12)

    axis[1].plot(
        np.array(data["time"]) / day,
        np.array(data["temperature"]) - 273,
        "-",
        color="steelblue",
        linewidth=2.0,
    )
    axis[1].set_xlabel("Time (days)", fontname="serif", size=12)
    axis[1].set_ylabel("Temperature (°C)", fontname="serif", size=12)

    apply_grey_theme(fig, axis.flatten(), transparent=True)

    plt.show()


if __name__ == "__main__":
    main()
