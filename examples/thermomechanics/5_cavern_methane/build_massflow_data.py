import numpy as np
from safeincave.Utils import day, year
import safeincave as sf
import matplotlib.pyplot as plt


def apply_grey_theme(fig, axes, transparent=True, grid_color="0.92", back_color="0.85"):
    fig.patch.set_facecolor("#212121ff")
    if transparent:
        fig.patch.set_alpha(0.0)
    for ax in axes:
        if ax is None:
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


def get_massflux_methane(scale=1.0):
    time_hist = (
        np.array(
            [
                0,
                5,
                15,
                25,
                35,
                45,
                55,
                65,
                75,
                85,
                95,
                105,
                115,
                125,
                135,
                145,
                155,
            ]
        )
        * day
    )
    mr1 = 30
    mr2 = -30
    flow_hist = np.array(
        [
            0.0,
            mr1,
            mr1,
            0.0,
            0.0,
            1.1 * mr2,
            1.1 * mr2,
            0.0,
            0.0,
            1.1 * mr1,
            1.1 * mr1,
            0.0,
            0.0,
            mr2,
            mr2,
            0.0,
            0.0,
        ]
    )
    time_hist = np.concatenate((time_hist, time_hist + time_hist[-1]))
    time_hist = np.concatenate((time_hist, time_hist + time_hist[-1]))
    flow_hist = np.concatenate((flow_hist, flow_hist))
    flow_hist = np.concatenate((flow_hist, flow_hist))

    # Limit it to 1 year
    mask = time_hist < 1 * year
    time_hist = time_hist[mask]
    flow_hist = flow_hist[mask]

    return time_hist, flow_hist


def main():

    time, flux = get_massflux_methane()

    # Initialize data
    data = {}
    data["time"] = [x.item() for x in time]
    data["flow"] = [x.item() for x in flux]
    sf.Utils.save_json(data, "input_massflux.json")

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    fig.subplots_adjust(
        top=0.935, bottom=0.155, left=0.138, right=0.980, hspace=0.35, wspace=0.38
    )
    fig.patch.set_alpha(0.0)

    ax.plot(
        np.array(data["time"]) / day,
        np.array(data["flow"]),
        "-",
        color="forestgreen",
        linewidth=2.0,
    )
    ax.set_xlabel("Time (days)", fontname="serif", size=12)
    ax.set_ylabel("Mass flow rate (kg/s)", fontname="serif", size=12)

    apply_grey_theme(fig, [ax], transparent=True)

    plt.show()


if __name__ == "__main__":
    main()
