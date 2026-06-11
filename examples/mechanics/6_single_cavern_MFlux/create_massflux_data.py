import numpy as np
from safeincave.Utils import hour
import safeincave as sf
from safeincave.CavernBC import Cavern_MassFlux
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


def get_massflux_history(scale=1.0):
    time_hist = (
        np.array(
            [
                -5.0,
                0.0,
                10,
                20,
                30,
                40,
                50,
                60,
                70,
                80,
                90,
                100,
                110,
                120,
                130,
                140,
                150,
                160,
                170,
                180,
                190,
            ]
        )
        * hour
        + 5 * hour
    )
    mr1 = 7e1
    mr2 = -7e1
    flow_hist = (
        np.array(
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
                mr1,
                mr1,
                0.0,
                0.0,
                mr2,
                mr2,
                0.0,
                0.0,
                1.1 * mr1,
                1.1 * mr1,
                0.0,
                0.0,
            ]
        )
        * scale
    )
    return time_hist, flow_hist


def main():
    # Read grid
    grid_path = os.path.join("..", "..", "..", "grids", "cavern_regular")
    grid = sf.GridHandlerGMSH("geom", grid_path)

    # Create mass flow rate history
    time_hist, flow_hist = get_massflux_history(15)

    fluid_name = "Methane"

    data = {"flow_rate": list(flow_hist), "time": list(time_hist)}
    sf.Utils.save_json(data, f"{fluid_name}.json")

    overburden = 10 * sf.Utils.MPa
    g = -9.81
    salt_density = 2200
    z_roof = 430
    z_floor = 205.189718
    z_ground = 660
    z_mid = (z_roof + z_floor) / 2
    p_mid = overburden + salt_density * abs(g) * (z_ground - z_mid)
    print(p_mid / sf.Utils.MPa)
    print(0.4 * p_mid / sf.Utils.MPa)

    cave = Cavern_MassFlux(
        grid=grid,
        cavern_name="Cavern",
        sym_scale=4,
        reference_point=[0.0, 0.0, z_mid],
        fluid="Methane",
        P_init=0.4 * p_mid,
        T_init=320.0,
        T_in=300.0,
        Mflux_values=flow_hist,
        time_values=time_hist,
        direction=2,
    )

    t = time_hist[0]
    t_final = time_hist[-1]
    dt = t_final / 50

    cave.calculate_volume()
    cave.calculate_initial_condition()
    cave.record_data(t)

    while t <= t_final:
        t += dt

        cave.update_cavern(t, dt)
        cave.record_data(t)
        cave.calculate_initial_condition()

    P_hist = np.array(cave.P_hist) / sf.Utils.MPa
    t_hist = np.array(cave.t_hist) / sf.Utils.day
    T_hist = np.array(cave.T_hist)
    density_hist = np.array(cave.density_hist)

    P_roof = P_hist + density_hist * abs(g) * (z_mid - z_roof) / sf.Utils.MPa
    P_floor = P_hist + density_hist * abs(g) * (z_mid - z_floor) / sf.Utils.MPa

    fig, axis = plt.subplots(2, 3, figsize=(14, 8))
    fig.subplots_adjust(
        top=0.935, bottom=0.155, left=0.119, right=0.980, hspace=0.35, wspace=0.38
    )
    fig.patch.set_alpha(0.0)

    axis[0, 0].plot(t_hist, P_hist, "-", color="steelblue", label="mid")
    axis[0, 0].plot(t_hist, P_floor, "-", color="0.3", label="floor")
    axis[0, 0].plot(t_hist, P_roof, "-", color="0.7", label="roof")
    axis[0, 0].set_ylabel("Pressure (MPa)", fontname="serif", size=12)
    axis[0, 0].set_xlabel("Time (day)", fontname="serif", size=12)
    axis[0, 0].grid(True)

    axis[1, 0].plot(t_hist, T_hist, "-")
    axis[1, 0].set_ylabel("Temperature (K)", fontname="serif", size=12)
    axis[1, 0].set_xlabel("Time (day)", fontname="serif", size=12)
    axis[1, 0].grid(True)

    axis[0, 1].plot(time_hist / sf.Utils.day, flow_hist, "-")
    axis[0, 1].set_ylabel("Flow rate (kg/m3)", fontname="serif", size=12)
    axis[0, 1].set_xlabel("Time (day)", fontname="serif", size=12)
    axis[0, 1].grid(True)

    print(P_hist.min(), P_hist.max())

    apply_grey_theme(fig, axis.flatten(), transparent=True)

    plt.show()


if __name__ == "__main__":
    main()
