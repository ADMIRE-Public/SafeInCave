import numpy as np
from safeincave.Utils import day, MPa
import safeincave as sf
import matplotlib.pyplot as plt


def apply_grey_theme(fig, axes, transparent=True, grid_color="0.92", back_color='0.85'):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax is not None:
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



def get_massflux_hydrogen():
    period = 312
    t_hist = np.linspace(0, period, 1000)
    P_hist = 9 + 6 * np.sin(2 * np.pi * 20 * t_hist / period)
    T_hist = 330 + 20 * np.sin(2 * np.pi * 20 * t_hist / period)
    return list(t_hist*day), list(P_hist*MPa), list(T_hist)


def get_PT_hydrogen():
    data = np.loadtxt("PT_H2.csv", delimiter=",", skiprows=1)
    t_hist = data[:, 0]*day
    P_hist = (data[:, 1] - 5)*MPa
    T_hist = data[:, 2]
    return list(t_hist), list(P_hist), list(T_hist)



def get_massflux_methane(scale=1.0):
    time_hist = np.array([
        0, 5, 15, 25, 35, 
        45, 55, 65, 75,
        85, 95, 105, 115,
        125, 135, 145, 155,
    ])*day
    mr1 = 7e1
    mr2 = -7e1
    flow_hist = np.array(
        [
            0.0, mr1, mr1, 0.0, 0.0,
            1.1*mr2, 1.1*mr2, 0.0, 0.0,
            mr1, mr1, 0.0, 0.0,
            mr2, mr2, 0.0, 0.0,
        ]
    ) * scale / 25
    time_ext = np.concatenate((time_hist, time_hist + time_hist[-1]))
    flow_ext = np.concatenate((flow_hist, flow_hist))

    return list(time_ext), list(flow_ext)


def main():

    # Initialize data
    data = {}
    data["Cavern_full"] = {"fluid": "Methane"}
    data["Cavern_half"] = {"fluid": "Hydrogen"}
    data["Cavern_quarter"] = {"fluid": "Water"}

    # Create Methane mass flow rate history
    time_hist, flow_hist = get_massflux_methane(15)
    data["Cavern_full"]["P_init"] = 6.95*MPa
    data["Cavern_full"]["T_init"] = 320.0
    data["Cavern_full"]["time"] = [x.item() for x in time_hist]
    data["Cavern_full"]["flow"] = [x.item() for x in flow_hist]

    # Create Hydrogen pressure history
    # t_hist, P_hist, T_hist = get_massflux_hydrogen()
    t_hist, P_hist, T_hist = get_PT_hydrogen()
    data["Cavern_half"]["time"] = [x.item() for x in t_hist]
    data["Cavern_half"]["P_hist"] = P_hist
    data["Cavern_half"]["T_hist"] = T_hist

    # Create Water mass flow rate history
    data["Cavern_quarter"]["P_init"] = 6.95*MPa
    data["Cavern_quarter"]["T_init"] = 320.0
    data["Cavern_quarter"]["time"] = [0.0, 312*day]
    data["Cavern_quarter"]["flow"] = [0.0, 0.0]

    sf.Utils.save_json(data, "input_cavern_data.json")


    fig, axis = plt.subplots(2, 3, figsize=(14, 8))
    fig.subplots_adjust(top=0.935, bottom=0.155, left=0.119, right=0.980, hspace=0.35, wspace=0.38)
    fig.patch.set_alpha(0.0)

    axis[0,0].plot(np.array(data["Cavern_half"]["time"])/day, np.array(data["Cavern_half"]["P_hist"])/MPa, "-", color="steelblue", linewidth=2.0)
    axis[0,0].set_xlabel("Time (days)", fontname="serif", size=12)
    axis[0,0].set_ylabel("Pressure (MPa)", fontname="serif", size=12)
    axis[0,0].set_title("Cavern_half - Hydrogen", fontname="serif", size=14)


    axis[1,0].plot(np.array(data["Cavern_half"]["time"])/day, np.array(data["Cavern_half"]["T_hist"]), "-", color="steelblue", linewidth=2.0)
    axis[1,0].set_xlabel("Time (days)", fontname="serif", size=12)
    axis[1,0].set_ylabel("Temperature (K)", fontname="serif", size=12)
    axis[1,0].set_title("Cavern_half - Hydrogen", fontname="serif", size=14)

    axis[0,1].plot(np.array(data["Cavern_full"]["time"])/day, np.array(data["Cavern_full"]["flow"]), "-", color="forestgreen", linewidth=2.0)
    axis[0,1].set_xlabel("Time (days)", fontname="serif", size=12)
    axis[0,1].set_ylabel("Mass flow rate (kg/s)", fontname="serif", size=12)
    axis[0,1].set_title("Cavern_full - Methane", fontname="serif", size=14)

    axis[0,2].plot(np.array(data["Cavern_quarter"]["time"])/day, np.array(data["Cavern_quarter"]["flow"]), "-", color="lightcoral", linewidth=2.0)
    axis[0,2].set_xlabel("Time (days)", fontname="serif", size=12)
    axis[0,2].set_ylabel("Mass flow rate (kg/s)", fontname="serif", size=12)
    axis[0,2].set_title("Cavern_quarter - Water", fontname="serif", size=14)


    apply_grey_theme(fig, axis.flatten(), transparent=True)

    plt.show()



    # fluid_name = "Methane"

    # data = {
    #       "flow_rate": list(flow_hist),
    #       "time": list(time_hist)
    # }
    # sf.Utils.save_json(data, f"{fluid_name}.json")
    
    # overburden = 10*sf.Utils.MPa
    # g = -9.81
    # salt_density = 2200
    # z_roof = 420
    # z_floor = 200
    # z_ground = 660
    # z_mid = (z_roof + z_floor) / 2
    # p_mid = overburden + salt_density*abs(g)*(z_ground - z_mid)
    # print(p_mid/sf.Utils.MPa)

    # cave = Cavern_MassFlux(
    #     grid = grid,
    #     cavern_name = "Cavern",
    #     sym_scale = 4,
    #     reference_point = [0.0, 0.0, z_mid],
    #     fluid = "Methane",
    #     P_init = 0.4*p_mid,
    #     T_init = 320.0,
    #     T_in = 300.0,
    #     Q_in = 0.0,
    #     Mflux_values = flow_hist,
    #     time_values = time_hist,
    #     direction = 2
    # )

    # t = time_hist[0]
    # t_final = time_hist[-1]
    # dt = t_final / 50

    # cave.calculate_volume()
    # cave.calculate_initial_condition()
    # cave.record_data(t)

    # while t <= t_final:

    #     t += dt

    #     cave.update_cavern(t, dt)
    #     cave.record_data(t)
    #     cave.calculate_initial_condition()


    # P_hist = np.array(cave.P_hist) / sf.Utils.MPa
    # V_hist = np.array(cave.V_hist)
    # M_hist = np.array(cave.M_hist)
    # t_hist = np.array(cave.t_hist) / sf.Utils.day
    # T_hist = np.array(cave.T_hist)
    # density_hist = np.array(cave.density_hist)

    # P_roof = P_hist + density_hist*abs(g)*(z_mid - z_roof) / sf.Utils.MPa
    # P_floor = P_hist + density_hist*abs(g)*(z_mid - z_floor) / sf.Utils.MPa


    # fig, axis = plt.subplots(2, 3, figsize=(14, 8))
    # fig.subplots_adjust(top=0.935, bottom=0.155, left=0.119, right=0.980, hspace=0.35, wspace=0.38)
    # fig.patch.set_alpha(0.0)
	
    # axis[0,0].plot(time_hist/day, flow_hist, "-", color="steelblue", label="mid")

    # axis[0,0].plot(t_hist, P_hist, "-", color="steelblue", label="mid")
    # axis[0,0].plot(t_hist, P_floor, "-", color="0.3", label="floor")
    # axis[0,0].plot(t_hist, P_roof, "-", color="0.7", label="roof")
    # axis[0,0].set_ylabel("Pressure (MPa)", fontname="serif", size=12)
    # axis[0,0].set_xlabel("Time (day)", fontname="serif", size=12)
    # axis[0,0].grid(True)

    # axis[1,0].plot(t_hist, T_hist, "-")
    # axis[1,0].set_ylabel("Temperature (K)", fontname="serif", size=12)
    # axis[1,0].set_xlabel("Time (day)", fontname="serif", size=12)
    # axis[1,0].grid(True)

    # print(P_hist.min(), P_hist.max())

    # apply_grey_theme(fig, axis.flatten(), transparent=True)

    plt.show()




if __name__ == "__main__":
    main()