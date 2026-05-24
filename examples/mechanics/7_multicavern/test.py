import numpy as np
from safeincave.Utils import read_json, day, MPa
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

def main():
    data = read_json("input_cavern_data.json")

    t_water = np.array(data["Cavern_quarter"]["time"]) / day
    f_water = np.array(data["Cavern_quarter"]["flow"])

    t_hydrogen = np.array(data["Cavern_half"]["time"]) / day
    p_hydrogen = np.array(data["Cavern_half"]["P_hist"]) /MPa

    t_methane = np.array(data["Cavern_full"]["time"]) / day
    f_methane = np.array(data["Cavern_full"]["flow"])

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
    fig.subplots_adjust(top=0.935, bottom=0.155, left=0.119, right=0.980, hspace=0.35, wspace=0.38)
    fig.patch.set_alpha(0.0)

    ax1.plot(t_water, f_water, ".-", color="steelblue", label="Water")
    ax1.set_xlabel("Time (days)", fontname="serif", size=12)
    ax1.set_ylabel("Mass flow rate (kg/s)", fontname="serif", size=12)
	
    ax2.plot(t_methane, f_methane, ".-", color="forestgreen", label="Methane")
    ax2.set_xlabel("Time (days)", fontname="serif", size=12)
    ax2.set_ylabel("Mass flow rate (kg/s)", fontname="serif", size=12)
	
    ax3.plot(t_hydrogen, p_hydrogen, ".-", color="lightcoral", label="Hydrogen")
    ax3.set_xlabel("Time (days)", fontname="serif", size=12)
    ax3.set_ylabel("Pressure (MPa)", fontname="serif", size=12)
    
    apply_grey_theme(fig, [ax1, ax2, ax3], transparent=True)
	
    plt.show()

if __name__ == "__main__":
    main()