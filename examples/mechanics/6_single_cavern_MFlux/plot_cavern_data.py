import safeincave as sf
import os
import numpy as np
import matplotlib.pyplot as plt

hour = 60*60
day = 24*hour
MPa = 1e6
ton = 1000
Mt = ton*1e6


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
			


def plot_cavern_data(axis, file_path, color, label):
    data = sf.Utils.read_json(file_path)

    time = np.array(data["Cavern"]["Data"]["Time"]) / sf.Utils.day
    pressure = np.array(data["Cavern"]["Data"]["Pressure"]) / sf.Utils.MPa
    temperature = np.array(data["Cavern"]["Data"]["Temperature"])
    volume = np.array(data["Cavern"]["Data"]["Volume"])
    volume_loss = 100*(volume[0] - volume) / volume[0]
    mass = np.array(data["Cavern"]["Data"]["Mass"]) / ton
    density = np.array(data["Cavern"]["Data"]["Density"])

    axis[0,0].plot(time, pressure, ".-", color=color, label=label)
    axis[0,0].set_xlabel("Time (days)", size=12, fontname="serif")
    axis[0,0].set_ylabel("Pressure (MPa)", size=12, fontname="serif")
    axis[0,0].grid(True)
    axis[0,0].legend(loc=0, shadow=True, fancybox=True)

    axis[1,0].plot(time, temperature, ".-", color=color, label=label)
    axis[1,0].set_xlabel("Time (days)", size=12, fontname="serif")
    axis[1,0].set_ylabel("Temperature (K)", size=12, fontname="serif")
    axis[1,0].grid(True)
    axis[1,0].legend(loc=0, shadow=True, fancybox=True)

    axis[0,1].plot(time, density, ".-", color=color, label=label)
    axis[0,1].set_xlabel("Time (days)", size=12, fontname="serif")
    axis[0,1].set_ylabel("Density (kg/m3)", size=12, fontname="serif")
    axis[0,1].grid(True)
    axis[0,1].legend(loc=0, shadow=True, fancybox=True)

    axis[0,2].plot(time, volume_loss, ".-", color=color, label=label)
    axis[0,2].set_xlabel("Time (days)", size=12, fontname="serif")
    axis[0,2].set_ylabel("Volume loss (%)", size=12, fontname="serif")
    axis[0,2].grid(True)
    axis[0,2].legend(loc=0, shadow=True, fancybox=True)

    axis[1,2].plot(time, mass, ".-", color=color, label=label)
    axis[1,2].set_xlabel("Time (days)", size=12, fontname="serif")
    axis[1,2].set_ylabel("Mass (ton)", size=12, fontname="serif")
    axis[1,2].grid(True)
    axis[1,2].legend(loc=0, shadow=True, fancybox=True)
	


def main_1():
    fig, axis = plt.subplots(2, 3, figsize=(14, 8))
    fig.subplots_adjust(top=0.935, bottom=0.155, left=0.119, right=0.980, hspace=0.35, wspace=0.38)
    fig.patch.set_alpha(0.0)

    file_path = os.path.join("output", "case_0", "cavern_data.json")
    plot_cavern_data(axis, file_path, color="steelblue", label="Methane")

    apply_grey_theme(fig, axis.flatten(), transparent=True)
    plt.show()

if __name__ == "__main__":
	main_1()