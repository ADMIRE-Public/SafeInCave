import numpy as np
import matplotlib.pyplot as plt
import json

minute = 60
hour = 60*minute
day = 24*hour
MPa = 1e6

def apply_grey_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.grid(True, color='0.92')
			ax.set_axisbelow(True)
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black')
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			ax.tick_params(axis='x', colors='black', which='both')
			ax.tick_params(axis='y', colors='black', which='both')
			ax.yaxis.label.set_color('black')
			ax.xaxis.label.set_color('black')
			ax.set_facecolor("0.85")

def save_json(data, file_name):
	with open(file_name, "w") as f:
	    json.dump(data, f, indent=4)

def get_pressure(schedule, p_0=0):
	time = [0]
	pressure = [p_0]

	t_schedule = schedule[:,[0, 2, 4, 6]]
	for i, dt in enumerate(t_schedule.flatten()):
		time.append(time[i] + dt)

	p_schedule = schedule[:,[1, 3, 5, 7]]
	for i, dp in enumerate(p_schedule.flatten()):
		p_i = pressure[i]
		pressure.append(p_i + dp)

	return np.array(time), -np.array(pressure)


def main_2():

	# schedule = np.array([
	#                     	[0*hour, 5*MPa, 5*MPa],
	#                     	[4*hour, 5*MPa, 10*MPa],
	#                     	[12*hour, 5*MPa, 10*MPa],
	#                     	[18*hour, 5*MPa, 8*MPa]
    # ])
	schedule = np.array([
	                    	[0*hour, 5*MPa, 6*MPa],
	                    	[2*hour, 5*MPa, 10*MPa],
	                    	[10*hour, 5*MPa, 10*MPa],
	                    	[12*hour, 5*MPa, 6*MPa],
	                    	[14*hour, 5*MPa, 6*MPa],
	                    	[16*hour, 5*MPa, 12*MPa],
	                    	[20*hour, 5*MPa, 12*MPa],
	                    	[22*hour, 5*MPa, 6*MPa],
	                    	[24*hour, 5*MPa, 6*MPa],
    ])
	time = schedule[:,0]
	sigma_radial = schedule[:,1]
	sigma_axial = schedule[:,2]

	# Build input_bc
	theta = 0.0
	input_bc = {}
	input_bc["Time"] = {"theta": theta, "timeList": list(time)}
	input_bc["sigma_radial"] = list(sigma_radial)
	input_bc["sigma_axial"] = list(sigma_axial)

	# Save input_bc
	save_json(input_bc, "input_bc_0.json")

	# Plot pressure schedule
	fig1, ax = plt.subplots(1, 1, figsize=(10, 4))
	fig1.subplots_adjust(top=0.945, bottom=0.115, left=0.065, right=0.980, hspace=0.35, wspace=0.225)

	ax.plot(time/hour, sigma_axial/MPa, "-", color="steelblue", label=r"$\sigma_{axial}$")
	ax.plot(time/hour, sigma_radial/MPa, "-", color="lightcoral", label=r"$\sigma_{radial}$")
	ax.set_xlabel("Time (days)", size=12, fontname="serif")
	ax.set_ylabel("Stress (MPa)", size=12, fontname="serif")
	ax.grid(True)

	apply_grey_theme(fig1, [ax], transparent=True)

	plt.show()

if __name__ == '__main__':
	main_2()