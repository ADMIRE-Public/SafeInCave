import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(os.path.join("..", "..", "libs"))
from ResultsHandler import ResultsReader
# from sls_maxwell import sls_maxwell

sec = 1.
minute = 60*sec
sec = 1.
hour = 60*minute
day = 24*hour
kPa = 1e3
MPa = 1e6
GPa = 1e9

def apply_white_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		ax.grid(True, color='#c4c4c4ff')
		ax.set_axisbelow(True)
		ax.spines['bottom'].set_color('black')
		ax.spines['top'].set_color('black')
		ax.spines['right'].set_color('black')
		ax.spines['left'].set_color('black')
		ax.tick_params(axis='x', colors='black', which='both')
		ax.tick_params(axis='y', colors='black', which='both')
		ax.yaxis.label.set_color('black')
		ax.xaxis.label.set_color('black')
		ax.set_facecolor("#e9e9e9ff")

def get_displacement(res, fun_point):
	# Axial strain
	res.load_points_of_interest(fun_point)
	x, y, z, t, field = res.get_results_over_all_times("Displacement")
	t = np.array(t)/day
	u = field[:,:,0]
	v = field[:,:,1]
	w = field[:,:,2]
	return t, x, y, z, u, v, w

def plot_initial_cavern(ax1, x, z):
	beta = 100
	time_index = -1
	ax1.plot(x, z, "--", color="0.0", label="Initial shape")

def plot_cavern(ax1, x, z, u, w, color_name, label_name):
	beta = 10
	time_index = -1
	dx = beta*u[time_index,:]
	dz = beta*w[time_index,:]

	ax1.plot(x+dx, z+dz, "-", color=color_name, label=label_name)

	ax1.set_xlabel("Coordinate X (m)", size=12, fontname="serif")
	ax1.set_ylabel("Coordinate Z (m)", size=12, fontname="serif")
	ax1.set_xlim(-50, 300)
	ax1.set_ylim(-477, -180)
	ax1.text(37, -475, f"*Amplification factor: {beta}x", size=10, fontname="serif", color="0.4")
	ax1.legend(loc=0, fancybox=True, shadow=True)

def sort_points(z, x, y, u, v, w):
	indices = sorted(range(len(z)), key=z.__getitem__)
	x = np.array([x[i] for i in indices])
	y = np.array([y[i] for i in indices])
	z = np.array([z[i] for i in indices])
	u = np.array([u[:,i] for i in indices]).transpose()
	v = np.array([v[:,i] for i in indices]).transpose()
	w = np.array([w[:,i] for i in indices]).transpose()
	return x, y, z, u, v, w

def trapezoidal_volume(x, y):
	"""
	This function calculates the volume of a solid of revolution (around y=0 axis) based on the trapezoidal rule.
	"""
	volume = 0.0
	n = len(x)
	for i in range(1, n):
		R = 0.5*(y[i] + y[i-1])
		A = np.pi*R**2
		d = x[i] - x[i-1]
		volume += A*d
	return volume

def wall(x, y, z):
	alpha = -1.
	R = 45*alpha
	x_0, z_0 = 0.0, 245*alpha
	x_1, z_1 = 0.0, 415*alpha
	if np.isclose(y, 0.0, atol=1e-12, rtol=0):
		x_wall = x_0 + np.sqrt((R**2 - (z - z_0)**2))
		if abs(z) < abs(z_0-R) and np.isclose(x, 0.0, atol=1e-12, rtol=0):
			return False
		elif abs(z_0-R) <= abs(z) and abs(z) <= abs(z_0) and np.isclose(x, x_0 + np.sqrt((R**2 - (z - z_0)**2)), atol=1e-12, rtol=0):
			return True
		elif abs(z_0) <= abs(z) and abs(z) <= abs(z_1) and np.isclose(x, abs(R), atol=1e-12, rtol=0):
			return True
		elif abs(z_1) <= abs(z) and abs(z) <= abs(z_1+R) and np.isclose(x, x_1 + np.sqrt((R**2 - (z - z_1)**2)), atol=1e-12, rtol=0):
			return True
		elif abs(z) > abs(z_1-R) and np.isclose(x, 0.0, atol=1e-12, rtol=0):
			return False
		else:
			return False
	else:
		return False

def plot_res(ax1, ax2, output_folder, color_name, label_name, plot_initial=False):
	res_A = ResultsReader(output_folder, pvd_file_index=1)
	t_day, x, y, z, u, v, w = get_displacement(res_A, wall)
	x, y, z, u, v, w = sort_points(z, x, y, u, v, w)

	if plot_initial:
		plot_initial_cavern(ax1, x, z)
	plot_cavern(ax1, x, z, u, w, color_name, label_name)

	vol_i = trapezoidal_volume(z, x)
	vol_f = trapezoidal_volume(z+w[-1,:], x+u[-1,:])
	print(f"Initial volume: {vol_i}")
	print(f"Final volume: {vol_f}")
	volumes = []
	for t_index in range(len(t_day)):
		vol = vol_i - trapezoidal_volume(z+w[t_index,:], x+u[t_index,:])
		volumes.append(vol/vol_i)
	volumes = np.array(volumes)

	ax2.plot(t_day, volumes*100, "-", color=color_name)
	ax2.set_xlabel("Time (days)", size=12, fontname="serif")
	ax2.set_ylabel(r"Volume loss (%)", size=12, fontname="serif")

def main():
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3))
	fig.subplots_adjust(top=0.970, bottom=0.16, left=0.114, right=0.98, hspace=0.21, wspace=0.340)

	output_folder = os.path.join("output", "case_0", "vtk")
	plot_res(ax1, ax2, output_folder, "steelblue", "Final shape", True)

	# apply_dark_theme(fig, [ax1, ax2], transparent=False)
	# apply_special_theme(fig, [ax1, ax2], transparent=True)
	apply_white_theme(fig, [ax1, ax2], transparent=True)

	plt.show()

if __name__ == '__main__':
	main()