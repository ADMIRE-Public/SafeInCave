import safeincave as sf
import safeincave.Utils as ut
from safeincave.Utils import GPa, MPa, day, hour, create_field_elems, create_field_nodes
import safeincave.HeatBC as heatBC
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def apply_grey_theme(fig, axes, transparent=True, grid_color="0.92", back_color='0.85'):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
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

def get_geometry_parameters(path_to_grid):
	f = open(os.path.join(path_to_grid, "geom.geo"), "r")
	data = f.readlines()
	ovb_thickness = float(data[10][len("ovb_thickness = "):-2])
	salt_thickness = float(data[11][len("salt_thickness = "):-2])
	hanging_wall = float(data[12][len("hanging_wall = "):-2])
	return ovb_thickness, salt_thickness, hanging_wall

def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cavern_overburden_coarse")
	grid = sf.GridHandlerGMSH("geom", grid_path)

	gas_density = 0.082
	salt_density = 2200
	ovb_density = 2800
	ovb_thickness, salt_thickness, hanging_wall = get_geometry_parameters(grid_path)
	cavern_roof = ovb_thickness + hanging_wall
	g = 9.81
	p_roof = 0 + salt_density*g*hanging_wall + ovb_density*g*ovb_thickness

	p_values = 3*[0.8*p_roof, 0.8*p_roof, 0.2*p_roof, 0.2*p_roof] + [0.8*p_roof]
	t_values = [20*day*i for i in range(13)]

	p_values = np.array(p_values)/MPa
	t_values = np.array(t_values)/day

	fig, ax = plt.subplots(1, 1, figsize=(5, 3.5))
	fig.subplots_adjust(top=0.92, bottom=0.155, left=0.138, right=0.980, hspace=0.35, wspace=0.295)

	ax.plot([-10, 0.0], [0.8*p_roof/MPa, 0.8*p_roof/MPa], "-", color="lightcoral", linewidth=2.0, label="Equilibrium stage")
	ax.plot(t_values, p_values, "-", color="steelblue", linewidth=2.0, label="Operation stage")
	ax.set_xlabel("Time (days)", fontname="serif", size=12)
	ax.set_ylabel("Gas pressure (MPa)", fontname="serif", size=12)
	ax.set_ylim(2.5, 22)
	ax.legend(loc=2, fancybox=True, shadow=True, ncol=2)
	ax.grid(True)

	apply_grey_theme(fig, [ax])

	plt.show()

if __name__ == '__main__':
	main()