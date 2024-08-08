import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
from ResultsHandler import convert_vtk_to_pandas
import matplotlib.pyplot as plt
import numpy as np
import meshio

minute = 60
hour = 60*minute
day = 24*hour
MPa = 1e6

def trapezoidal_volume(x, y):
	""" This function calculates the volume of a solid of revolution (around y=0 axis) based on the trapezoidal rule. """
	volume = 0.0
	area = 0.0
	n = len(x)
	for i in range(1, n):
		R = 0.5*(y[i] + y[i-1])
		A = np.pi*R**2
		d = x[i] - x[i-1]
		area += R*d
		volume += A*d
	return volume

def reorder_data(df_coord, u, v, w, wall_ind):
	# Initial cavern shape
	x0 = df_coord.iloc[wall_ind]["x"]
	y0 = df_coord.iloc[wall_ind]["y"]
	z0 = df_coord.iloc[wall_ind]["z"]
	# Reorder all coordinates according to coordinate z
	sorted_z0_ind = z0.sort_values().index
	x0 = x0[sorted_z0_ind]
	y0 = y0[sorted_z0_ind]
	z0 = z0[sorted_z0_ind]
	# Reorder all displacements according to coordinate z
	u = u.iloc[wall_ind].loc[sorted_z0_ind]
	v = v.iloc[wall_ind].loc[sorted_z0_ind]
	w = w.iloc[wall_ind].loc[sorted_z0_ind]
	return x0, y0, z0, u, v, w

# Define folders
results_folder = os.path.join("output", "case_0", "operation", "vtk")
mesh_folder = os.path.join("..", "..", "grids", "cavern_irregular")

# Read displacement results
pvd_path = os.path.join(results_folder, "displacement")
pvd_file = "displacement.pvd"
df_coord, u, v, w = convert_vtk_to_pandas(pvd_path, pvd_file)

# Get indices of wall profile
mesh = meshio.read(os.path.join(mesh_folder, "geom.msh"))
wall_ind = np.unique(mesh.cells["line"].flatten())

# Get reordered data over cavern wall
x0, y0, z0, u, v, w = reorder_data(df_coord, u, v, w, wall_ind)

# Get times
times = u.columns.values
t_final = times[-1]
t_initial = times[1]

# Compute cavern volumes over time
vol_0 = trapezoidal_volume(z0.values, x0.values)
volumes = []
for t in times[1:]:
	z = z0.values + w[t].values
	x = x0.values + u[t].values
	vol = trapezoidal_volume(z, x)
	volumes.append(100*abs(vol_0 - vol)/vol_0)

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3))
fig.subplots_adjust(top=0.985, bottom=0.145, left=0.070, right=0.990, hspace=0.35, wspace=0.260)

# Plot cavern shape
expansion_factor = 50
xi = x0 + 0*expansion_factor*u[t_initial]
yi = y0 + 0*expansion_factor*v[t_initial]
zi = z0 + 0*expansion_factor*w[t_initial]
xf = x0 + expansion_factor*u[t_final]
yf = y0 + expansion_factor*v[t_final]
zf = z0 + expansion_factor*w[t_final]
ax1.plot(xi, zi, "-", color="black", linewidth=2.0, label="Initial shape")
ax1.plot(-xi, zi, "-", color="black", linewidth=2.0)
ax1.plot(xf, zf, "-", color="#377eb8", linewidth=2.0, label=f"Final shape")
ax1.plot(-xf, zf, "-", color="#377eb8", linewidth=2.0)
ax1.set_xlabel("x (m)", size=12, fontname="serif")
ax1.set_ylabel("z (m)", size=12, fontname="serif")
ax1.legend(loc=1, shadow=True, fancybox=True)
ax1.axis("equal")
ax1.grid(True, color='0.92')
ax1.set_facecolor("0.85")

# Plot cavern volumetric closure
ax2.plot(times[2:]/hour, volumes[1:] - 0*volumes[1], ".-", color="#377eb8", linewidth="2.0")
ax2.set_xlabel("Time (h)", size=12, fontname="serif")
ax2.set_ylabel("Cavern closure (%)", size=12, fontname="serif")
ax2.grid(True, color='0.92')
ax2.set_facecolor("0.85")

plt.show()
