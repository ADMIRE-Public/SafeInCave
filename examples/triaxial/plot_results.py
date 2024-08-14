import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from ResultsHandler import read_vector_from_points

hour = 60*60
MPa = 1e6

# Read displacement results
pvd_path = os.path.join("output", "case_0", "operation", "vtk", "displacement")
pvd_file = "displacement.pvd"
df, u, v, w = read_vector_from_points(pvd_path, pvd_file)

point_P = df[(df["z"]==1) & (df["x"]==1) & (df["y"]==1)].index[0]

Lz = df["z"].values.max()
Lx = df["x"].values.max()

u_P = u.iloc[point_P].values[1:]
w_P = w.iloc[point_P].values[1:]

eps_h = -u_P / Lx
eps_v = -w_P / Lz

t = w.iloc[point_P].index.values[1:]/hour

# Read input file
with open("input_file.json", "r") as j_file:
	input_file = json.load(j_file)

# Extract vertical stresses
time_list = np.array(input_file["time_settings"]["time_list"])/hour
sigma_v_list = np.array(input_file["boundary_conditions"]["TOP"]["values"])/MPa
sigma_v = np.interp(t, time_list, sigma_v_list)

# Plot pressure schedule
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3))
fig.subplots_adjust(top=0.970, bottom=0.160, left=0.085, right=0.990, hspace=0.35, wspace=0.225)

ax1.plot(t, eps_v*100, ".-", color="#377eb8", label=r"$\varepsilon_v$")
ax1.plot(t, eps_h*100, ".-", color="#ff7f00", label=r"$\varepsilon_h$")
ax1.set_xlabel("Time (h)", size=12, fontname="serif")
ax1.set_ylabel("Total strain (%)", size=12, fontname="serif")
ax1.legend(loc=0, shadow=True, fancybox=True)
ax1.grid(True, color="0.92")
ax1.set_facecolor("0.85")
fig.patch.set_alpha(0.0)

ax2.plot(eps_v*100, sigma_v, ".-", color="#377eb8")
ax2.set_xlabel("Vertical strain (%)", size=12, fontname="serif")
ax2.set_ylabel("Vertical stress (MPa)", size=12, fontname="serif")
ax2.grid(True, color="0.92")
ax2.set_facecolor("0.85")
fig.patch.set_alpha(0.0)

plt.show()
