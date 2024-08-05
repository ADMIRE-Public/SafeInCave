import matplotlib.pyplot as plt
import numpy as np
import json

hour = 60*60
MPa = 1e6

# Read input file
with open("input_file.json", "r") as j_file:
	input_file = json.load(j_file)

# Extract time and stresses
time_list = np.array(input_file["time_settings"]["time_list"])/hour
over_burden = np.array(input_file["boundary_conditions"]["Top"]["values"])/MPa
side_burden = np.array(input_file["boundary_conditions"]["East"]["values"])/MPa
gas_pressure = np.array(input_file["boundary_conditions"]["Cavern"]["values"])/MPa

# Plot loading schedule
fig, ax = plt.subplots(1, 1, figsize=(7, 3))
fig.subplots_adjust(top=0.970, bottom=0.155, left=0.080, right=0.980, hspace=0.35, wspace=0.225)

ax.plot(time_list, over_burden, ".-", color="#377eb8", label="Overburden")
ax.plot(time_list, side_burden, ".-", color="#ff7f00", label="Sideburden")
ax.plot(time_list, gas_pressure, ".-", color="#4daf4a", label="Sideburden")
ax.set_xlabel("Time (h)", size=12, fontname="serif")
ax.set_ylabel("Stress (MPa)", size=12, fontname="serif")
ax.legend(loc=0, shadow=True, fancybox=True)
ax.grid(True, color="0.92")
ax.set_facecolor("0.85")
fig.patch.set_alpha(0.0)
ax.set_ylim(0, 13)

plt.show()