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
sigma_v = np.array(input_file["boundary_conditions"]["TOP"]["values"])/MPa
sigma_h = np.array(input_file["boundary_conditions"]["EAST"]["values"])/MPa

# Plot loading schedule
fig, ax = plt.subplots(1, 1, figsize=(7, 3.5))
fig.subplots_adjust(top=0.970, bottom=0.135, left=0.080, right=0.980, hspace=0.35, wspace=0.225)

ax.plot(time_list, sigma_v, ".-", color="#377eb8", label=r"$\sigma_v$")
ax.plot(time_list, sigma_h, ".-", color="#ff7f00", label=r"$\sigma_h$")
ax.set_xlabel("Time (h)", size=12, fontname="serif")
ax.set_ylabel("Stress (MPa)", size=12, fontname="serif")
ax.legend(loc=0, shadow=True, fancybox=True)
ax.grid(True, color="0.92")
ax.set_facecolor("0.85")
fig.patch.set_alpha(0.0)
ax.set_ylim(0, 13)

plt.show()