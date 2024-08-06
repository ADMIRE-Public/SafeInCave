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

salt_density = input_file["boundary_conditions"]["East"]["density"]
h2_density = input_file["boundary_conditions"]["Cavern"]["density"]
g = input_file["body_force"]["gravity"]
H = 660
H_roof = 430
H_sump = 205.2

# Plot loading schedule
fig = plt.figure()
fig.subplots_adjust(top=0.970, bottom=0.155, left=0.080, right=0.980, hspace=0.35, wspace=0.225)
ax = fig.add_subplot(projection='3d')
ax.set_xlabel("Time (h)")
ax.set_ylabel("Pressure (MPa)")
ax.set_zlabel("Height (m)")

# p_top = side_burden[0]
# p_bot = p_top - g*salt_density*H/MPa
# t0 = time_list[0]
# x = np.array([t0, t0, t0, t0, t0, t0])
# y = np.array([0, p_top, p_bot, p_top, 0, 0])
# z = np.array([0, 0, 0, H, H, 0])
# ax.plot(x, y, z, ".-")

# p_top = side_burden[-1]
# p_bot = p_top - g*salt_density*H/MPa
# tf = time_list[-1]
# x = np.array([tf, tf, tf, tf, tf, tf])
# y = np.array([0, p_top, p_bot, p_top, 0, 0])
# z = np.array([0, 0, 0, H, H, 0])
# ax.plot(x, y, z, ".-")

xb, yb, zb = [], [], []
for i in [0, -1]:
	p_top = side_burden[i]
	tf = time_list[i]
	p_bot = p_top - g*salt_density*H/MPa
	x = np.array([tf, tf, tf, tf, tf])
	y = np.array([0, p_bot, p_top, 0, 0])
	z = np.array([0, 0, H, H, 0])
	xb.append(x)
	yb.append(y)
	zb.append(z)
	ax.plot(x, y, z, ".-", color="steelblue")
	if i != 0:
		ax.plot([x[0], xb[i-1][0]], [y[0], yb[i-1][0]], [z[0], zb[i-1][0]], "--", color="0.5", linewidth=1.0)
		ax.plot([x[1], xb[i-1][1]], [y[1], yb[i-1][1]], [z[1], zb[i-1][1]], "--", color="0.5", linewidth=1.0)
		ax.plot([x[2], xb[i-1][2]], [y[2], yb[i-1][2]], [z[2], zb[i-1][2]], "--", color="0.5", linewidth=1.0)
		ax.plot([x[3], xb[i-1][3]], [y[3], yb[i-1][3]], [z[3], zb[i-1][3]], "--", color="0.5", linewidth=1.0)
		ax.plot([x[4], xb[i-1][4]], [y[4], yb[i-1][4]], [z[4], zb[i-1][4]], "--", color="0.5", linewidth=1.0)
ax.plot(x, y, z, ".-", color="steelblue", label="Sideburden")

xb, yb, zb = [], [], []
for i in range(len(gas_pressure)):
	p_roof = gas_pressure[i]
	t_mid = time_list[i]
	p_bot = p_roof - g*h2_density*(H_roof - H_sump)/MPa
	x = np.array([t_mid, t_mid, t_mid, t_mid, t_mid, t_mid])
	y = np.array([0, p_roof, p_bot, p_roof, 0, 0])
	z = np.array([H_sump, H_sump, H_sump, H_roof, H_roof, H_sump])
	xb.append(x)
	yb.append(y)
	zb.append(z)
	ax.plot(x, y, z, ".-", color="lightcoral")
	if i > 0:
		ax.plot([x[0], xb[i-1][0]], [y[0], yb[i-1][0]], [z[0], zb[i-1][0]], "--", color="0.5", linewidth=1.0)
		ax.plot([x[1], xb[i-1][1]], [y[1], yb[i-1][1]], [z[1], zb[i-1][1]], "--", color="0.5", linewidth=1.0)
		ax.plot([x[2], xb[i-1][2]], [y[2], yb[i-1][2]], [z[2], zb[i-1][2]], "--", color="0.5", linewidth=1.0)
		ax.plot([x[3], xb[i-1][3]], [y[3], yb[i-1][3]], [z[3], zb[i-1][3]], "--", color="0.5", linewidth=1.0)
		ax.plot([x[4], xb[i-1][4]], [y[4], yb[i-1][4]], [z[4], zb[i-1][4]], "--", color="0.5", linewidth=1.0)
		ax.plot([x[5], xb[i-1][5]], [y[5], yb[i-1][5]], [z[5], zb[i-1][5]], "--", color="0.5", linewidth=1.0)
ax.plot(x, y, z, ".-", color="lightcoral", label="Gas pressure")

ax.legend(loc=0, shadow=True, fancybox=True)



# ax.plot(time_list, over_burden, ".-", color="#377eb8", label="Overburden")
# ax.plot(time_list, side_burden, ".-", color="#ff7f00", label="Sideburden")
# ax.plot(time_list, gas_pressure, ".-", color="#4daf4a", label="Gas pressure")
# ax.set_xlabel("Time (h)", size=12, fontname="serif")
# ax.set_ylabel("Stress (MPa)", size=12, fontname="serif")
# ax.legend(loc=0, shadow=True, fancybox=True)
# ax.grid(True, color="0.92")
# ax.set_facecolor("0.85")
# fig.patch.set_alpha(0.0)
# ax.set_ylim(0, 13)

plt.show()