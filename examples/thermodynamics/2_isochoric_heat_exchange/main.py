import numpy as np
import matplotlib.pyplot as plt
from safeincave.Thermodynamics import CavernThermodynamics

# Define conversion units
minute = 60
hour = 60*minute
day = 24*hour
cm = 1e-2
MPa = 1e6
ton = 1000
Mton = ton*1e6

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



def heat(Q_min, Q_max):
    time = np.linspace(0, 100, 100)
    Q_list = Q_min + (Q_max - Q_min)*time/time[-1]
    return time, Q_list



def calculate_thermodynamics(fluid_name, Q_list, P_init=10*MPa, T_init=300):
    model = CavernThermodynamics(fluid_name)
    P_atm = 101325.0
    rho_init, _, _ = model.rho_u_h(P_init + P_atm, T_init)
    T, P, rho = [T_init], [P_init], [rho_init]
    for i in range(1, len(Q_list)):
        P1, T1, rho1 = model.solve( dm = 0.0,
                                    Q_in = Q_list[i],
                                    T_in = 0.0,
                                    P0 = P[i-1] + P_atm,
                                    T0 = T[i-1],
                                    V0 = 1,
                                    V1 = 1)
        P.append(P1 - P_atm)
        T.append(T1)
        rho.append(rho1)
    return np.array(P), np.array(T), np.array(rho)
	
	

def main():
    fig, axis = plt.subplots(2, 4, figsize=(16, 8))
    fig.subplots_adjust(top=0.935, bottom=0.155, left=0.05, right=0.980, hspace=0.35, wspace=0.38)
    ax1, ax2, ax3, ax4 = axis[0,:]
    ax5, ax6, ax7, ax8 = axis[1,:]
	
    time, Q_heat_up = heat(7e5, 7e5)
    P, T, rho = calculate_thermodynamics("Water", Q_heat_up)
	
    ax1.plot(time, Q_heat_up/1e3, "-", linewidth=2.0, color="black")
    ax1.set_xlabel("Time (s)", fontname="serif", size=12)
    ax1.set_ylabel("Heat (kJ)", fontname="serif", size=12)
	
    ax2.plot(time, P/MPa, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax2.set_xlabel("Time", fontname="serif", size=12)
    ax2.set_ylabel("Pressure (MPa)", fontname="serif", size=12)
	
    ax3.plot(time, T-273, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax3.set_xlabel("Time", fontname="serif", size=12)
    ax3.set_ylabel("Temperature (°C)", fontname="serif", size=12)
	
    ax4.plot(time, rho, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax4.set_xlabel("Time", fontname="serif", size=12)
    ax4.set_ylabel("Density (kg/m3)", fontname="serif", size=12)
    ax4.set_ylim(1000.4, 1001.4)
	

    time, Q_cool_down = heat(-7e5, -7e5)
    P, T, rho = calculate_thermodynamics("Water", Q_cool_down)
	
    ax5.plot(time, Q_cool_down/1e3, "-", linewidth=2.0, color="black")
    ax5.set_xlabel("Time (s)", fontname="serif", size=12)
    ax5.set_ylabel("Heat (kJ)", fontname="serif", size=12)
	
    ax6.plot(time, P/MPa, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax6.set_xlabel("Time", fontname="serif", size=12)
    ax6.set_ylabel("Pressure (MPa)", fontname="serif", size=12)
	
    ax7.plot(time, T-273, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax7.set_xlabel("Time", fontname="serif", size=12)
    ax7.set_ylabel("Temperature (°C)", fontname="serif", size=12)
	
    ax8.plot(time, rho, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax8.set_xlabel("Time", fontname="serif", size=12)
    ax8.set_ylabel("Density (kg/m3)", fontname="serif", size=12)
    ax8.set_ylim(1000.4, 1001.4)
	
    apply_grey_theme(fig, axis.flatten())
	
    plt.show()



if __name__ == "__main__":
	main()

