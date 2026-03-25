import numpy as np
import matplotlib.pyplot as plt

#%% PARAMETERS

# ambient 
g = 9.81            # Gravity [m/s^2]
rhow = 1025         # Water density [kg/m^3]

# floater and pontoons parameters
draft = np.array([0.10,0.20,0.30,0.40])    # length under water [m]
Diameter = np.linspace (0.1,0.35,6)        # floater diameter [m]
D_pon = 0.1                                # pontoon diameter [m]

# tower parameters
zhub = 1.2                                 # hub hight [m]

# maximum pitch parameters
theta_deg = np.deg2rad(10)                 # pitch angle [deg]
T = 85                                     # thust single rotor[N]

#--------------------------------------------
# Calculate the length of the pontoons
#--------------------------------------------
M = 2*T*zhub #overturning moment [Nm]
c55 = M/theta #restoring [Nm/rad]
L_p = np.sqrt((4*c55/(3*rhow*g*np.pi*(Diameter/2)**2))-(Diameter/2)**2/3) #lenght pontoon [m]
param1 = c55/(rhow*g)
param2 = 3*np.pi*((Diameter/2)**4)/4
param3 = np.pi*((Diameter/2)**2)/4
L_p = np.sqrt((param1 - param2)/param3)

plt.figure()
plt.plot(Diameter,L_p, marker='o', label=f'Pontoons diameter = {D_pon} m')
plt.xlabel('Diameter (m)')
plt.ylabel('Pontoons lenght (m)')
plt.title(f'Pontoons length to obtain pitch of {theta_deg} deg')
plt.grid()

#---------------------------------------------
# Create a grid to avoid for loop
#---------------------------------------------
D_grid, d_grid = np.meshgrid(Diameter, draft)
Lp_grid, d_grid = np.meshgrid(L_p, draft)

#--------------------------------------------
# Calculate the total mass
#--------------------------------------------
mtot = rhow * d_grid * np.pi * 3 *(D_grid / 2)**2 #total mass excluding the pontoon
mtot_pontoon = mtot + rhow *( 3*(Lp_grid-D_grid) * np.pi* (D_pon/2)**2 )    #total mass including the pontoon
                                                                            #(Lp_grid-D_grid) because Lp is distance
                                                                            # between the centers two floaters  
#--------------------------------------------
# Plotting
#--------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

for i, dr in enumerate(draft):
    ax1.plot(Diameter, mtot[i, :], marker='o',linestyle= '--', label=f'Draft={draft[i]} m ')
    ax2.plot(Diameter, mtot_pontoon[i, :], marker='o', label=f'Draft={draft[i]}')

ax1.set_ylabel('$m_{tot}$ (kg)')
ax1.set_title('Total Mass with NO pontoons')
ax1.legend(loc='upper left', bbox_to_anchor=(1, 1))
ax1.grid(True, linestyle='--', alpha=0.5)

ax2.set_xlabel('Floater diameter (m)')
ax2.set_ylabel('$m_{tot}$ (kg)')
ax2.set_title(f'Total Mass with pontoons of D = {D_pon}m')
ax2.legend(loc='upper left', bbox_to_anchor=(1, 1))
ax2.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()
