"""
This code determines the matrices of the floating wind turbine and the relative
eigen frequencies.
"""

#%% IMPORT PACKAGES
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.linalg import eig

#%% AMBIENT AND FLOATER PARAMETERS
# The hub height of the turbines (base of tower + height of floaters above water
# should be limited to 1.2 m)!!!!!!!

# ambient
g = 9.81            # Gravity [m/s^2]
h = 2               # Water depth [m]
rhow = 1025         # Water density [kg/m^3]
rhoa = 1.29         # Air density [kg/m^3]
rho_pvc = 1400      # PVC density [kg/m^3]

# tower and rotor parameters
z_hub = 1.2        # tower height above SWL [kg]
D_tow = 38.2e-2    # tower outer diameter [m]
D_tow_in = 37.2e-2 # tower inner diameter [m]

D_rot = 1.2        # rotor diameter [m]
m_rot = 8       # rotor + nacelle [kg]


# floater parameters

m_bf = 5           # Back Floater mass [kg]
m_ff = 5           # Front Floater mass [kg]
h_bf = 0.5         # Back floater hight, including ballast [m]
h_ff = 0.5         # Front floater hight, including ballast [m] 
D_b = 0.2         # Back floater Diameter [m]
D_b_in = 0.19     # Back floater inner diameter [m]
D_f = 0.2          # Front floater Diameter [m]
D_f_in = 0.19      # Front floater inner diameter [m]
A_ff_solid = (np.pi*(D_f)**2)/4
A_bf_solid = (np.pi*(D_b)**2)/4



m_bb = 10          # Back Ballast mass [kg]

h_bb = 0.2         # Back Ballast height [m]
h_fb = 0.2         # Front Ballast height [m] 

edge = 1           # distance between columns: equi lateral triangle

# mooring system
K_moor = 12        # mooring stiffness [N/m]
z_moor = 0.3 

draft = 0.3 # Draft [m]


# Mass of tower
h_out = h_ff - draft
h_left = z_hub - h_out # Length of the tower [m]\
m_tow = rho_pvc * np.pi*(((D_tow/2)**2) - ((D_tow_in/2)**2))*h_left

# MASS CALCULATION

L_pontoon = 2.15 # Length of the pontoon [m]
D_pontoon = 0.1 # Diameter of the pontton [m]
D_pontoon_in = 0.09 # Inner diameter of pontoon [m]


m_tot = rhow*((A_ff_solid*draft) + (2*A_bf_solid*draft) + (3*(L_pontoon*np.pi*(D_pontoon/2)**2)))

m_pontoon = rho_pvc*np.pi*((D_pontoon/2)**2-(D_pontoon_in/2)**2)*L_pontoon

m_back = (m_tot - (3*m_pontoon))/(2 + ((D_f**2)/(D_b**2))) # For a single back floater

m_bb = (m_back - (m_rot + m_tow + m_bf))

m_front = m_back*((np.pi*((D_f/2)**2))/(2*(np.pi*((D_b/2)**2)))) # Front side mass [kg] (considering only the OD)

m_fb = m_front - m_ff  # Calculation of front ballast to maintain the same draft                       



A_back = np.pi * D_b**2 / 4                   # Sectional area back cylinders
                                              # individually [m^2]

A_front = np.pi * D_f**2 / 4                  # Sectional area front cylinder [m^2]

A_total = 2 * A_back + A_front                # total sectional area [m^2]



# printing
if draft >= (h_bf + h_bb):
    sys.exit(f'\nWE ARE SINKING: draft ({draft} m) is bigger than height of floater')
else:
    if h_bb >= draft:
        print('\nThe draft is higher than the height of the ballast')
    else:
        print(f'\nThe draft is {draft} m ')
        print(f'\nYou have {h_out} m of the floater above water')
        print(f'The maximum length of the tower is {h_left} m')

#%% COORDINATE SYSTEM
#   Determine positions of columns based on center of the back floaters.
#   The x axis is pointing from the front floater to the ones in the back,
#   in the direction of the wind. The y axis follows the right hand rule.    

# The back floaters are N° 1 and N° 2, while the front is N°3
y_1 = - edge / 2
y_2 = + edge / 2
y_3 = 0
x_1 = 0
x_2 = 0
x_3 = - np.sqrt(3)*edge/2

#%% CM COORDINATES    

z_CM_tow = (z_hub - h_out) / 2        # z center of mass tower

z_CM_bf = h_out - h_bf / 2            # z center of back floater
z_CM_ff = h_out - h_bf / 2            # z center of front floater

z_CM_bb = (h_bb / 2) - draft          # z center of back ballast
z_CM_fb = (h_fb / 2) - draft          # z center of front ballast

# z position of center of mass
z_CM_tot = (2 * z_CM_tow * m_tow + 2 * z_hub * m_rot + 2 * z_CM_bf * m_bf +\
            2 * z_CM_bb * m_bb + z_CM_ff * m_ff + z_CM_fb * m_fb) / m_tot

print(f'The center of mass is at {z_CM_tot} m')

# y position center of mass
x_CM_tot = (m_back * x_1 + m_back * x_2 + m_front *x_3) / m_tot

print(f'The center of mass around y is at {x_CM_tot} m')

# x position center of mass
y_CM_tot = (m_back * y_1 + m_back * y_2 + m_front *y_3) / m_tot

print(f'The center of mass around x is at {y_CM_tot} m')

#%% CB COORDINATES

# z position of center of buoyancy
z_CB_tot = - draft / 2

print(f'The center of bouyancy is at {z_CB_tot} m')

Vol_1_sub = draft * A_back
Vol_2_sub = draft * A_back
Vol_3_sub = draft * A_front

# x position of center of mass
x_CB_tot = (Vol_1_sub * x_1 + Vol_2_sub * x_2 + Vol_3_sub * x_3) /\
                    (Vol_1_sub + Vol_2_sub + Vol_3_sub)
                    
print(f'The center of bouyancy around y is at {x_CB_tot} m')
      
# y position of center of mass              
y_CB_tot = (Vol_1_sub * y_1 + Vol_2_sub * y_2 + Vol_3_sub * y_3) /\
                    (Vol_1_sub + Vol_2_sub + Vol_3_sub)

print(f'The center of bouyancy around x is at {y_CB_tot} m')

#%% INERTIA CALCULATION FOR MASS  AND ADDED MASS MATRICES

# turbine: considered as point mass
I0_XX_turb = 0 
I0_YY_turb = m_rot * y_1 **2  # first approximation (actually is an incline tower)
I0_ZZ_turb = m_rot * z_hub**2
I0_XY_turb = 0
I0_XZ_turb = 0
I0_YZ_turb_1 = m_rot * y_1 * z_hub
I0_YZ_turb_2 = m_rot * y_2 * z_hub

# tower: it's a thick-walled cylinder
I0_XX_tower = m_tow * (3 * ((D_tow / 2)**2 + (D_tow_in / 2)**2)+ h_left**2) / 12
I0_YY_tower = m_tow * (3 * ((D_tow / 2)**2 + (D_tow_in / 2)**2)+ h_left**2) / 12 +\
              m_tow * y_1**2
I0_ZZ_tower = 0.5 *m_tow * ((D_tow / 2)**2 + (D_tow_in / 2)**2) + m_tow * z_CM_tow**2 
I0_XY_tower = 0
I0_XZ_tower = 0
I0_YZ_tower_1 = m_rot * y_1 * z_CM_tow
I0_YZ_tower_2 = m_rot * y_2 * z_CM_tow

# back floater: it's a thick-walled cylinder
I0_XX_backfloat =  m_bf * (3 * ((D_b / 2)**2 + (D_b_in / 2)**2)+ h_bf**2) / 12
I0_YY_backfloat =  m_bf * (3 * ((D_b / 2)**2 + (D_b_in / 2)**2)+ h_bf**2) / 12 +\
                   m_bf * y_1**2
I0_ZZ_backfloat = 0.5 *m_bf * ((D_b / 2)**2 + (D_b_in / 2)**2) + m_bf * z_CM_bf**2 
I0_XY_backfloat = 0
I0_XZ_backfloat = 0
I0_YZ_backfloat_1 = m_bf * y_1 * z_CM_bf
I0_YZ_backfloat_2 = m_bf * y_2 * z_CM_bf

# front floater: it's a thick-walled cylinder
I0_XX_frontfloat =  m_ff * (3 * ((D_f / 2)**2 + (D_f_in / 2)**2)+ h_ff**2) / 12 +\
                   m_ff * x_3 **2
I0_YY_frontfloat =  m_ff * (3 * ((D_f / 2)**2 + (D_f_in / 2)**2)+ h_ff**2) / 12
I0_ZZ_frontfloat = 0.5 *m_ff * ((D_f / 2)**2 + (D_f_in / 2)**2) + m_ff * z_CM_ff**2 
I0_XY_frontfloat = 0
I0_XZ_frontfloat = m_ff * x_3 * z_CM_ff
I0_YZ_frontfloat = 0

# back ballast: it's a solid cylinder
I0_XX_backball =  m_bb * (3 * (D_b / 2)**2+ h_bb**2) / 12
I0_YY_backball =  m_bb * (3 * (D_b / 2)**2+ h_bb**2) / 12 + m_bb * y_1**2
I0_ZZ_backball = 0.5 *m_bb * (D_b / 2)**2 + m_bb * z_CM_bb**2
I0_XY_backball = 0
I0_XZ_backball = 0
I0_YZ_backball_1 = m_bb * y_1 * z_CM_bb
I0_YZ_backball_2 = m_bb * y_2 * z_CM_bb

# front ballast: it's a solid cylinder
I0_XX_frontball =  m_fb * (3 * (D_f / 2)**2+ h_fb**2) / 12 + m_fb * x_3**2
I0_YY_frontball =  m_fb * (3 * (D_f / 2)**2+ h_fb**2) / 12
I0_ZZ_frontball = 0.5 *m_fb * (D_f / 2)**2 + m_fb * z_CM_fb**2
I0_XY_frontball = 0
I0_XZ_frontball = m_fb * x_3 * z_CM_fb
I0_YZ_frontball = 0

# sum of inertia components
I0_XX = 2*(I0_XX_turb + I0_XX_backball + I0_XX_backfloat + I0_XX_tower) +\
        I0_XX_frontfloat + I0_XX_frontball
I0_YY = 2*(I0_YY_turb + I0_YY_backball + I0_YY_backfloat + I0_YY_tower) +\
        I0_YY_frontfloat + I0_YY_frontball
I0_ZZ = 2*(I0_ZZ_turb + I0_ZZ_backball + I0_ZZ_backfloat + I0_ZZ_tower) +\
        I0_ZZ_frontfloat + I0_ZZ_frontball
I0_XY = 2*(I0_XY_turb + I0_XY_backball + I0_XY_backfloat + I0_XY_tower) +\
        I0_XY_frontfloat + I0_XY_frontball
I0_XZ = 2*(I0_XZ_turb + I0_XZ_backball + I0_XZ_backfloat + I0_XZ_tower) +\
        I0_XZ_frontfloat + I0_XZ_frontball
I0_YZ = I0_YZ_turb_1 + I0_YZ_backball_1 + I0_YZ_backfloat_1 + I0_YZ_tower_1 +\
        I0_YZ_turb_2 + I0_YZ_backball_2 + I0_YZ_backfloat_2 + I0_YZ_tower_2+\
        I0_YZ_frontfloat + I0_YZ_frontball
#%% MASS MATRIX

# surge
M11 = m_tot                   # surge-surge
M12 = 0                       # surge-sway
M13 = 0                       # surge-heave
M14 = 0                       # surge-roll
M15 = m_tot * z_CM_tot        #surge-pitch
M16 = - m_tot * y_CM_tot      # surge-yaw

# sway
M21 = M12                     # sway-surge
M22 = m_tot                   # sway-sway
M23 = 0                       # sway-heave
M24 = -m_tot * z_CM_tot       # sway-roll
M25 = 0                       # sway-pitch
M26 = m_tot * x_CM_tot        # sway-yaw

# heave
M31 = M13                     # heave-surge
M32 = M23                     # heave-sway
M33 = m_tot                   # heave-heave
M34 = m_tot * y_CM_tot        # heave-roll
M35 = - m_tot * x_CM_tot      # heave-pitch
M36 = 0                       # heave-yaw

# roll
M41 = M14                     # roll-surge
M42 = M24                     # roll-sway
M43 = M34                     # roll-heave
M44 = I0_YY + I0_ZZ           # roll-roll 
M45 = - I0_XY                 # roll-pitch
M46 = - I0_XZ                 # roll-yaw

# pitch
M51 = M15                    # pitch-surge
M52 = M25                    # pitch-sway
M53 = M35                    # pitch-heave
M54 = M45                    # pitch-roll
M55 = I0_XX + I0_ZZ          # pitch-pitch
M56 = - I0_YZ                # pitch-yaw

# yaw
M61 = M16                    # yaw-surge
M62 = M26                    # yaw-sway
M63 = M36                    # yaw-heave
M64 = M46                    # yaw-roll
M65 = M56                    # yaw-pitch
M66 = I0_XX + I0_YY          # yaw-yaw


M = np.array([  [M11, M12, M13, M14, M15, M16],
                [M21, M22, M23, M24, M25, M26],
                [M31, M32, M33, M34, M35, M36],
                [M41, M42, M43, M44, M45, M46],
                [M51, M52, M53, M54, M55, M56],
                [M61, M62, M63, M64, M65, M66]])


#%% ADDED MASS MATRIX

A = np.zeros_like(M)

#%% DAMPING MATRIX

B = np.zeros_like(M)

#%% MOMENTS OF WATERPLANE AREA

# Linear moments
I_x = 2 * A_back * x_1 +  A_front * x_3

I_y = A_back *y_1 + A_back * y_2 + A_front * y_3

# Second-order moments
I_xx = 2 * np.pi * D_b**4 / 64 + np.pi * D_f**4 / 64 +\
                    2 * A_back * x_1**2 +\
                        A_front * x_3**2
                        
I_yy = 2 * np.pi * D_b**4 / 64 + np.pi * D_f**4 / 64 +\
                    A_back * y_1**2 + A_back * y_2**2 +\
                    A_front * y_3**2
                    
I_xy = 2 * np.pi * D_b**4 / 64 + np.pi * D_f**4 / 64 +\
                    A_back * y_1 * x_1 + A_back * y_2 * x_2 +\
                    A_front * y_3 * x_3

#%% HYDROSTATIC STIFFNESS MATRIX

# surge
C_11 = 0       # surge-surge
C_12 = 0       # surge-sway
C_13 = 0       # surge-heave
C_14 = 0       # surge-roll
C_15 = 0       # surge-pitch
C_16 = 0       # surge-yaw

# sway
C_21 = 0       # sway-surge
C_22 = 0       # sway-sway
C_23 = 0       # sway-heave
C_24 = 0       # sway-roll
C_25 = 0       # sway-pitch
C_26 = 0       # sway-yaw


# heave
C_31 = 0       # heave-surge
C_32 = 0       # heave-sway
C_33 = 2 * rhow * g * A_back + rhow * g * A_front   # heave-heave
C_34 = rhow * g * I_y       # heave-roll
C_35 = - rhow * g * I_x     # heave-pitch
C_36 = 0       # heave-yaw

# roll
C_41 = 0       # roll-surge
C_42 = 0       # roll-sway
C_43 = C_34    # roll-heave
C_44 =  rhow * g * I_yy + m_tot * g * (z_CB_tot - z_CM_tot) # roll-roll
C_45 = - rhow * g * I_xy   # roll-pitch
C_46 = 0       # roll-yaw

# pitch
C_51 = 0       # pitch-surge
C_52 = 0       # pitch-sway
C_53 = C_35    # pitch-heave
C_54 = C_45    # pitch-roll
C_55 = m_tot * g * (z_CB_tot - z_CM_tot) + rhow * g * I_xx   # pitch-pitch
C_56 = 0       # pitch-yaw

# yaw
C_61 = 0       # yaw-surge
C_62 = 0       # yaw-sway
C_63 = 0       # yaw-heave
C_64 = - m_tot * g * (x_CB_tot - x_CM_tot)   # yaw-roll 
C_65 = - m_tot * g * (y_CB_tot - y_CM_tot)   # yaw-pitch
C_66 = 0       # yaw-yaw

C_hydro = np.array([[C_11, C_21, C_31, C_41, C_51, C_61],
                    [C_12, C_22, C_32, C_42, C_52, C_62],
                    [C_13, C_23, C_33, C_43, C_53, C_63],
                    [C_14, C_24, C_34, C_44, C_54, C_64],
                    [C_15, C_25, C_35, C_45, C_55, C_65],
                    [C_16, C_26, C_36, C_46, C_56, C_66]])

#%% MOORING MATRIX

C_mooring = np.zeros_like(C_hydro)
#%% TOTAL STIFFNESS MATRIX

C_tot = C_hydro + C_mooring

#%% NATURAL FREQUENCIES EVALUATION


CoMA = np.linalg.solve(M + A, C_tot)
eigVal, eigVec = np.linalg.eig(CoMA)

# Natural frequencies
omega_nat = np.sqrt(eigVal)
f_nat = omega_nat / (2 * np.pi)

# Display natural frequencies
print(f'Surge period: {f_nat[0]:.2f} [Hz]')
print(f'Sway period: {f_nat[1]:.2f} [Hz]')
print(f'Heave period: {f_nat[2]:.2f} [Hz]')
print(f'Roll period: {f_nat[3]:.2f} [Hz]')
print(f'Pitch period: {f_nat[4]:.2f} [Hz]')
print(f'Yaw period: {f_nat[5]:.2f} [Hz]')





