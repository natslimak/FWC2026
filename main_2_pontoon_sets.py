"""
This code determines the matrices of the floating wind turbine and the relative
eigen frequencies by setting a specific draft.
"""

#%% IMPORT PACKAGES
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.linalg import eig
from scipy.optimize import fsolve

#%% AMBIENT AND TOWER PARAMETERS (EXLUDING MASS)
# The hub height of the turbines (base of tower + height of floaters above water
# should be limited to 1.2 m)!!!!!!!

# ambient
g = 9.81            # Gravity [m/s^2]
h = 2               # Water depth [m]
rhow = 1025         # Water density [kg/m^3]
rhoa = 1.29         # Air density [kg/m^3]
rho_pvc = 1400      # PVC density [kg/m^3]
rho_gf = 2460       # glass fiber density [kg/m^3]
rho_pla = 1240       # PLA density [kg/m^3]


# tower and rotor parameters
z_hub = 1.2        # tower height above SWL [kg]
D_tow = 0.15    # tower outer diameter [m]
D_tow_in = 0.14    # tower inner diameter [m]
D_rot = 1.2        # rotor diameter [m]

# floater and ballast parameters
h_bf = 0.5         # Back floater height, including ballast [m]
h_ff = 0.5         # Front floater hight, including ballast [m] 
h_bb = 0.2         # Back Ballast height [m]
h_fb = 0.2         # Front Ballast height [m] 
D_b = 0.23          # Back floater Diameter [m]
D_b_in = 0.22      # Back floater inner diameter [m]
D_f = 0.23          # Front floater Diameter [m]
D_f_in = 0.22    # Front floater inner diameter [m]

# pontoons parameters
D_pon = 0.049        # Outer iameter of the pontton [m]
D_pon_in = 0.041  # Inner diameter of pontoon [m]

# heave plate
D_hp = 0.3         # Diameter heave plate [m]
h_hp = 0.03       # Height heave plate [m]
 
# draft
draft = 0.3        # draft [m]

# maximum thrust and pitch allowed
T = 85             # max thrust per rotor [N]
max_pitch = 10     # maximum pitch angle [deg]

#%% FLOATER AND PONTOON FRONTAL AREAS
A_back = np.pi * D_b**2 / 4              # Sectional area back cylinder [m^2]

A_front = np.pi * D_f**2 / 4             # Sectional area front cylinder [m^2]

A_pon = np.pi * D_pon**2 / 4             # Sectional area pontoon [m^2]     

A_hp = np.pi * D_hp**2 / 4               # Sectional area heave plate [m^2]

#%% LENGTH OF PONTOONS CALCULATION

# maximum overturning moment
max_M = 2 * T * z_hub                    # [Nm]

# maximum pitch-pitch component of hydrostatic stiffness matrix
c_55 = max_M / np.deg2rad(max_pitch)     # [Nm/rad]

# get length of pontoons to resist maximum overturning moment
def solve_equation(func, initial_guess):
    solution = fsolve(func, initial_guess)
    return solution[0]

def equation(L):
    I_xx_floaters = 3 * np.pi * D_b**4 / 64 + A_front * 3 / 4 * L**2
    # I_xx_pontoons = 3* np.pi * D_pon **4 / 64 + 2 *  A_pon  * 3 / 16 * L**2
    # return c_55 - rhow * g * (I_xx_floaters - I_xx_pontoons)
    return c_55 - rhow * g * (I_xx_floaters)


# The difference between the edge length and the L_pon is that the length of 
# the pontoons does not include the two half diameters of the distance between
# one floater center and another one
edge_length = solve_equation(equation, 2)
L_pon = edge_length - D_b

print(f'The pontoons are {L_pon:2} m long')

#%% MASS CALCULATION

# displaced volume: floater+ballast+pontoons+heave plate
Vol_displaced =  A_front * (draft - h_hp) + (2 * A_back * (draft - h_hp)) + 3 * L_pon * A_pon + 3 * h_hp * A_hp

# total mass
m_tot = Vol_displaced * rhow
print(f'The total mass of the system is {m_tot} kg')

# rotor mass 
m_rot = 6                        # rotor + nacelle [kg]

# tower mass
h_out = h_ff - draft + h_hp      # part of the floater above water [m]
h_left = z_hub - h_out           # Length of the tower from floater to rotor [m]

h_pontoon_above = 0.1  # from the water level to the centre of pontoon above water


m_tow = rho_pvc * np.pi * ((D_tow / 2)**2 - (D_tow_in / 2)**2) * h_left  # [kg]

# mass pontoons
m_pon = rho_pvc * np.pi * ((D_pon / 2)**2 - (D_pon_in / 2)**2) * L_pon   # [kg] 

# mass heave plate
m_hp =  rho_pla * A_hp * h_hp

# mass floater in the back (without ballast)
m_bf = rho_pvc * np.pi * ((D_b / 2)**2 - (D_b_in / 2)**2) * h_bf

# mass floater in the front (without ballast)
m_ff = rho_pvc * np.pi * ((D_f / 2)**2 - (D_f_in / 2)**2) * h_ff

# remaining mass without pontoons and heave plates
m_remaining = m_tot - (6 * m_pon) # Considering 6 pontoons (2 sets of 3 pontoons --> one set below the water and one set above the water)


divisor = ((A_hp * h_hp) + (A_back *(draft - h_hp)))/((A_hp * h_hp) + (A_front * (draft - h_hp)))

m_front = m_remaining / (1 + (2/divisor)) # mass in the front floater (including ballast)

m_back = (m_remaining - m_front)/2 # For one floater component

m_bb = m_back - m_tow - m_rot - m_bf - m_hp # For one floater component
print (f'The mass of the back ballast is {m_bb} kg')

m_fb = m_front - m_ff - m_hp # For one floater component
print (f'The mass of the front ballast is {m_fb} kg')


# # mass back and front (back is floater+ballast+rotor+tower, front is floater+ballast)
# m_back = m_remaining / (1 + A_front /( 2 * A_back))
# m_front = m_remaining - m_back

# # mass ballast in the back
# m_bb = (m_back - 2 * m_rot - 2* m_tow - 2* m_bf) / 2

# # mass ballast in the front
# m_fb = m_front - m_ff

# total mass in the back in each column
m_tot_back = m_rot + m_tow +m_bf + m_bb + m_hp

# total mass in the front
m_tot_front = m_front

#%% COORDINATE SYSTEM
#   Determine positions of columns based on center of the back floaters.
#   The x axis is pointing from the front floater to the ones in the back,
#   in the direction of the wind. The y axis follows the right hand rule.    

# FLOATERS
# The back floaters are N° 1 and N° 2, while the front is N°3
y_1 = - edge_length / 2
y_2 = + edge_length / 2
y_3 = 0
x_1 = 0
x_2 = 0
x_3 = - np.sqrt(3) * edge_length / 2

# PONTOONS (these remain the same for both sets of pontoons, the ones above the water and the ones below the water)
# N° 1 is between floater 1 and 2; N° 2 is between floater 2 and 3
# y_p1 = 0
# y_p2 = (edge_length / 2) * np.sin(np.deg2rad(30))
# y_p3 = - y_p2
# x_p1 = 0
# x_p2 = -(edge_length / 2) * np.sin(np.deg2rad(60))
# x_p3 = x_p2

y_p1 = 0
y_p2 = (edge_length / 2) * (1 - np.cos(np.deg2rad(60)))
y_p3 = - y_p2
x_p1 = 0
x_p2 = (edge_length / 2) * np.sin(np.deg2rad(60))
x_p3 = x_p2

#%% CM COORDINATES    

z_CM_tow = (z_hub - h_out) / 2        # z center of mass tower

z_CM_bf = h_out - h_bf / 2            # z center of back floater
z_CM_ff = h_out - h_bf / 2            # z center of front floater

z_CM_bb = (h_bb / 2) - draft          # z center of back ballast
z_CM_fb = (h_fb / 2) - draft          # z center of front ballast

z_CM_pon_below = - D_pon / 2                # z center of pontoon, considering them - below the water level

z_CM_pon_above =  h_pontoon_above  # z center of pontoon, considering them - above the water level
                                      # immediately below the water level
z_CM_hp = - draft + h_hp / 2          # z center of heave plates

# z position of center of mass
z_CM_tot = (2 * z_CM_tow * m_tow + 2 * z_hub * m_rot + 2 * z_CM_bf * m_bf +\
            2 * z_CM_bb * m_bb + z_CM_ff * m_ff + z_CM_fb * m_fb +\
            3 * z_CM_pon_below * m_pon + 3 * z_CM_pon_above * m_pon + 3 * z_CM_hp * m_hp) / m_tot

print(f'The center of mass is at {z_CM_tot} m')

# y position center of mass
x_CM_tot = (m_tot_back * (x_1 + x_2) + m_tot_front  * x_3 +\
            2* m_pon *(x_p1 + x_p2 + x_p3)) / m_tot

print(f'The center of mass around y is at {x_CM_tot} m')

# x position center of mass
y_CM_tot = (m_tot_back * (y_1 + y_2) + m_tot_front * y_3 +\
            2* m_pon *(y_p1 + y_p2 + y_p3)) / m_tot

print(f'The center of mass around x is at {y_CM_tot} m')

#%% CB COORDINATES

# z position of center of buoyancy
z_CB_tot = - draft / 2

print(f'The center of bouyancy is at {z_CB_tot} m')

# submerged volumes
Vol_1_sub = A_back * (draft - h_hp)              # volume floater 1
Vol_2_sub = A_back * (draft - h_hp)              # volume floater 2
Vol_3_sub = A_front * (draft - h_hp)             # volume floater 3
Vol_hp = A_hp * h_hp                             # volume heave plates
Vol_pon = np.pi * ((D_pon / 2)**2) * L_pon       # volume pontoons

# x position of center of mass
x_CB_tot = (Vol_1_sub * x_1 + Vol_2_sub * x_2 + Vol_3_sub * x_3 +\
            Vol_pon * (x_p1 + x_p2 + x_p3) + Vol_hp * (x_1 + x_2 + x_3)) /\
            Vol_displaced
                    
print(f'The center of bouyancy around y is at {x_CB_tot} m')
      
# y position of center of mass              
y_CB_tot = (Vol_1_sub * y_1 + Vol_2_sub * y_2 + Vol_3_sub * y_3 +\
            Vol_pon * (y_p1 + y_p2 + y_p3) + Vol_hp * (y_1 + y_2 + y_3)) /\
            Vol_displaced

print(f'The center of bouyancy around x is at {y_CB_tot} m')

#%% INERTIA CALCULATION FOR MASS AND ADDED MASS MATRICES

# turbine: considered as point mass. first approximation (actually is an incline tower)
I0_XX_turb =  m_rot * (y_1**2 + z_hub**2)
I0_YY_turb = m_rot * (x_1**2 + z_hub**2)
I0_ZZ_turb = m_rot * (x_1**2 + y_1**2)
I0_XY_turb_1 = - m_rot * y_1 * x_1
I0_XY_turb_2 = - m_rot * y_2 * x_2
I0_XZ_turb_1 = - m_rot * x_1 * z_hub
I0_XZ_turb_2 = - m_rot * x_2 * z_hub
I0_YZ_turb_1 =  - m_rot * y_1 * z_hub
I0_YZ_turb_2 = - m_rot * y_2 * z_hub

# tower: it's a thick-walled cylinder
I0_XX_tower = m_tow * (3 * ((D_tow / 2)**2 + (D_tow_in / 2)**2)+ h_left**2) / 12 +\
                m_tow * (y_1**2 + z_CM_tow**2)
I0_YY_tower = m_tow * (3 * ((D_tow / 2)**2 + (D_tow_in / 2)**2)+ h_left**2) / 12 +\
              m_tow * (x_1**2 + z_CM_tow**2)
I0_ZZ_tower = 0.5 *m_tow * ((D_tow / 2)**2 + (D_tow_in / 2)**2) +\
              m_tow * (x_1**2 + y_1**2)  
I0_XY_tower_1 = - m_tow * x_1 * y_1
I0_XY_tower_2 = - m_tow * x_2 * y_2
I0_XZ_tower_1 = - m_tow * x_1 * z_CM_tow
I0_XZ_tower_2 = - m_tow * x_2 * z_CM_tow
I0_YZ_tower_1 = - m_tow * y_1 * z_CM_tow
I0_YZ_tower_2 = - m_tow * y_2 * z_CM_tow

# back floater: it's a thick-walled cylinder
I0_XX_backfloat =  m_bf * (3 * ((D_b / 2)**2 + (D_b_in / 2)**2)+ h_bf**2) / 12 +\
                   m_bf * (y_1**2 + z_CM_bf**2) 
I0_YY_backfloat =  m_bf * (3 * ((D_b / 2)**2 + (D_b_in / 2)**2)+ h_bf**2) / 12 +\
                   m_bf * (x_1**2 + z_CM_bf**2)
I0_ZZ_backfloat = 0.5 *m_bf * ((D_b / 2)**2 + (D_b_in / 2)**2) +\
                   m_bf * (x_1**2 + y_1**2)
I0_XY_backfloat_1 = - m_bf * x_1 * y_1
I0_XY_backfloat_2 = - m_bf * x_2 * y_2
I0_XZ_backfloat_1 = - m_bf * x_1 * z_CM_bf
I0_XZ_backfloat_2 = - m_bf * x_2 * z_CM_bf
I0_YZ_backfloat_1 = - m_bf * y_1 * z_CM_bf
I0_YZ_backfloat_2 = - m_bf * y_2 * z_CM_bf

# front floater: it's a thick-walled cylinder
I0_XX_frontfloat =  m_ff * (3 * ((D_f / 2)**2 + (D_f_in / 2)**2)+ h_ff**2) / 12 +\
                    m_ff * (y_3**2 + z_CM_ff**2) 
I0_YY_frontfloat =  m_ff * (3 * ((D_f / 2)**2 + (D_f_in / 2)**2)+ h_ff**2) / 12 +\
                    m_ff * (x_3**2 + z_CM_ff**2) 
I0_ZZ_frontfloat = 0.5 *m_ff * ((D_f / 2)**2 + (D_f_in / 2)**2) +\
                    m_ff * (x_3**2 + y_3**2) 
I0_XY_frontfloat = - m_ff * x_3 * y_3
I0_XZ_frontfloat = - m_ff * x_3 * z_CM_ff
I0_YZ_frontfloat = - m_ff * y_3 * z_CM_ff

# back ballast: it's a solid cylinder
I0_XX_backball =  m_bb * (3 * (D_b / 2)**2+ h_bb**2) / 12 +\
                  m_bb * (y_1**2 + z_CM_bb**2)  
I0_YY_backball =  m_bb * (3 * (D_b / 2)**2+ h_bb**2) / 12 +\
                  m_bb * (x_1**2 + z_CM_bb**2)  
I0_ZZ_backball = 0.5 *m_bb * (D_b / 2)**2 +\
                  m_bb * (x_1**2 + y_1**2)  
I0_XY_backball_1 = - m_bb * x_1 * y_1
I0_XY_backball_2 = - m_bb * x_2 * y_2
I0_XZ_backball_1 = - m_bb * x_1 * z_CM_bb
I0_XZ_backball_2 = - m_bb * x_2 * z_CM_bb
I0_YZ_backball_1 = - m_bb * y_1 * z_CM_bb
I0_YZ_backball_2 = - m_bb * y_2 * z_CM_bb

# front ballast: it's a solid cylinder
I0_XX_frontball =  m_fb * (3 * (D_f / 2)**2 + h_fb**2) / 12 +\
                   m_fb * (y_3**2 + z_CM_fb**2) 
I0_YY_frontball =  m_fb * (3 * (D_f / 2)**2 + h_fb**2) / 12 +\
                   m_fb * (x_3**2 + z_CM_fb**2) 
I0_ZZ_frontball = 0.5 *m_fb * (D_f / 2)**2 +\
                   m_fb * (x_3**2 + y_3**2) 
I0_XY_frontball = - m_fb * x_3 * y_3
I0_XZ_frontball = - m_fb * x_3 * z_CM_fb
I0_YZ_frontball = - m_fb * y_3 * z_CM_fb

# back heave plate: it's a solid cylinder
I0_XX_backhp =  m_hp * (3 * (D_hp / 2)**2 + h_hp**2) / 12 +\
                m_hp * (y_1**2 + z_CM_hp**2)
I0_YY_backhp =  m_hp * (3 * (D_hp / 2)**2 + h_hp**2) / 12 +\
                m_hp * (x_1**2 + z_CM_hp**2)
I0_ZZ_backhp = 0.5 *m_hp * (D_hp / 2)**2 +\
                m_hp * (x_1**2 + y_1**2)
I0_XY_backhp_1 = - m_hp * x_1 * y_1
I0_XY_backhp_2 = - m_hp * x_2 * y_2
I0_XZ_backhp_1 = - m_hp * x_1 * z_CM_hp
I0_XZ_backhp_2 = - m_hp * x_2 * z_CM_hp
I0_YZ_backhp_1 = - m_hp * y_1 * z_CM_hp
I0_YZ_backhp_2 = - m_hp * y_2 * z_CM_hp

# front heave plate: it's a solid cylinder
I0_XX_fronthp =  m_hp * (3 * (D_hp / 2)**2 + h_hp**2) / 12 +\
                 m_hp * (y_3**2 + z_CM_hp**2)
I0_YY_fronthp =  m_hp * (3 * (D_hp / 2)**2+ h_hp**2) / 12 +\
                 m_hp * (x_3**2 + z_CM_hp**2)
I0_ZZ_fronthp = 0.5 *m_hp * (D_hp / 2)**2 +\
                 m_hp * (x_3**2 + y_3**2)
I0_XY_fronthp = - m_hp * x_3 * y_3
I0_XZ_fronthp = - m_hp * x_3 * z_CM_hp
I0_YZ_fronthp = - m_hp * y_3 * z_CM_hp

# Pontoon 1: it's a horizontal thick-walled cylinder
I0_XX_pontoon_1 =  m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12 +\
                   m_pon * (y_p1**2 + z_CM_pon_below**2) 
I0_YY_pontoon_1 =  m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12 +\
                   m_pon * (x_p1**2 + z_CM_pon_below**2)  
I0_ZZ_pontoon_1 = 0.5 *m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2) +\
                   m_pon * (x_p1**2 + y_p1**2)  
I0_XY_pontoon_1 = - m_pon * x_p1 * y_p1
I0_XZ_pontoon_1 = - m_pon * x_p1 * z_CM_pon_below
I0_YZ_pontoon_1 = - m_pon * y_p1 * z_CM_pon_below

theta_2 = +30
# Pontoon 2: it's a horizontal thick-walled cylinder
I0_XX_pontoon_2 = (m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) * (np.sin(np.deg2rad(theta_2))**2) +\
                  (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2)) * (np.cos(np.deg2rad(theta_2))**2) +\
                  m_pon * (y_p2**2 + z_CM_pon_below**2)
I0_YY_pontoon_2 = (m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) * (np.sin(np.deg2rad(theta_2))**2) +\
                  (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2)) * (np.cos(np.deg2rad(theta_2))**2) +\
                  m_pon * (x_p2**2 + z_CM_pon_below**2) 
I0_ZZ_pontoon_2 = m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2) + L_pon**2) / 12 +\
                  m_pon * (x_p2**2 + y_p2**2)  
I0_XY_pontoon_2 =  ((m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) -\
                   (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2))) *\
                   np.sin(np.deg2rad(theta_2)) * np.cos(np.deg2rad(theta_2)) -\
                   m_pon * x_p2 * y_p2
I0_XZ_pontoon_2 = - m_pon * x_p2 * z_CM_pon_below
I0_YZ_pontoon_2 = - m_pon * y_p2 * z_CM_pon_below

theta_3 = -30
# Pontoon 3: it's a horizontal thick-walled cylinder
I0_XX_pontoon_3 = (m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) * (np.sin(np.deg2rad(theta_3))**2) +\
                  (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2)) * (np.cos(np.deg2rad(theta_3))**2) +\
                  m_pon * (y_p3**2 + z_CM_pon_below**2)
I0_YY_pontoon_3 = (m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) * (np.sin(np.deg2rad(theta_3))**2) +\
                  (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2)) * (np.cos(np.deg2rad(theta_3))**2) +\
                  m_pon * (x_p3**2 + z_CM_pon_below**2) 
I0_ZZ_pontoon_3 = m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2) + L_pon**2) / 12 +\
                  m_pon * (x_p3**2 + y_p3**2)  
I0_XY_pontoon_3 = ((m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) -\
                   (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2))) *\
                   np.sin(np.deg2rad(theta_3)) * np.cos(np.deg2rad(theta_3)) -\
                   m_pon * x_p3 * y_p3
I0_XZ_pontoon_3 = - m_pon * x_p3 * z_CM_pon_below
I0_YZ_pontoon_3 = - m_pon * y_p3 * z_CM_pon_below


# Pontoon 1 (upper set): it's a horizontal thick-walled cylindee
I0_XX_pontoon_1_up =  m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12 +\
                   m_pon * (y_p1**2 + z_CM_pon_above**2) 
I0_YY_pontoon_1_up =  m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12 +\
                   m_pon * (x_p1**2 + z_CM_pon_above**2)  
I0_ZZ_pontoon_1_up = 0.5 *m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2) +\
                   m_pon * (x_p1**2 + y_p1**2)  
I0_XY_pontoon_1_up = - m_pon * x_p1 * y_p1
I0_XZ_pontoon_1_up = - m_pon * x_p1 * z_CM_pon_above
I0_YZ_pontoon_1_up = - m_pon * y_p1 * z_CM_pon_above

theta_2 = +30
# Pontoon 2 (upper set): it's a horizontal thick-walled cylinder
I0_XX_pontoon_2_up = (m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) * (np.sin(np.deg2rad(theta_2))**2) +\
                  (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2)) * (np.cos(np.deg2rad(theta_2))**2) +\
                  m_pon * (y_p2**2 + z_CM_pon_above**2)
I0_YY_pontoon_2_up = (m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) * (np.sin(np.deg2rad(theta_2))**2) +\
                  (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2)) * (np.cos(np.deg2rad(theta_2))**2) +\
                  m_pon * (x_p2**2 + z_CM_pon_above**2) 
I0_ZZ_pontoon_2_up = m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2) + L_pon**2) / 12 +\
                  m_pon * (x_p2**2 + y_p2**2)  
I0_XY_pontoon_2_up =  ((m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) -\
                   (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2))) *\
                   np.sin(np.deg2rad(theta_2)) * np.cos(np.deg2rad(theta_2)) -\
                   m_pon * x_p2 * y_p2
I0_XZ_pontoon_2_up = - m_pon * x_p2 * z_CM_pon_above
I0_YZ_pontoon_2_up = - m_pon * y_p2 * z_CM_pon_above

theta_3 = -30
# Pontoon 3 (upper set): it's a horizontal thick-walled cylinder
I0_XX_pontoon_3_up = (m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) * (np.sin(np.deg2rad(theta_3))**2) +\
                  (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2)) * (np.cos(np.deg2rad(theta_3))**2) +\
                  m_pon * (y_p3**2 + z_CM_pon_above**2)
I0_YY_pontoon_3_up = (m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) * (np.sin(np.deg2rad(theta_3))**2) +\
                  (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2)) * (np.cos(np.deg2rad(theta_3))**2) +\
                  m_pon * (x_p3**2 + z_CM_pon_above**2) 
I0_ZZ_pontoon_3_up = m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2) + L_pon**2) / 12 +\
                  m_pon * (x_p3**2 + y_p3**2)  
I0_XY_pontoon_3_up = ((m_pon * (3 * ((D_pon / 2)**2 + (D_pon_in / 2)**2)+ L_pon**2) / 12) -\
                   (0.5 * m_pon * ((D_pon / 2)**2 + (D_pon_in / 2)**2))) *\
                   np.sin(np.deg2rad(theta_3)) * np.cos(np.deg2rad(theta_3)) -\
                   m_pon * x_p3 * y_p3
I0_XZ_pontoon_3_up = - m_pon * x_p3 * z_CM_pon_above
I0_YZ_pontoon_3_up = - m_pon * y_p3 * z_CM_pon_above

# sum of inertia components
I0_XX = 2*(I0_XX_turb + I0_XX_tower + I0_XX_backball + I0_XX_backfloat + I0_XX_backhp) +\
        I0_XX_frontfloat + I0_XX_frontball + I0_XX_fronthp +\
        I0_XX_pontoon_1 + I0_XX_pontoon_2 + I0_XX_pontoon_3 +\
        I0_XX_pontoon_1_up + I0_XX_pontoon_2_up + I0_XX_pontoon_3_up
I0_YY = 2*(I0_YY_turb + I0_YY_tower + I0_YY_backball + I0_YY_backfloat + I0_YY_backhp) +\
        I0_YY_frontfloat + I0_YY_frontball + I0_YY_fronthp +\
        I0_YY_pontoon_1 + I0_YY_pontoon_2 + I0_YY_pontoon_3 +\
        I0_YY_pontoon_1_up + I0_YY_pontoon_2_up + I0_YY_pontoon_3_up
I0_ZZ = 2*(I0_ZZ_turb + I0_ZZ_tower + I0_ZZ_backball + I0_ZZ_backfloat + I0_ZZ_backhp) +\
        I0_ZZ_frontfloat + I0_ZZ_frontball + I0_ZZ_fronthp +\
        I0_ZZ_pontoon_1 + I0_ZZ_pontoon_2 + I0_ZZ_pontoon_3 +\
        I0_ZZ_pontoon_1_up + I0_ZZ_pontoon_2_up + I0_ZZ_pontoon_3_up
I0_XY = I0_XY_turb_1 + I0_XY_tower_1 + I0_XY_backball_1 + I0_XY_backfloat_1 + I0_XY_backhp_1 +\
        I0_XY_turb_2 + I0_XY_tower_2 + I0_XY_backball_2 + I0_XY_backfloat_2 + I0_XY_backhp_2 +\
        I0_XY_frontfloat + I0_XY_frontball + I0_XY_fronthp +\
        I0_XY_pontoon_1 + I0_XY_pontoon_2 + I0_XY_pontoon_3 +\
        I0_XY_pontoon_1_up + I0_XY_pontoon_2_up + I0_XY_pontoon_3_up
I0_XZ = I0_XZ_turb_1 + I0_XZ_tower_1 + I0_XZ_backball_1 + I0_XZ_backfloat_1 + I0_XZ_backhp_1 +\
        I0_XZ_turb_2 + I0_XZ_tower_2 + I0_XZ_backball_2 + I0_XZ_backfloat_2 + I0_XZ_backhp_2 +\
        I0_XZ_frontfloat + I0_XZ_frontball + I0_XZ_fronthp +\
        I0_XZ_pontoon_1 + I0_XZ_pontoon_2 + I0_XZ_pontoon_3 +\
        I0_XZ_pontoon_1_up + I0_XZ_pontoon_2_up + I0_XZ_pontoon_3_up
I0_YZ = I0_YZ_turb_1 + I0_YZ_tower_1 + I0_YZ_backball_1 + I0_YZ_backfloat_1 + I0_YZ_backhp_1 +\
        I0_YZ_turb_2 + I0_YZ_tower_2 + I0_YZ_backball_2 + I0_YZ_backfloat_2 + I0_YZ_backhp_2 +\
        I0_YZ_frontfloat + I0_YZ_frontball + I0_YZ_fronthp +\
        I0_YZ_pontoon_1 + I0_YZ_pontoon_2 + I0_YZ_pontoon_3 +\
        I0_YZ_pontoon_1_up + I0_YZ_pontoon_2_up + I0_YZ_pontoon_3_up

# force the symmetric terms around the X axis (I0_XY and I0_YZ) to be zero,
# otherwise they could have some floating error
I0_YZ = 0
I0_XY = 0
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

A = M

# Heave coupled components of added mass matrix changing because of heave plate
m_a_heave_plate = (1/3) * rhow * (D_hp**3) # For one heave plate 
                                           # (According to L.Tao et.al.--> Eq. 24)
m_a_total = 3 * m_a_heave_plate            # For three heave plates, assuming
                                           # that the heave plates are seperated
                                           # such that they don't interact with
                                           # each other

A33 = m_a_total                         # Heave-heave
A34 = m_a_total * (y_1 + y_2 + y_3)     # Heave-roll
A35 = - m_a_total * (x_1 + x_2 + x_3)   # Heave-pitch
A43 = A34                               # Roll-heave
A53 = A35                               # Pitch-heave

# Rewrite surge and sway components with columns and pontoons contribution
A11 = 2 * rhow * np.pi * D_b**2 / 4 * (draft - h_hp) +\
      rhow * np.pi * D_f**2 / 4 * (draft - h_hp) +\
      rhow * np.pi * D_pon**2 / 4 * L_pon * (1 + 2 * np.sin(np.deg2rad(60))**2)
A22 = 2 * rhow * np.pi * D_b**2 / 4 * (draft - h_hp) +\
      rhow * np.pi * D_f**2 / 4 * (draft - h_hp) +\
      rhow * np.pi * D_pon**2 / 4 * L_pon * (0 + 2 * np.sin(np.deg2rad(30))**2)
A33 += 3 * rhow * np.pi * D_pon**2 / 4 * L_pon


A[2,2] = A33
A[2,3] = A34
A[2,4] = A35
A[3,2] = A43
A[4,2] = A53
#%% DAMPING MATRIX

B = np.zeros_like(M)

#%% MOMENTS OF WATERPLANE AREA

# Linear moments
I_x = 2 * A_back * x_1 +  A_front * x_3

I_y = A_back * y_1 + A_back * y_2 + A_front * y_3

# Second-order moments
I_xx = 2 * np.pi * D_b**4 / 64 + np.pi * D_f**4 / 64 +\
                    2 * A_back * x_1**2 +\
                        A_front * x_3**2
                        
I_yy = 2 * np.pi * D_b**4 / 64 + np.pi * D_f**4 / 64 +\
                    A_back * y_1**2 + A_back * y_2**2 +\
                    A_front * y_3**2
                    
I_xy = A_back * y_1 * x_1 + A_back * y_2 * x_2 +\
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

#%% MOORING STIFFNESS MATRIX

# surge
C_mooring_11 = 1.26538165e+02   # surge-surge
C_mooring_12 = -2.61478562e-01  # surge-sway
C_mooring_13 = 1.21770107e+00   # surge-heave
C_mooring_14 = 2.15752957e-01   # surge-roll
C_mooring_15 = 5.80170115e+00   # surge-pitch
C_mooring_16 = -7.84234243e-01  # surge-yaw

# sway
C_mooring_21 = -8.88178420e-15  # sway-surge
C_mooring_22 = 1.26776086e+02   # sway-sway
C_mooring_23 = 2.62712856e+00   # sway-heave
C_mooring_24 = -8.82064317e+00  # sway-roll
C_mooring_25 = 3.02903601e+00   # sway-pitch
C_mooring_26 = 2.55817695e-02   # sway-yaw

# heave
C_mooring_31 = 1.47573058e+00   # heave-surge
C_mooring_32 = -1.35542875e-01  # heave-sway
C_mooring_33 = 1.67204365e+02   # heave-heave
C_mooring_34 = 1.61298010e+00   # heave-roll
C_mooring_35 = 3.19502978e+00   # heave-pitch
C_mooring_36 = -2.01957624e+00  # heave-yaw

# roll
C_mooring_41 = -0.00000000e+00  # roll-surge
C_mooring_42 = -8.56427108e+00  # roll-sway
C_mooring_43 = -2.90021400e-01  # roll-heave
C_mooring_44 = 2.92091256e+00   # roll-roll
C_mooring_45 = -3.85455962e-01  # roll-pitch
C_mooring_46 = 4.13293631e-03   # roll-yaw

# pitch
C_mooring_51 = 8.53327833e+00   # pitch-surge
C_mooring_52 = -3.00750034e-02  # pitch-sway
C_mooring_53 = 1.05623851e-01   # pitch-heave
C_mooring_54 = 3.63873979e-02   # pitch-roll
C_mooring_55 = 2.54547112e+00   # pitch-pitch
C_mooring_56 = -8.70129687e-02  # pitch-yaw

# yaw
C_mooring_61 = -3.46944695e-17  # yaw-surge
C_mooring_62 = 3.41095414e-02   # yaw-sway
C_mooring_63 = -1.78247970e-02  # yaw-heave
C_mooring_64 = -2.25648954e-03  # yaw-roll
C_mooring_65 = -8.07955243e-04  # yaw-pitch
C_mooring_66 = 4.24727171e+00   # yaw-yaw

# matrix assembly
C_mooring = np.array([[C_mooring_11, C_mooring_21, C_mooring_31, C_mooring_41, C_mooring_51, C_mooring_61],
                      [C_mooring_12, C_mooring_22, C_mooring_32, C_mooring_42, C_mooring_52, C_mooring_62],
                      [C_mooring_13, C_mooring_23, C_mooring_33, C_mooring_43, C_mooring_53, C_mooring_63],
                      [C_mooring_14, C_mooring_24, C_mooring_34, C_mooring_44, C_mooring_54, C_mooring_64],
                      [C_mooring_15, C_mooring_25, C_mooring_35, C_mooring_45, C_mooring_55, C_mooring_65],
                      [C_mooring_16, C_mooring_26, C_mooring_36, C_mooring_46, C_mooring_56, C_mooring_66]])


C_tot = C_hydro + C_mooring

#%% NATURAL FREQUENCIES EVALUATION


CoMA = np.linalg.solve(M + A, C_tot)
eigVal, eigVec = np.linalg.eig(CoMA)

for i in range(6):
    print(f"\nMode {i}")
    print(eigVec[:, i])

# eliminate floating errors
eigVal[np.abs(eigVal) < 1e-10] = 0

# Natural frequencies
omega_nat = np.sqrt(eigVal)
f_nat = omega_nat / (2 * np.pi)

# Display natural frequencies
print(f'Surge period: {f_nat[3]:.2f} [Hz]')
print(f'Sway period: {f_nat[4]:.2f} [Hz]')
print(f'Heave period: {f_nat[0]:.2f} [Hz]')
print(f'Roll period: {f_nat[5]:.2f} [Hz]')
print(f'Pitch period: {f_nat[1]:.2f} [Hz]')
print(f'Yaw period: {f_nat[2]:.2f} [Hz]')





