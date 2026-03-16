'''This code is for a fast dimensioning of a semi-sub triangular floater
to better understand the imput mesurement you can refer to the README file'''

# %%

from scipy.linalg import eig
import numpy as np

## ambient parameters

g = 9.81            # Gravity [m/s^2]
h = 2               # Water depth [m]
rhow = 1025         # Water density [kg/m^3]
rhoa = 1.29         # Air density [kg/m^3]

# tower and rotor parameter
mt = 3              # Tower mass [kg]
z_hub = 0.8         # tower length [m]

Drot = 1            # Rotor diameter [m]
Vrated = 11.4       # Rated wind speed [m/s]
Trated = 900        # Thrust at rated [N]

mn = 3          # Nacelle mass [kg]
mr = 1            # Rotor mass [kg]

# floater parameters

m_bf = 5          # Back Floater mass [kg]
m_ff = 5            # Front Floater mass [kg]
h_bf = 1            # Back floater hight [m]
h_ff = 1            # Front floater hight [m] 
D_b = 0.5           # Back floater Diameter [m]
D_f = 0.6           # Front floater Diameter [m]

m_bb = 8            # Back Ballast mass [kg]
m_fb = 20            # Front Ballast mass [kg]

h_bb = 0.5             # Back Ballast height [m]
h_fb = 0.5             # Front Ballast height [m] 


edge = 1           # distance between columns: equi lateral triangle

# mooring parameter 

Kmoor = 12          # Mooring stiffness [N/m]
zmoor = -70         # Mooring fairlead [m]

#others

#Cm = 1.0            # Added mass coefficient
#Cd = 0.6            # Drag coefficient

 
#%% Calculate the draft 

m_back =  2 * (mt + mn + mr + m_bb + m_bf)   # Back side mass [kg] 

m_front = m_fb + m_ff                        # Front side mass [kg] 

m_tot = m_back + m_front                     # Total mass

A_back = np.pi * D_b**2 / 4                  # Sectional area back cylinders [m^2]

A_front = np.pi * D_f**2 / 4                 # Sectional area front cylinder [m^2]

A_total = 2 * A_back + A_front               # total sectional area [m^2]

draft = m_tot / (A_total*rhow) #draft [m]    # draft [m]

print(f'The draft is {draft} m ')

h_out =  (h_bf +  h_bb) - draft              # height of the floater outside water

#%% Determine positions of columns based on centroid. The x axis is pointing
#   from the front floater to the ones in the back, in the direction of the
#   wind. The y axis follows the right hand rule.    

# The back floaters are N° 1 and N° 2, while the front is N°3
y_1 = - edge / 2
y_2 = + edge / 2
y_3 = 0
x_1 = + edge / 3
x_2 = x_1
x_3 =  - 2 * edge / 3

#%% Calculate center of mass coordinates

z_CM_tow = (z_hub / 2 + h_out)        # z center of mass tower

z_CM_bf = h_out - h_bf / 2            # z center of back floater
z_CM_ff = h_out - h_bf / 2            # z center of front floater

z_CM_bb = (h_bb / 2) - draft          # z center of back ballast
z_CM_fb = (h_fb / 2) - draft          # z center of front ballast

# z position of center of mass
z_CM_tot = (2 * z_CM_tow * mt + 2 * (z_hub + h_out) * (mr + mn) + 2 * z_CM_bf * m_bf +\
            2 * z_CM_bb * m_bb + z_CM_ff * m_ff + z_CM_fb * m_fb) / m_tot

# x position center of mass
x_CM_tot = ((mr + mn + mt + m_bf + m_bb) * x_1 + (mr + mn + mt + m_bf + m_bb) * x_2 + (m_ff + m_fb) *x_3) / m_tot

# y position center of mass
y_CM_tot = ((mr + mn + mt + m_bf + m_bb) * y_1 + (mr + mn + mt + m_bf + m_bb) * y_2 + (m_ff + m_fb) *y_3) / m_tot

print(f'The center of mass is at {z_CM_tot} m')

# y position center of mass
x_CM_tot = ((mr + mn + mt + m_bf + m_bb) * x_1 + (mr + mn + mt + m_bf + m_bb) * x_2\
            + (m_ff + m_fb) *x_3) / m_tot

print(f'The center of mass around y is at {x_CM_tot} m')

# x position center of mass
y_CM_tot = ((mr + mn + mt + m_bf + m_bb) * y_1 + (mr + mn + mt + m_bf + m_bb) * y_2\
            + (m_ff + m_fb) *y_3) / m_tot

print(f'The center of mass around x is at {y_CM_tot} m')

#%% Calculate center of bouyancy coordinates

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

#%% Hydrostatic restoring matrix
# surge
C_11 = 0       # surge-surge
C_12 = 0       # surge-heave
C_13 = 0       # surge-pitch
C_14= 0        # surge-yaw

# heave
C_21 = 0       # heave-surge
C_22 = 2 * rhow * g * A_back + rhow * g * A_front   # heave-heave

I_x = -rhow * g * (2 * A_back * x_1 +  A_front * x_3)
C_23 = I_x     # heave-pitch
C_24 = 0       # heave-yaw

# pitch
C_31 = 0       # pitch-surge
C_32 = C_23    # pitch-heave

I_xx = rhow * g * (2 * np.pi * D_b**4 / 64 + np.pi * D_f**4 / 64 +\
                    2 * A_back * x_1**2 +\
                        A_front * x_3**2)
C_33 = m_tot * g * (z_CB_tot - z_CM_tot) + I_xx   # pitch-pitch
C_34 = 0       # pitch-yaw
  
# yaw
C_41 = 0       # yaw-surge
C_42 = 0       # yaw-heave
C_43 = - m_tot * g * (y_CB_tot - y_CM_tot)   # yaw-pitch (if buyonacy and mass are not aligned in y direction there is a restoring moment in yaw)
C_44 = 0       # yaw-yaw


C_hydro = np.array([[C_11, C_21, C_31, C_41],
                    [C_12, C_22, C_32, C_42],
                    [C_13, C_23, C_33, C_43],
                    [C_14, C_24, C_34, C_44]])


# Mass matrix calculation
I0turb = (mn+mr)*z_hub**2
I0tower =  mt * z_CM_tow**2
I0float =  m_ff * z_CM_ff**2 + 2 * m_bf * z_CM_bf**2

IOtot = I0float + I0tower + I0turb
M11 = m_tot # surge-surge
M12 = 0 # surge-heave
M13 = m_tot * z_CM_tot # surge-pitch
M14 = -m_tot * y_CM_tot # surge-yaw
M21 = 0 # heave-surge
M22 = m_tot # heave-heave
M23 = -m_tot * x_CM_tot # heave-pitch
M24 = 0 # heave-yaw
M31 = m_tot * z_CM_tot # pitch-surge
M32 = -m_tot * x_CM_tot # pitch-heave
M33 = IOtot # pitch-pitch
M34 = 0 # pitch-yaw (considering xg = yg = 0)
M41 = -m_tot * y_CM_tot # yaw-surge
M42 = 0 # yaw-heave
M43 = 0 # yaw-pitch (considering xg = yg = 0)
M44 = 2*IOtot # yaw-yaw 

M = np.array([[M11, M12, M13, M14],
                    [M21, M22, M23, M24],
                    [M31, M32, M33, M34],
                    [M41, M42, M43, M44]])


#-----------------------------------------------
# Only surge and pitch are considered in added mass matrices
#-----------------------------------------------
# Added mass matrix
D1 = 240   # back left column based on previous year design
D2 = 240  # back right column based on previous year design
D3 = 240  # front column based on previous year design

z_bot = -draft
Cm = 2
#From A5 on offshore there was used only one cylinder (one diameter), here is the sum of the three cylinders
# --- A11 (surge added mass)
A11 = -rhow * (np.pi/4) * Cm * (z_bot * D1**2 + z_bot * D2**2 + z_bot * D3**2)

# --- A15 (surge-pitch coupling)
A15 = -rhow * (np.pi/4) * Cm * (z_bot * D1**2 * (z_bot/2) + z_bot * D2**2 * (z_bot/2) +z_bot * D3**2 * (z_bot/2))

#   symmetry
A51 = -rhow * (np.pi/4) * Cm * (z_bot * D1**2 * (z_bot/2) + z_bot * D2**2 * (z_bot/2) +z_bot * D3**2 * (z_bot/2))

# --- A55 (pitch added inertia)
A55 = -rhow * (np.pi/4) * Cm * (z_bot * D1**2 * (z_bot**2/3) +z_bot * D2**2 * (z_bot**2/3) +z_bot * D3**2 * (z_bot**2/3))

A = np.array([[A11, A15],[A51, A55]])
# restoring matrix
#C = np.array([[0. ,0.],[0. , m_total*g*(zCM_tot - zCB) + rhow*g*IAA]])

# Damping matrix
#B = np.array([[B11, 0],[0 ,0]])

print("M: ",M)
print("A: ",A)