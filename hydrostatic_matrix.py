'''This code is for a fast dimensioning of a semi-sub triangular floater
to better understand the imput mesurement you can refer to the README file'''

# %%

#imports
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

I_y = A_back *y_1 + A_back * y_2 + A_front * y_3
I_x = 2 * A_back * x_1 +  A_front * x_3

C_34 = rhow * g * I_y       # heave-roll
C_35 = - rhow * g * I_x     # heave-pitch
C_36 = 0       # heave-yaw

# roll
C_41 = 0       # roll-surge
C_42 = 0       # roll-sway
C_43 = C_34    # roll-heave

I_yy = 2 * np.pi * D_b**4 / 64 + np.pi * D_f**4 / 64 +\
                    A_back * y_1**2 + A_back * y_2**2 +\
                    A_front * y_3**2
I_xy = 2 * np.pi * D_b**4 / 64 + np.pi * D_f**4 / 64 +\
                    A_back * y_1 * x_1 + A_back * y_2 * x_2 +\
                    A_front * y_3 * x_3

C_44 =  rhow * g * I_yy + m_tot * g * (z_CB_tot - z_CM_tot) # roll-roll
C_45 = - rhow * g * I_xy   # roll-pitch
C_46 = 0       # roll-yaw
 
# pitch
C_51 = 0       # pitch-surge
C_52 = 0       # pitch-sway
C_53 = C_35    # pitch-heave
C_54 = C_45    # pitch-roll

I_xx = 2 * np.pi * D_b**4 / 64 + np.pi * D_f**4 / 64 +\
                    2 * A_back * x_1**2 +\
                        A_front * x_3**2
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