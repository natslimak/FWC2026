'''This code is for a fast dimensioning of a semi-sub triangular floater
to better understand the imput mesurement you can refer to the README file'''

# %%

from scipy.linalg import eig
import numpy as np
import sys

## ambient parameters

g = 9.81            # Gravity [m/s^2]
h = 2               # Water depth [m]
rhow = 1025         # Water density [kg/m^3]
rhoa = 1.29         # Air density [kg/m^3]

# tower and rotor parameter
mt = 3              # Tower mass [kg]
z_hub = 1.2         # tower length [m]
D_t = 0.1           # tower diameter [m]

Drot = 1            # Rotor diameter [m]
Vrated = 11.4       # Rated wind speed [m/s]
Trated = 900        # Thrust at rated [N]

mn = 5.5          # Nacelle mass [kg]
mr = 1            # Rotor mass [kg]

# floater parameters

m_bf = 1          # Back Floater mass [kg]
m_ff = 1            # Front Floater mass [kg]
h_bf = 1            # Back floater hight [m]
h_ff = 1            # Front floater hight [m] 
D_b = 0.2           # Back floater Diameter [m]
D_f = 0.2           # Front floater Diameter [m]

m_bb = 8            # Back Ballast mass [kg]

edge = 1           # distance between columns: equi lateral triangle

# mooring parameter 

Kmoor = 12          # Mooring stiffness [N/m]
zmoor = -70         # Mooring fairlead [m]

rho_ballast = 1025         # Ballast density [kg/m^3]


# Calculate the draft 

m_back = 2* (mt + mn + mr + m_bb + m_bf)   # Back side mass [kg] 
 
A_back_floater = 2 * np.pi * D_b**2 /4 # Total area back floater [m^2]
A_back_single_floater = np.pi * D_b**2 /4

draft = m_back / (A_back_floater*rhow) #draft [m]


if draft >= h_bf: 

    sys.exit('WE ARE SINKING :(')

else:
    print ('WE ARE FLOATING :)')
    print ('The draft with', m_bb,'kg of ballast in the back is',draft,'m')
    



A_front_floater =  np.pi * D_f**2 /4

m_front = (rhow * A_front_floater * draft) 

m_fb = m_front - m_ff

print ('We need',m_fb ,'kg of ballast in front')

m_tot = m_back + m_front

print('Total weigth of the structure is', m_tot, 'kg')


#calculate the hight of the tower and the ballast

h_t = z_hub-(h_bf-draft) # Tower hight [m]
h_bb = m_bb/(rho_ballast* np.pi * D_b**2 /4) # back balast hight [m]
h_fb = m_fb/(rho_ballast* np.pi * D_f**2 /4) # back balast hight [m]

if  h_fb >= draft: 

    sys.exit('We need an heavier weigth')

print ('The tower shold be',h_t,'m')


h_out =  h_bf - draft              # height of the floater outside water (height of the ballast is included in height of the floater)

#%% Determine positions of columns based on centroid. The x axis is pointing
#   from the front floater to the ones in the back, in the direction of the
#   wind. The y axis follows the right hand rule.    


# The back floaters are N° 1 and N° 2, while the front is N°3
y_1 = - edge / 2
y_2 = + edge / 2
y_3 = 0
x_1 = 0
x_2 = 0
x_3 =  np.sqrt(3)*edge/2

#%% Calculate center of mass coordinates

z_CM_tow = (h_t / 2) + h_out        # z center of mass tower

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

Vol_1_sub = draft * A_back_single_floater
Vol_2_sub = draft * A_back_single_floater
Vol_3_sub = draft * A_front_floater

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
C_22 = 2 * rhow * g * A_back_single_floater + rhow * g * A_front_floater   # heave-heave

I_x = -rhow * g * (2 * A_back_single_floater * x_1 +  A_front_floater * x_3)
C_23 = I_x     # heave-pitch
C_24 = 0       # heave-yaw

# pitch
C_31 = 0       # pitch-surge
C_32 = C_23    # pitch-heave

I_xx = rhow * g * (2 * np.pi * D_b**4 / 64 + np.pi * D_f**4 / 64 +\
                    2 * A_back_single_floater * x_1**2 +\
                        A_front_floater * x_3**2)
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
#-----------------------------------------------
# All six degree of freedom
#-----------------------------------------------

#Inertia moment
#-----------------------------------------------
# Calculate all the componentsin the CM and transport it 
#  to the middle point of the back side(floating point)
#-----------------------------------------------
I0_XX_turb = 0 
I0_YY_turb = (mn+mr)*(edge/2)**2 
I0_ZZ_turb = (mn+mr)*z_hub**2

I0_XX_tower = mt*(3*(D_t/2)**2+ h_t**2)/12
I0_YY_tower = mt*(3*(D_t/2)**2+ h_t**2)/12 + mt * (edge/2)**2
I0_ZZ_tower = 0.5 *mt*(D_t/2)**2+ mt * z_CM_tow**2 

I0_XX_backfloat =  m_bf*(3*(D_b/2)**2+ h_bf**2)/12 
I0_YY_backfloat =  m_bf*(3*(D_b/2)**2+ h_bf**2)/12 + m_bf * (edge/2)**2
I0_ZZ_backfloat = 0.5 *m_bf*(D_b/2)**2+ m_bf * z_CM_bf**2 

I0_XX_backBallast =  m_bb*(3*(D_b/2)**2+ h_bb**2)/12 
I0_YY_backBallast =  m_bb*(3*(D_b/2)**2+ h_bb**2)/12 + m_bb * (edge/2)**2
I0_ZZ_backBallast = 0.5 *m_bb*(D_b/2)**2+ m_bb * z_CM_bb**2 

I0_XX_frontfloat =  m_ff*(3*(D_f/2)**2+ h_ff**2)/12 + m_ff * (np.sqrt(3)*edge/2)**2
I0_YY_frontfloat =  m_ff*(3*(D_f/2)**2+ h_ff**2)/12 
I0_ZZ_frontfloat = 0.5 *m_ff*(D_f/2)**2+ m_ff * z_CM_ff**2 

I0_XX_frontBallast =  m_fb*(3*(D_f/2)**2+ h_fb**2)/12 + m_fb * (np.sqrt(3)*edge/2)**2
I0_YY_frontBallast =  m_fb*(3*(D_f/2)**2+ h_fb**2)/12 
I0_ZZ_frontBallast = 0.5 *m_fb*(D_f/2)**2+ m_fb * z_CM_fb**2 

#-----------------------------------------------
# Sum all the inertia components 
#  Note: all the back components are doubled here
#-----------------------------------------------

I0_XX = 2*(I0_XX_turb + I0_XX_backBallast + I0_XX_backfloat + I0_XX_tower) +I0_XX_frontfloat + I0_XX_frontBallast 
I0_YY = 2*(I0_YY_turb + I0_YY_backBallast + I0_YY_backfloat + I0_YY_tower) +I0_YY_frontfloat + I0_YY_frontBallast 
I0_ZZ = 2*(I0_ZZ_turb + I0_ZZ_backBallast + I0_ZZ_backfloat + I0_ZZ_tower) +I0_ZZ_frontfloat + I0_ZZ_frontBallast 


#-----------------------------------------------
# Assemble the M matrice
#  Note: M is simmetrical
#-----------------------------------------------

M11 = m_tot # surge-surge
M12 = 0 # surge-sway
M13 = 0 # surge-heave
M14 = 0 # surge-roll
M15 = m_tot * z_CM_tot#surge-pitch
M16 = -m_tot * y_CM_tot #surge-ya

M21 = M12 # sway-surge
M22 = m_tot # sway-sway
M23 = 0 # sway-heave
M24 = -m_tot * z_CM_tot # sway-roll
M25 = 0 #sway-pitch
M26 = m_tot * x_CM_tot#sway-yaw

M31 = M13 # heave-surge
M32 = M23 # heave-sway
M33 = m_tot # heave-heave
M34 = m_tot*y_CM_tot # heave-roll
M35 = -m_tot*x_CM_tot  #heave-pitch
M36 = 0 #heave-yaw

M41 = M14 # roll-surge
M42 = M24 # roll-sway
M43 = M34 # roll-heave
M44 = I0_YY + I0_ZZ # roll-roll 
M45 = 0 #roll-pitch
M46 = m_tot*z_CM_tot*x_CM_tot #roll-yaw

M51 = M15#pitch-surge
M52 = M25#pitch-sway
M53 = M35#pitch-heave
M54 = M45#pitch-roll
M55 = I0_XX + I0_ZZ #pitch-pitch
M56 = 0#pitch-yaw

M61 = M16#yaw-surge
M62 = M26#yaw-sway
M63 = M36#yaw-heave
M64 = M46#yaw-roll
M65 = M56#yaw-pitch
M66 = I0_XX + I0_YY#yaw-yaw


M = np.array([  [M11, M12, M13, M14, M15, M16],
                [M21, M22, M23, M24, M25, M26],
                [M31, M32, M33, M34, M35, M36],
                [M41, M42, M43, M44, M45, M46],
                [M51, M52, M53, M54, M55, M56],
                [M61, M62, M63, M64, M65, M66]])



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