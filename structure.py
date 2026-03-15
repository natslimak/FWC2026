'''This code is for a fast dimensioning of a semi-sub triangular floater
to better understand the imput mesurement you can refer to the README file'''

# %%

#imports
import numpy as np
from scipy.linalg import eig
import sys

## ambient parameters

g = 9.81            # Gravity [m/s^2]
h = 320             # Water depth [m]
rhow = 1025         # Water density [kg/m^3]
rhoa = 1.29         # Air density [kg/m^3]

# tower and rotor parameter
mt = 5        # Tower mass [kg]
z_hub = 1.2            # tower hight [m]
d_t = 0.1       #tower diameter

Drot = 126          # Rotor diameter [m]
Vrated = 11.4       # Rated wind speed [m/s]
Trated = 900      # Thrust at rated [N]

mn = 6.5          # Nacelle mass [kg]
mr = 1          # Rotor mass [kg]

# floater parameters

m_bf = 1            # Back Floater mass [kg]
m_ff = 1            # Front Floater mass [kg]
h_bf = 1            # Back floater hight [m]
h_ff = 1            # Front floater hight [m] 
D_b = 0.2           # Back floater Diameter [m]
D_f = 0.2           # Front floater Diameter [m]

m_bb = 0             # Back Ballast mass [kg]
rho_ballast = 1025         # Ballast density [kg/m^3]
#h_bb = 0.5             # Back Ballast hight [m]
#h_fb = 0.5             # Front Ballast hight [m] 


d_b = 1             # Distance between back floater [m]
d_f = 1             # distance between the center of the back row and front floater [m]


# mooring parameter 

Kmoor = 12       # Mooring stiffness [N/m]
zmoor = -70         # Mooring fairlead [m]

#others

#Cm = 1.0            # Added mass coefficient
#Cd = 0.6            # Drag coefficient


####################### CALCULATION #################################################
# 
# Calculate the draft 

m_back = 2* (mt + mn + mr + m_bb + m_bf)   # Back side mass [kg] 
 
A_back_floater = 2 * np.pi * D_b**2 /4 # Total area back floater [m^2]

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

mtot = m_back + m_front

print('Total weigth of the structure is', mtot, 'kg')

#calculate the hight of the tower and the ballast

h_t = z_hub-(h_bf-draft) # Tower hight [m]
h_bb = m_bb/(rho_ballast* np.pi * D_b**2 /4) # back balast hight [m]
h_fb = m_fb/(rho_ballast* np.pi * D_f**2 /4) # back balast hight [m]

if  h_fb >= draft: 

    sys.exit('We need an heavier weigth')

print ('The tower shold be',h_t,'m')
# Calculate center of mass coordinates

zCMt = h_t/2 + (h_bf-draft) # z center of mass tower

zCM_bf = (h_bf/2) - draft   # z center of back floater
zCM_ff = (h_ff/2) - draft   # z center of front floater

zCM_bb = (h_bb/2) - draft   # z center of back ballast
zCM_fb = (h_fb/2) - draft   # z center of front ballast

zCMtot = ( m_bb * zCM_bb + m_fb * zCM_ff + m_bf * zCM_bf + m_ff * zCM_ff + mt * zCMt + (mn + mr) * z_hub  ) / mtot

xCMtot = (m_front * d_f ) / mtot


#Inertia moment

I0_XX_turb = (mn+mr)*(d_f/3)**2
I0_YY_turb = (mn+mr)*(d_b/2)**2
I0_ZZ_turb = (mn+mr)*z_hub**2
I0_XY_turb = (mn+mr)*d_f/3 * d_b/2

I0_XX_tower = mt*(3*(d_t/2)**2+ h_t**2)/12 + mt * (d_f/3)**2
I0_YY_tower = mt*(3*(d_t/2)**2+ h_t**2)/12 + mt * (d_b/2)**2
I0_ZZ_tower = 0.5 *mt*(d_t/2)**2+ mt * zCMt**2 
I0_XY_tower = mt*d_f/3 * d_b/2

I0_XX_backfloat =  m_bf*(3*(D_b/2)**2+ h_bf**2)/12 + m_bf * (d_f/3)**2
I0_YY_backfloat =  m_bf*(3*(D_b/2)**2+ h_bf**2)/12 + m_bf * (d_b/2)**2
I0_ZZ_backfloat = 0.5 *m_bf*(D_b/2)**2+ m_bf * zCM_bf**2 
I0_XY_backfloat = m_bf*d_f/3 * d_b/2

I0_XX_backBallast =  m_bb*(3*(D_b/2)**2+ h_bb**2)/12 + m_bb * (d_f/3)**2
I0_YY_backBallast =  m_bb*(3*(D_b/2)**2+ h_bb**2)/12 + m_bb * (d_b/2)**2
I0_ZZ_backBallast = 0.5 *m_bb*(D_b/2)**2+ m_bb * zCM_bb**2 
I0_XY_backBallast = m_bb*d_f/3 * d_b/2

I0_XX_frontfloat =  m_ff*(3*(D_f/2)**2+ h_ff**2)/12 + m_ff * (d_f*2/3)**2
I0_YY_frontfloat =  m_ff*(3*(D_f/2)**2+ h_ff**2)/12 
I0_ZZ_frontfloat = 0.5 *m_ff*(D_f/2)**2+ m_ff * zCM_ff**2 
I0_XY_frontfloat = m_ff * d_f*2/3 * d_b/2

I0_XX_frontBallast =  m_fb*(3*(D_f/2)**2+ h_fb**2)/12 + m_fb * (d_f*2/3)**2
I0_YY_frontBallast =  m_fb*(3*(D_f/2)**2+ h_fb**2)/12 
I0_ZZ_frontBallast = 0.5 *m_fb*(D_f/2)**2+ m_fb * zCM_fb**2 
I0_XY_frontBallast = m_ff * d_f*2/3 * d_b/2

I0_XX = I0_XX_turb + I0_XX_backBallast + I0_XX_backfloat +I0_XX_frontfloat + I0_XX_frontBallast + I0_XX_tower
I0_YY = I0_YY_turb + I0_YY_backBallast + I0_YY_backfloat +I0_YY_frontfloat + I0_YY_frontBallast + I0_YY_tower
I0_ZZ = I0_ZZ_turb + I0_ZZ_backBallast + I0_ZZ_backfloat +I0_ZZ_frontfloat + I0_ZZ_frontBallast + I0_ZZ_tower
I0_XY = I0_XY_turb + I0_XY_backBallast + I0_XY_backfloat +I0_XY_frontfloat + I0_XY_frontBallast + I0_XY_tower
I0_YX = I0_XY

# build mass matrix

M11 = mtot
M12 = 0
M13 = 0
M14 = mtot * zCMtot

M21 = 0
M22 = mtot
M23 = (mt + m_bf + m_bb + mn+ mr) * d_b
M24 = -2*(mt + m_bf + m_bb + mn + mr ) * d_f/3 - (m_ff + m_fb) * 2*d_f/3

M31 = 0
M32 = M23 
M33 = I0_YY + I0_ZZ
M34 = -I0_YX

M41 = M14
M42 = M24
M43 = M34
M44 =I0_XX + I0_ZZ

M = np.array ([[M11, M12, M13, M14],
               [M21, M22, M23, M24],
               [M31, M32, M33, M34],
               [M41,M42,M43, M44]])

# Added mass matrix


A = np.array([[0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0]])

# restoring matrix
#C = np.array([[0. ,0.],[0. , m_total*g*(zCM_tot - zCB) + rhow*g*IAA]])

# Damping matrix
#B = np.array([[B11, 0],[0 ,0]])

print("M: ",M)
print("A: ",A)
#print("C: ",C)
'''
# hydrodynamic stiffness
zCB = - draft/2
IAA = np.pi * (D_f /2)**4/4 
Chst = np.array([[0. ,0.],[0. , mtot*g*(zCB-zCMtot) + rhow*g*IAA]])



# Mooring restoring matrix
Cmoor = np.array([[Kmoor, Kmoor*zmoor], [Kmoor*zmoor, Kmoor*(zmoor**2)]])

C = Chst + Cmoor
print("C: ",C)

# angle for rated trust 

tetha = np.rad2deg ((Trated*z_hub)/Chst[1,1])

## ----  NATURAL FREQUENCIES ----

MA = np.add(A,M)

eigVal, eigVec = eig(C,MA)

eigVal = np.real(eigVal)


# Natural frequencies
omeganat = np.sqrt(eigVal) 
fnat = omeganat/2/np.pi

print("Natural frequencies (rad/s): ", omeganat)

#SparBuoyData["fnat"] = fnat

# Natural periods
Tnat = 1./fnat

# Display surge and pitch natural periods
print(f'Surge period: {Tnat[0]:.2f} [s]')
print(f'Pitch period: {Tnat[1]:.2f} [s]')

'''