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
z_hub = 1.2            # tower lenght [m]

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
D_f = 0.1           # Front floater Diameter [m]

m_bb = 0             # Back Ballast mass [kg]

h_bb = 0.5             # Back Ballast hight [m]
h_fb = 0.5             # Front Ballast hight [m] 


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
 
A_back_floater = 2 * np.pi * D_b  # Total area back floater [m^2]

draft = m_back / (A_back_floater*rhow) #draft [m]

if draft >= h_bf: 

    sys.exit('WE ARE SINKING :(')

else:
    print ('WE ARE FLOATING :)')
    print ('The draft with', m_bb,'kg of ballast in the back is',draft,'m')

A_front_floater =  np.pi * D_f

m_front = (rhow * A_front_floater * draft) 

m_fb = m_front - m_ff

print ('We need',m_fb ,'kg of ballast in front')

mtot = m_back + m_front

print('Total weigth of the structure is', mtot, 'kg')
# Calculate center of mass coordinates

zCMt = (z_hub+h_bf-draft)/2 # z center of mass tower

zCM_bf = (h_bf/2) - draft   # z center of back floater
zCM_ff = (h_ff/2) - draft   # z center of front floater

zCM_bb = (h_bb/2) - draft   # z center of back ballast
zCM_fb = (h_fb/2) - draft   # z center of front ballast

zCMtot = ( m_bb * zCM_bb + m_fb * zCM_ff + m_bf * zCM_bf + m_ff * zCM_ff + mt * zCMt + (mn + mr) * z_hub  ) / mtot

xCMtot = (m_front * d_f ) / mtot

#print(xCMtot)

I0turb = (mn+mr)*z_hub**2
I0tower =  mt * zCMt**2
I0float =  m_ff * zCM_ff**2 + 2 * m_bf * zCM_bf**2

IOtot = I0float + I0tower + I0turb

M11 = mtot
M12 = 0
M13 = 0
M14 = 2* (mt * zCMt + m_bf * zCM_bf) 

M21 = 0
M22 = mtot
M23 = mtot * 
M24

M31
M32
M33
M34

M41
M42
M43
M44

M = np.array ([[mtot, mtot*zCMtot ],[mtot*zCMtot, IOtot]])

# Added mass matrix
z_bot= - draft
Cm = 2
A11 = - rhow * (np.pi/4) * Cm * (z_bot) * Dspar**2
A21 = A12 = -rhow*np.pi/4*Dspar*2*Cm*z_bot*2/2
A22 =  - rhow*np.pi/4*Dspar*2*Cm*z_bot*3/3

A = np.array([[0, 0],[0, 0]])

# restoring matrix
#C = np.array([[0. ,0.],[0. , m_total*g*(zCM_tot - zCB) + rhow*g*IAA]])

# Damping matrix
#B = np.array([[B11, 0],[0 ,0]])

print("M: ",M)
print("A: ",A)
#print("C: ",C)

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

