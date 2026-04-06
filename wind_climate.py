import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({
    "text.usetex": False,
    "lines.linewidth": 1.2,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 18,
    "axes.labelsize": 20,
    "axes.titlesize": 22,
    "legend.fontsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "lines.linewidth": 2.5,
    "figure.dpi": 300
})

# ----------------------------
# WIND speed
# ----------------------------

df = pd.read_csv("wind\windSpeed.csv") # this files are at 10 m
percentile = df["perc"]        # percentiles
wind_speed = df["val"]         # wind speed [m/s]
mean_speed = wind_speed.mean()
print("Mean wind speed (aprox):", mean_speed) # Calculate the mean wind spee

# ----------------------------
# Wind Direction
# ----------------------------
df = pd.read_csv("wind\windSpeedRose.csv")
directions_deg = df["center_degrees"]
prob = df["value"]  # 
prob = prob / prob.sum()
theta = np.deg2rad(directions_deg)
theta = np.append(theta, theta[0])
prob = np.append(prob, prob[0])
mean_direction_deg = np.arctan2(np.sum(prob * np.sin(theta)), np.sum(prob * np.cos(theta)))
mean_direction_deg = np.rad2deg(mean_direction_deg) % 360
print("Mean wind direction (aprox):", mean_direction_deg)

# ----------------------------
# Wind Power Density
# ----------------------------

df = pd.read_csv("wind\powerDensity.csv")
percentile = df["perc"]
power_density = df["val"]


# ----------------------------
# Power 
# ----------------------------
rho= 1.225
D= 1.2
A= np.pi * (D/2)**2
Cp= 0.2

# ----------------------------
# Log Law? 
# ----------------------------
z_ref = 10 # reference height (10 m)
z = 1.2 # height of the turbine hub 
z0 = 0.003  # roughness length (for open terrain)
V12 = wind_speed * (np.log(z / z0) / np.log(z_ref / z0)) 
mean_v12 = V12.mean()
print("Mean wind speed (aprox):", mean_v12) # Calculate the mean wind speed (approximate)

# Wind Power
P_wind = 0.5 * rho * A * (V12**3)
# Power extracted by the turbine
P_turbine = Cp * P_wind
# ----------------------------
# Add to DataFrame
# ----------------------------
df["V_1.2m"] = V12
df["P_wind_W"] = P_wind
df["P_turbine_W"] = P_turbine

print(df.head())

# ----------------------------
# Mean Power
# ----------------------------
print("\nMean Wind Power (W):", np.mean(P_wind))
print("Mean Turbine Power (W):", np.mean(P_turbine))
# ----------------------------
# AEP
# ----------------------------
p = df["perc"].values / 100
pdf = np.diff(np.insert(p, 0, 0))
pdf = pdf / np.sum(pdf)
P = df["P_turbine_W"].values
AEP_Wh = np.sum(P * pdf * 8760)
AEP_kWh = AEP_Wh / 1000
print("AEP (Wh/year):", AEP_Wh)
print("AEP (kWh/year):", AEP_kWh)

# ----------------------------
# Turbulence intensity
# ----------------------------
iref = 0.14
b= 5.6
sigma= iref *(0.75* V12 + 5.6)
TI= sigma/V12
df["TI"] = TI
mean_TI= TI.mean()
print("Mean Turbulence Intensity:", mean_TI)

# ----------------------------
# PLOT WIND speed
# ----------------------------

plt.figure(figsize=(10,5))
plt.plot(percentile, wind_speed, marker='o')
plt.xlabel("Percentile (%)")
plt.ylabel("Wind speed (m/s)")
plt.grid(True)
plt.savefig('wind_speed.pdf', format='pdf', bbox_inches='tight')

# ----------------------------
# PLOT WIND ROSE
# ----------------------------
plt.figure(figsize=(7,7))
ax = plt.subplot(111, polar=True)
ax.plot(theta, prob, marker='o')
ax.fill(theta, prob, alpha=0.3)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
plt.title("Wind Rose ")
plt.savefig('wind_direction.pdf', format='pdf', bbox_inches='tight')


# ----------------------------
# PLOT Power Density 
# ----------------------------
plt.figure(figsize=(10,5))
plt.plot(percentile, power_density, marker='o')

plt.title("Wind Power Density vs Percentile (10 m)")
plt.xlabel("Percentile (%)")
plt.ylabel("Wind Power Density (W/m²)")
plt.grid(True)
plt.savefig('power_density.pdf', format='pdf', bbox_inches='tight')

plt.grid(True)
# ----------------------------
# wind at V1.2 m
# ----------------------------
plt.figure(figsize=(10,5))
plt.plot(df["perc"], df["V_1.2m"], marker='o')

plt.title("Wind Speed at 1.2 m (corrected)")
plt.xlabel("Percentile (%)")
plt.ylabel("Wind speed (m/s)")
plt.grid(True)
# ----------------------------
# Wind Poer
# ----------------------------
plt.figure(figsize=(10,5))
plt.plot(df["V_1.2m"], df["P_wind_W"], marker='o')

plt.title("Wind Power Available")
plt.xlabel("Wind speed (m/s)")
plt.ylabel("Power (W)")
plt.grid(True)

# ----------------------------
# Turbine Power
# ----------------------------
plt.figure(figsize=(10,5))
plt.plot(df["V_1.2m"], df["P_turbine_W"], marker='o')

plt.title("Turbine Power Output (Cp model)")
plt.xlabel("Wind speed (m/s)")
plt.ylabel("Power (W)")
plt.grid(True)
plt.show()




