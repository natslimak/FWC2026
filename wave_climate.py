# ----------------------
# IMPORTS
# ----------------------
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------
# FUNCTIONS
# -----------------------

# JONSWAP spectrum function
def jonswap_spectrum(f, Hs, Tp, gamma=3.3):
    fp = 1.0 / Tp

    # Sigma definition
    sigma = np.where(f <= fp, 0.07, 0.09)

    # Core spectrum (PM-like part)
    S_pm = 0.3125 * Hs**2 * Tp * (f / fp)**(-5) * \
           np.exp(-1.25 * (f / fp)**(-4))

    # Peak enhancement term
    exponent = np.exp(-0.5 * ((f / fp - 1) / sigma)**2)
    peak_enhancement = gamma**exponent

    # Final JONSWAP spectrum
    S = S_pm * peak_enhancement

    return S

# ----------------------
# INPUT DATA
# ----------------------

ds = xr.open_dataset("wave/wave_climate_1982_2026.nc")

TIME = pd.DatetimeIndex(ds.coords['time'])

VHM0 = np.array(ds["VHM0"][:, 0, 0])  # sea surface wave significant height
VHM0_WW = np.array(ds["VHM0_WW"][:, 0, 0])  # sea surface wind wave significant height
VMDR = np.array(ds["VMDR"][:, 0, 0]) # Sea surface wave from direction
VMDR_WW = np.array(ds["VMDR_WW"][:, 0, 0]) # Sea surface wind wave from direction
VPED = np.array(ds["VPED"][:, 0, 0]) # Sea surface wave from direction at variance spectral density maximum
VTM01_WW = np.array(ds["VTM01_WW"][:, 0, 0]) # Sea surface wind wave mean period
VTM02 = np.array(ds["VTM02"][:, 0, 0]) # Sea surface wave mean period from variance spectral density second frequency moment
VTM10 = np.array(ds["VTM10"][:, 0, 0]) # Sea surface wave mean period from variance spectral density inverse frequency moment
VTPK = np.array(ds["VTPK"][:, 0, 0]) # Sea surface wave period at variance spectral density maximum

print("TIME:", TIME)

# -----------------------
# CLEANING DATA
# -----------------------

print("NaNs in VTPK:", np.isnan(VTPK).sum())
print("Zeros in VTPK:", np.sum(VTPK == 0))
print("Min Tp:", np.nanmin(VTPK))
print("Max Tp:", np.nanmax(VTPK))

valid_mask = (~np.isnan(VTPK)) 

VTPK_clean = VTPK[valid_mask]
VHM0_clean = VHM0[valid_mask]
TIME_clean = TIME[valid_mask]

# ------------------------
# PEAK FREQUENCY ANALYSIS
# ------------------------

fp = 1 / VTPK_clean  # Peak frequency in Hz

# Plot of unfiltered peak frequency over time (1982-2026)
plt.figure()
plt.plot(TIME_clean, fp)
plt.xlabel("Time")
plt.ylabel("Peak Frequency (Hz)")
plt.title("Peak Wave Frequency over Time")
plt.grid()

# Averaging over ~1 month to see trends more clearly
fp_smooth = pd.Series(fp).rolling(24*30).mean()  # ~monthly smoothing

# Plot of monthly smoothed peak frequency over time
plt.figure(figsize=(10,4))
plt.plot(TIME_clean, fp_smooth)
plt.xlabel("Time")
plt.ylabel("Peak Frequency (Hz)")
plt.title("Smoothed Peak Wave Frequency")
plt.grid()

# Statistical analysis of peak frequency
fp_mean = np.mean(fp) # Typical sea state
fp_p90 = np.percentile(fp, 90) # More engergetic / extreme sea state
weights = VHM0_clean**2  # wave energy ∝ Hs²; waves that contribute more to the energy should have more weight
fp_weighted = np.sum(fp * weights) / np.sum(weights)

print(f"Mean Peak Frequency: {fp_mean:.3f} Hz")
print(f"90th Percentile Peak Frequency: {fp_p90:.3f} Hz")
print(f"Energy-Weighted Mean Peak Frequency: {fp_weighted:.3f} Hz")

# ------------------------
# SPECTRUm AND FREE SURFACE SIMULATION (EXAMPLE - Only with first valid data point)
# ------------------------

# Time for simulation (seconds)
t = np.linspace(0, 600, 2000)

# Frequency domain
f = np.linspace(0.02, 1, 200)
df = f[1] - f[0]

Hs = VHM0_clean[0]
Tp = VTPK_clean[0]

# Spectrum
S = jonswap_spectrum(f, Hs, Tp)

# Convert to amplitudes
a = np.sqrt(2 * S * df)

# Random phases
phi = np.random.uniform(0, 2*np.pi, len(f))

# Angular frequency
omega = 2 * np.pi * f

# Build surface elevation
eta = np.zeros_like(t)

for j in range(len(f)):
    eta += a[j] * np.cos(omega[j]*t + phi[j])

plt.figure()
plt.plot(f, S)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Energy Density")
plt.title("Wave Spectrum (JONSWAP)")
plt.grid()

plt.figure()
plt.plot(t, eta)
plt.xlabel("Time (s)")
plt.ylabel("Surface Elevation (m)")
plt.title("Simulated Free Surface Elevation")
plt.grid()
plt.show()

