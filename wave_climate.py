"""This file analyzes the wave climate at the competition's site."""
#%% IMPORTS
from pathlib import Path

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

#%% Plot commands
# Size
mpl.rcParams['figure.figsize'] = (16, 10)

# Font size of label, title, and legend
mpl.rcParams['font.size'] = 25
mpl.rcParams['xtick.labelsize'] = 25
mpl.rcParams['ytick.labelsize'] = 25
mpl.rcParams['axes.labelsize'] = 25
mpl.rcParams['axes.titlesize'] = 25
mpl.rcParams['legend.fontsize'] = 25

# Lines and markers
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['scatter.marker'] = 'd'
plt_marker = 'd'

# Latex font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'

# Export
mpl.rcParams['savefig.bbox'] = "tight"
#%%  FUNCTIONS

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

#%%  INPUT DATA

ROOT = Path(__file__).parent
WAVE_FILE_PATH = ROOT /  'wave' / 'wave_climate_1982_2026.nc'
ds = xr.open_dataset(WAVE_FILE_PATH)

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

# The time is every 3 hours
print("TIME:", TIME)


#%%  CLEANING DATA

print("NaNs in VTPK:", np.isnan(VTPK).sum())
print("Zeros in VTPK:", np.sum(VTPK == 0))
print("Min Tp:", np.nanmin(VTPK))
print("Max Tp:", np.nanmax(VTPK))

valid_mask = (~np.isnan(VTPK)) 

VTPK_clean = VTPK[valid_mask]
VHM0_clean = VHM0[valid_mask]
TIME_clean = TIME[valid_mask]

#%%  PEAK FREQUENCY ANALYSIS
Hs_smooth = pd.Series(VHM0_clean).rolling(24*30).mean()

# Plot of unfiltered significant wave height over time (1982-2026)
fig, ax = plt.subplots(1, 1)
ax.plot(TIME_clean, Hs_smooth, linewidth=2, color='r', alpha=1, zorder=3,
        label='Monthly smoothed')
ax.plot(TIME_clean, VHM0_clean, linewidth=1, color='b', alpha=0.7, zorder=2,
        label='Unfiltered')
ax.set_xlabel("Time [year]")
ax.set_ylabel("Significant height [m]")
# ax.set_title("Significant height over Time")
ax.grid(which='major', alpha=0.4, zorder=1)
ax.set_ylim([0, 4])
ax.set_xlim([TIME_clean[0], TIME_clean[-1]])
ax.legend(loc='upper right', frameon=True)
ax.minorticks_on()
ax.tick_params(direction='in', which='major', length=10,
               right=True, top =True, left=True, bottom=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True,
               labelright=False)
ax.tick_params(direction='in', which='minor', length=5, bottom=True,
               top=True, left=True, right=True)


fp = 1 / VTPK_clean  # Peak frequency in Hz
# Averaging over ~1 month to see trends more clearly
fp_smooth = pd.Series(fp).rolling(24*30).mean()

# Plot of unfiltered peak frequency over time (1982-2026)
fig, ax = plt.subplots(1, 1)
ax.plot(TIME_clean, fp, linewidth=1, color='b', alpha=0.7, zorder=2,
        label='Unfiltered')
ax.plot(TIME_clean, fp_smooth, linewidth=2, color='r', alpha=1, zorder=3,
        label='Monthly smoothed')
ax.set_xlabel("Time [year]")
ax.set_ylabel("Peak Frequency [Hz]")
# ax.set_title("Peak Wave Frequency over Time")
ax.set_ylim([0, 0.7])
ax.set_xlim([TIME_clean[0], TIME_clean[-1]])
ax.legend(loc='upper right', frameon=True)
ax.grid(which='major', alpha=0.4, zorder=1)
ax.minorticks_on()
ax.tick_params(direction='in', which='major', length=10,
               right=True, top =True, left=True, bottom=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True,
               labelright=False)
ax.tick_params(direction='in', which='minor', length=5, bottom=True,
               top=True, left=True, right=True)

# Statistical analysis of peak frequency
fp_mean = np.mean(fp)          # Typical sea state
fp_p90 = np.percentile(fp, 90) # More engergetic / extreme sea state
weights = VHM0_clean**2        # wave energy ∝ Hs²; waves that contribute more to the energy should have more weight
fp_weighted = np.sum(fp * weights) / np.sum(weights)

print(f"Mean Peak Frequency: {fp_mean:.3f} Hz")
print(f"90th Percentile Peak Frequency: {fp_p90:.3f} Hz")
print(f"Energy-Weighted Mean Peak Frequency: {fp_weighted:.3f} Hz")

#%% SPECTRUM AND FREE SURFACE SIMULATION

# Time for simulation
t = np.linspace(0, 600, 2000) # [s]

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

fig, ax = plt.subplots(1, 1)
ax.plot(f, S, linewidth=2, color='b', alpha=1, zorder=3)
ax.set_xlabel("Frequency [Hz]")
ax.set_ylabel("Energy density")
ax.set_title("JONSWAP Wave Spectrum")
ax.grid(which='major', alpha=0.5, zorder=1)
ax.grid(which='minor', alpha=0.2, zorder=1)
ax.set_ylim([0, 3])
ax.set_xlim([f[0], f[-1]])
ax.minorticks_on()
ax.tick_params(direction='in', which='major', length=10,
               right=True, top =True, left=True, bottom=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True,
               labelright=False)
ax.tick_params(direction='in', which='minor', length=5, bottom=True,
               top=True, left=True, right=True)

fig, ax = plt.subplots(1, 1)
ax.plot(t, eta, linewidth=2, color='b', alpha=1, zorder=3)
ax.set_xlabel("Time [s]")
ax.set_ylabel("Surface elevation [m]")
ax.set_title("Simulated Free Surface Elevation")
ax.grid(which='major', alpha=0.5, zorder=1)
ax.grid(which='minor', alpha=0.2, zorder=1)
ax.set_ylim([-1.5, 1.5])
ax.set_xlim([t[0], t[-1]])
ax.minorticks_on()
ax.tick_params(direction='in', which='major', length=10,
               right=True, top =True, left=True, bottom=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True,
               labelright=False)
ax.tick_params(direction='in', which='minor', length=5, bottom=True,
               top=True, left=True, right=True)


