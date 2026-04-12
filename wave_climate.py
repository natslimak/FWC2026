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

# Filter data just in summer
ds['time'] = pd.to_datetime(ds['time'].values)

# Filter summer months
summer_ds = ds.where(ds['time.month'].isin([6, 7, 8, 9]), drop=True)


TIME = pd.DatetimeIndex(ds.coords['time'])
TIME_SUMMER = pd.DatetimeIndex(summer_ds.coords['time'])

VHM0 = np.array(summer_ds["VHM0"][:, 0, 0])  # sea surface wave significant height
VHM0_WW = np.array(summer_ds["VHM0_WW"][:, 0, 0])  # sea surface wind wave significant height
VMDR = np.array(summer_ds["VMDR"][:, 0, 0]) # Sea surface wave from direction
VMDR_WW = np.array(summer_ds["VMDR_WW"][:, 0, 0]) # Sea surface wind wave from direction
VPED = np.array(summer_ds["VPED"][:, 0, 0]) # Sea surface wave from direction at variance spectral density maximum
VTM01_WW = np.array(summer_ds["VTM01_WW"][:, 0, 0]) # Sea surface wind wave mean period
VTM02 = np.array(summer_ds["VTM02"][:, 0, 0]) # Sea surface wave mean period from variance spectral density second frequency moment
VTM10 = np.array(summer_ds["VTM10"][:, 0, 0]) # Sea surface wave mean period from variance spectral density inverse frequency moment
VTPK = np.array(summer_ds["VTPK"][:, 0, 0]) # Sea surface wave period at variance spectral density maximum

#%%  CLEANING DATA

print("NaNs in VTPK:", np.isnan(VTPK).sum())
print("Zeros in VTPK:", np.sum(VTPK == 0))
print("Min Tp:", np.nanmin(VTPK))
print("Max Tp:", np.nanmax(VTPK))

valid_mask = (~np.isnan(VTPK)) 

VTPK_clean = VTPK[valid_mask]
VHM0_clean = VHM0[valid_mask]
TIME_clean = TIME_SUMMER[valid_mask]

#%%  PEAK FREQUENCY ANALYSIS
Hs_smooth = pd.Series(VHM0_clean).rolling(24*30).mean()
fp = 1 / VTPK_clean  # Peak frequency in Hz
T_smooth = pd.Series(VTPK_clean).rolling(24*30).mean()
fp_smooth = pd.Series(fp).rolling(24*30).mean()

# Plot of unfiltered significant wave height over time (1982-2026)
fig, ax = plt.subplots(1, 1)
ax.plot(TIME_clean, Hs_smooth, linewidth=2, color='r', alpha=1, zorder=3,
        label='Monthly smoothed')
ax.plot(TIME_clean, VHM0_clean, linewidth=1, color='b', alpha=0.7, zorder=2,
        label='Unfiltered')
ax.set_xlabel("Time [year]")
ax.set_ylabel("Significant height [m]")
ax.grid(which='major', alpha=0.5, zorder=1)
ax.grid(which='minor', alpha=0.2, zorder=1)
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

# distribution of significant wave height
fig, ax = plt.subplots(1, 1)
bins = np.linspace(min(VHM0_clean), max(VHM0_clean), 50)
ax.hist(VHM0_clean, bins=bins, density=False, alpha=0.7, color='b',
        align='mid', linewidth=0.5, edgecolor="black", zorder=2,
        label='Unfiltered')
ax.hist(Hs_smooth, bins=bins, density=False, alpha=0.7, color='r',
        align='mid', linewidth=0.5, edgecolor="black", zorder=3,
        label='Monthly smoothed')
ax.set_ylabel("N° of observations [-]")
ax.set_xlabel("Significant height [m]")
ax.grid(which='major', alpha=0.5, zorder=1)
ax.grid(which='minor', alpha=0.2, zorder=1)
ax.set_xlim([-0.1, 4])
ax.set_ylim([0, None])
ax.legend(loc='upper right', frameon=True)
ax.minorticks_on()
ax.tick_params(direction='in', which='major', length=10,
               right=True, top =True, left=True, bottom=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True,
               labelright=False)
ax.tick_params(direction='in', which='minor', length=5, bottom=True,
               top=True, left=True, right=True)

# Plot of unfiltered peak frequency over time (1982-2026)
fig, ax = plt.subplots(1, 1)
ax.plot(TIME_clean, fp, linewidth=1, color='b', alpha=0.7, zorder=2,
        label='Unfiltered')
ax.plot(TIME_clean, fp_smooth, linewidth=2, color='r', alpha=1, zorder=3,
        label='Monthly smoothed')
ax.set_xlabel("Time [year]")
ax.set_ylabel("Peak Frequency [Hz]")
ax.set_ylim([0, 0.7])
ax.set_xlim([TIME_clean[0], TIME_clean[-1]])
ax.legend(loc='upper right', frameon=True)
ax.grid(which='major', alpha=0.5, zorder=1)
ax.grid(which='minor', alpha=0.2, zorder=1)
ax.minorticks_on()
ax.tick_params(direction='in', which='major', length=10,
               right=True, top =True, left=True, bottom=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True,
               labelright=False)
ax.tick_params(direction='in', which='minor', length=5, bottom=True,
               top=True, left=True, right=True)

# distribution of peak period
fig, ax = plt.subplots(1, 1)
bins = np.linspace(min(VTPK_clean), max(VTPK_clean), 50)
ax.hist(VTPK_clean, bins=bins, density=False, alpha=0.7, color='b',
        align='mid', linewidth=0.5, edgecolor="black", zorder=2,
        label='Unfiltered')
ax.hist(T_smooth, bins=bins, density=False, alpha=0.7, color='r',
        align='mid', linewidth=0.5, edgecolor="black", zorder=3,
        label='Monthly smoothed')
ax.set_ylabel("N° of observations [-]")
ax.set_xlabel("Wave peak period [s]")
ax.grid(which='major', alpha=0.5, zorder=1)
ax.grid(which='minor', alpha=0.2, zorder=1)
# ax.set_xlim([-0.1, 4])
ax.set_ylim([0, None])
ax.legend(loc='upper right', frameon=True)
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

#%% JOINT DISTRIBUTION


Hs = VHM0_clean   # significant wave height
Tp = VTPK_clean   # peak period


Hs_bins = np.arange(0, np.max(Hs) + 0.5, 0.5)
Tp_bins = np.arange(0, np.max(Tp) + 1.0, 1.0)

H, Hs_edges, Tp_edges = np.histogram2d(Hs, Tp, bins=[Hs_bins, Tp_bins])

P = H / np.sum(H)

Hs_centers = 0.5 * (Hs_edges[:-1] + Hs_edges[1:])
Tp_centers = 0.5 * (Tp_edges[:-1] + Tp_edges[1:])

f = np.linspace(0.01, 1.0, 1000)

spectra = []
weights = []
max_peak_f = []

start_Tp_idx = 3   # skip very low Tp values

for i in range(len(Hs_centers[:1])):  # currently only lowest Hs bin
    for j in range(start_Tp_idx, min(5, len(Tp_centers))):

        prob = P[i, j]

        # Skip non-occurring sea states
        if prob == 0:
            continue

        # Representative values
        Hs_rep = Hs_centers[i]
        Tp_rep = Tp_centers[j]

        # Compute spectrum
        S = jonswap_spectrum(f, Hs_rep, Tp_rep)

        # Store results
        spectra.append(S)
        weights.append(prob)

        # Store peak frequency
        peak_idx = np.argmax(S)
        max_peak_f.append(f[peak_idx])

# weighted mean spectrum
spectra = np.array(spectra)
weights = np.array(weights)
S_mean = np.sum(spectra.T * weights, axis=1)

# Joint distribution
fig, ax = plt.subplots(1, 1)
pcm = ax.pcolormesh(Hs_edges, Tp_edges, P.T, cmap='plasma')
ax.set_xlabel("Hs [m]")
ax.set_ylabel("Tp [s]")
ax.grid(which='major', alpha=0.5, zorder=1)
ax.set_xticks(Hs_edges)  # every 2 bins
ax.set_yticks(Tp_edges[::2])
ax.set_title("Joint Distribution")
fig.colorbar(pcm, label="Probability", ax=ax)
ax.minorticks_on()
ax.tick_params(direction='in', which='major', length=10,
               right=True, top =True, left=True, bottom=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True,
               labelright=False)
ax.tick_params(direction='in', which='minor', length=5, bottom=True,
               top=True, left=True, right=True)

# Example spectra
fig, ax = plt.subplots(1, 1)
ax.plot(f, S_mean, 'b', linewidth=2, alpha=1, zorder=2)
for max_f in max_peak_f:
    ax.axvline(max_f, color='r', linewidth=2, alpha=1, zorder=3, linestyle='--',
               label=f'Peak frequency {max_f:.2} Hz')
ax.set_ylabel("S(f)")
ax.set_xlabel("Frequency [Hz]")
ax.set_title('Weighted mean JONSWAP spectrum')
ax.legend(loc='upper right', frameon=False)
# ax.set_ylim([0, 0.7])
ax.set_xlim([f[0], f[-1]])
ax.grid(which='major', alpha=0.5, zorder=1)
ax.grid(which='minor', alpha=0.2, zorder=1)
ax.minorticks_on()
ax.tick_params(direction='in', which='major', length=10,
               right=True, top =True, left=True, bottom=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True,
               labelright=False)
ax.tick_params(direction='in', which='minor', length=5, bottom=True,
               top=True, left=True, right=True)
#%% FREE SURFACE SIMULATION

# Time array
dt = 0.5                         # time step [s]
T_total = 600                    # total duration [s]
t = np.arange(0, T_total, dt)

# Frequency spacing
df = f[1] - f[0]

# Random phases
phi = 2 * np.pi * np.random.rand(len(f))

# Amplitudes from spectrum
A = np.sqrt(2 * S_mean * df)

# Build free surface elevation
eta = np.zeros_like(t)

for i in range(len(f)):
    eta += A[i] * np.cos(2 * np.pi * f[i] * t + phi[i])


fig, ax = plt.subplots(1, 1)
ax.plot(t, eta, linewidth=2, color='b', alpha=1, zorder=3)
ax.set_xlabel("Time [s]")
ax.set_ylabel("Surface elevation [m]")
ax.set_title("Simulated Free Surface Elevation")
ax.grid(which='major', alpha=0.5, zorder=1)
ax.grid(which='minor', alpha=0.2, zorder=1)
ax.set_ylim([-0.3, 0.3])
ax.set_xlim([t[0], t[-1]])
ax.minorticks_on()
ax.tick_params(direction='in', which='major', length=10,
               right=True, top =True, left=True, bottom=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True,
               labelright=False)
ax.tick_params(direction='in', which='minor', length=5, bottom=True,
               top=True, left=True, right=True)

plt.show()



