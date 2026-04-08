import xarray as xr
import numpy as np
import pandas as pd

ds = xr.open_dataset("wave/wave_climate_1982_2026.nc")

TIME = pd.DatetimeIndex(ds.coords['time'])

VHM0 = np.array(ds["VHM0"][:, 0, 0])  # sea surface wave significant height
VHM0_WW = np.array(ds["VHM0_WW"][:, 0, 0])  # sea surface wind wave significant height
VMDR = np.array(ds["VMDR"][:, 0, 0]) # Sea surface wave from direction
VMDR_WW = np.array(ds["VHM0_WW"][:, 0, 0]) # Sea surface wind wave from direction
VPED = np.array(ds["VMDR_WW"][:, 0, 0]) # Sea surface wave from direction at variance spectral density maximum
VTM01_WW = np.array(ds["VTM01_WW"][:, 0, 0]) # Sea surface wind wave mean period
VTM02 = np.array(ds["VTM02"][:, 0, 0]) # Sea surface wave mean period from variance spectral density second frequency moment
VTM10 = np.array(ds["VTM10"][:, 0, 0]) # Sea surface wave mean period from variance spectral density inverse frequency moment
VTPK = np.array(ds["VTPK"][:, 0, 0]) # Sea surface wave period at variance spectral density maximum