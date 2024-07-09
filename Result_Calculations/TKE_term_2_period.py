
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# /Users/loganknudsen/Documents/UMD_Research_Local_File_Saves/
ocean_data=xr.open_dataset("BBL_w_O_updated_diagnostics_TKE_terms.nc")
KE_spatial_average_time_series = ocean_data.k.mean(dim=['xF','yC','zC'])
AGSP_spatial_average_time_series = ocean_data.AGSP.mean(dim=['xF','yC','zC'])
GSP_spatial_average_time_series = ocean_data.GSP.mean(dim=['xC','yF','zC'])
BFLUX_spatial_average_time_series = ocean_data.BFLUX.mean(dim=['xF','yC','zC'])
E_spatial_average_time_series = ocean_data.E.mean(dim=['xC','yC','zC'])
PWORK_spatial_average_time_series = ocean_data.PWORK.mean(dim=['xC','yC','zC'])
inertial_period = E_spatial_average_time_series.time/pd.Timedelta("1s")*(1e-4)/(2*np.pi)
dkdt = KE_spatial_average_time_series.differentiate(coord='sec',edge_order=2)

ts = 2081
tf = 2481
plt.plot(inertial_period[ts:tf],1e4*dkdt[ts:tf]/KE_spatial_avyerage_time_series[ts:tf],color="k",label='$\dfrac{1}{k}\dfrac{\partial k}{\partial t}$')
plt.plot(inertial_period[ts:tf],1e4*AGSP_spatial_average_time_series.values[ts:tf]/KE_spatial_average_time_series.values[ts:tf],color='red',label='AGSP/k')
plt.plot(inertial_period[ts:tf],GSP_spatial_average_time_series.values[ts:tf]/KE_spatial_average_time_series.values[ts:tf]*1e4,color='green',label='GSP/k')
plt.plot(inertial_period[ts:tf],BFLUX_spatial_average_time_series.values[ts:tf]/KE_spatial_average_time_series.values[ts:tf]*1e4,color='blue',label='BFLUX/k')
plt.plot(inertial_period[ts:tf],-1*E_spatial_average_time_series.values[ts:tf]/KE_spatial_average_time_series.values[ts:tf]*1e4,color='pink',label='$-\epsilon/k$')
plt.plot(inertial_period[ts:tf],1e4*PWORK_spatial_average_time_series.values[ts:tf]/KE_spatial_average_time_series.values[ts:tf],color='orange')
plt.axhline(color="violet")
plt.xlabel("$t\dfrac{2\pi}{f}$")
plt.ylabel("TKE growth rate/f")
plt.legend(loc='lower left')
plt.figure(figsize=(12,5))
plt.savefig("TKEterms.png")
plt.show()
