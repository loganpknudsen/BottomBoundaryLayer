"""Figure 16: time-integrated TKE budget terms for Sim7 (k, ASP = AGSP + WSP, GSP, BP, -epsilon),
integrated from T = 5; symmetric-log scale with linear threshold 1e-15.

Source: Result_Calculations/TKE_terms.ipynb, cells 3, 11 (time crop), 15-19, 28-34 and 37
(origin/Base_Case numbering; cell 33 in the submission commit 27f8414). Same code as Fig 14.
Run: Sim7 (S = 1, Ri_inv = 0.75, delta = 0.8), March output TKE_terms_..._Ri_inv_0.75_delta_0.8.nc.
"""
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import SIM_DIR, save  # noqa: E402

# TODO(decision): which notebook version the preprint figure follows.
#   "base_case"  (origin/Base_Case cell 37): adds a grey "Residual" curve, no explicit ticks,
#                x-limits fixed at [4.95, 20.01], which were set for Sim9; no x-limits are recorded
#                for Sim7, so this version cannot reproduce Fig 16's axis without a further choice.
#                It also integrates from T = 0 (time crop commented out); the crop to T = 5
#                (ts = 1000) is kept here as instructed.
#   "submission" (27f8414 cell 33): no residual curve, x-ticks every 2.5 and y-ticks +-1e-6,
#                +-1e-10, +-1e-14, x-limits from the cropped record, which was
#                ocean_data.isel(time=slice(ts, tf + 1000)). The preprint axis (5 to 27.5) implies
#                ts = 1000, tf = 4500; tf is inferred from the figure, not recorded in the code.
# The preprint Fig 16 has no residual curve and those exact ticks, i.e. it matches "submission".
VERSION = None
if VERSION not in ("base_case", "submission"):
    raise NotImplementedError('Fig 16: follow the "base_case" or the "submission" notebook version?')

file = "TKE_terms_Sinf_1.0_Ri_inv_0.75_delta_0.8.nc"
ocean_data = xr.open_dataset(SIM_DIR / file)

# Run parameters (Sim7), cell 3
theta = 0.1
f = 10 ** (-4)
S = 1.0
fstar = f * np.cos(theta) * (1 + S**2) ** (0.5)

# Time crop (cell 11): integration starts at output ts (T = 5)
ts = 1000
tf = 4500
ocean_data = ocean_data.isel(time=slice(ts, tf + 1000)) if VERSION == "submission" else ocean_data.isel(time=slice(ts, None))


def with_coords(da):
    seconds = da.time / pd.Timedelta("1s")
    return da.assign_coords(intp=('time', (seconds * fstar / (2 * np.pi)).data), sec=('time', seconds.data))


AGSP = with_coords(ocean_data.AGSP)
WSP = with_coords(ocean_data.WSP)
GSP = with_coords(ocean_data.GSP)
BFLUX = with_coords(ocean_data.BFLUX)
E = with_coords(ocean_data.E)
integrated_AGSP = AGSP.isel(time=slice(0, -1)).cumulative_integrate(coord="sec")
integrated_WSP = WSP.isel(time=slice(0, -1)).cumulative_integrate(coord="sec")
integrated_GSP = GSP.isel(time=slice(0, -1)).cumulative_integrate(coord="sec")
integrated_BFLUX = BFLUX.isel(time=slice(0, -1)).cumulative_integrate(coord="sec")
integrated_E = (-1 * E.isel(time=slice(0, -1))).cumulative_integrate(coord="sec")
KE_turb_ip = with_coords(ocean_data.k)

start = 0
end = -1
fig = plt.figure(figsize=(7, 5))
plt.plot(KE_turb_ip.intp[start:end], KE_turb_ip[start:end] - KE_turb_ip[start], label="$k$", color="k")
plt.plot(integrated_AGSP.intp[start:end], (integrated_AGSP[start:end] + integrated_WSP[start:end]) - (integrated_AGSP[start] + integrated_WSP[start]), label=r"$\int_0^tASP\,dt'$", color="red")
plt.plot(integrated_GSP.intp[start:end], integrated_GSP[start:end] - integrated_GSP[start], label=r"$\int_0^tGSP\,dt'$", color="green")
plt.plot(integrated_BFLUX.intp[start:end], integrated_BFLUX[start:end] - integrated_BFLUX[start], label=r"$\int_0^tBP\,dt'$", color="blue")
plt.plot(integrated_E.intp[start:end], integrated_E[start:end] - integrated_E[start], label=r"$-\int_0^t\epsilon\,dt'$", color="pink")
if VERSION == "base_case":
    resid = KE_turb_ip - integrated_AGSP - integrated_GSP - integrated_WSP - integrated_BFLUX - integrated_E
    plt.plot(resid.intp[start:end], resid[start:end] - resid[start], label="Residual", color="grey", alpha=0.9)
plt.xlabel("Modified Inertial Period")
plt.ylabel(r"k $(m^2 s^{-2})$")
plt.yscale("symlog", linthresh=1e-15)
if VERSION == "submission":
    plt.xticks([5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30])
    plt.yticks([1e-6, -1e-6, 1e-10, -1e-10, 1e-14, -1e-14])
plt.legend(loc="lower left", prop={'size': 7.5}, ncol=2, handlelength=0.5)
if VERSION == "submission":
    plt.xlim([KE_turb_ip.intp[start] - 0.05, KE_turb_ip.intp[end] + 0.01])
else:
    # TODO(decision): x-limits for Sim7 are not recorded in the Base_Case notebook.
    raise NotImplementedError("Fig 16 (base_case version): x-limits for Sim7 are not recorded in the notebook")

save(fig, __file__)
