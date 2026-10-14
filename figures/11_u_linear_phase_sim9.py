"""Figure 11: perturbation u for Sim9 over the linear growth phase (T = 10, 10.5, 15, 15.5),
with isopycnals; symmetric-log colour scale, linear threshold 1e-6.

Source: Result_Calculations/velocity_visualization.ipynb (origin/Base_Case), cells 1, 5 and 23.
Run: Sim9 (S = 1, Ri_inv = 0.95, delta = 0.8), March output flow_fields_..._Ri_inv_0.95_delta_0.8.nc.
"""
import pathlib
import sys

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import SIM_DIR, save  # noqa: E402

# Run parameters (Sim9)
theta = 0.1
f = 1e-4
Sinf = 1.0
Ri_inv = 0.95
delta = 0.8
fstar = f * np.cos(theta) * (1 + Sinf**2) ** (0.5)

file = "flow_fields_Sinf_1.0_Ri_inv_0.95_delta_0.8.nc"
ocean_data = xr.open_dataset(SIM_DIR / file).chunk(time=1)
seconds = ocean_data.time / pd.Timedelta("1s")
ocean_data = ocean_data.assign_coords(time=('time', seconds.data))

upper_limit = -1
isopycnal_numbers = 15
t1 = 200
t2 = 210
t3 = 300
t4 = 310
v_max = 1e-4
lin_thresh = 1e-6
# Colormap: "bwr" (Base_Case notebook; the preprint used RdYlBu). Chosen 2026-10-14.
c_map = "bwr"

fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 4), layout="constrained")
panels = [(axs[0, 0], t1, " ", "z [$m$]"), (axs[0, 1], t2, " ", " "),
          (axs[1, 0], t3, "x [$m$]", "z [$m$]"), (axs[1, 1], t4, "x [$m$]", " ")]
for ax, ti, xlab, ylab in panels:
    im = ocean_data.u.isel(time=ti, zC=slice(0, upper_limit)).plot(
        ax=ax, label=" ", add_colorbar=False, cmap=c_map,
        norm=colors.SymLogNorm(linthresh=lin_thresh, vmin=-v_max, vmax=v_max))
    (ocean_data.B + ocean_data.ba).isel(time=ti, yC=0, zC=slice(0, upper_limit)).plot.contour(
        ax=ax, colors="k", levels=isopycnal_numbers, alpha=0.75)
    ax.set(title=r"T=" + str(round(ocean_data.time[ti].values * fstar / (2 * np.pi), 2)), xlabel=xlab, ylabel=ylab)
fig.colorbar(im, ax=axs, label="u-velocity [$m\\,s^{-1}$]",
             ticks=[1e-2, -1e-2, 1e-3, -1e-3, 1e-4, -1e-4, 1e-5, -1e-5, 1e-6, -1e-6, 0])

save(fig, __file__)
