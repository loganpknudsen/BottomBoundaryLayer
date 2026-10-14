"""Figure 13: (a) perturbation u for Sim9 as the flow reaches finite amplitude (T = 17, 17.5, 18,
18.5), with isopycnals and Ri <= 0.25 contours at T = 17, 17.5; (b) fraction of BBL grid points
with Ri <= 0.25 against time.

Source: Result_Calculations/velocity_visualization.ipynb (origin/Base_Case), cells 1, 5, 6, 8, 10,
14 and 28. Ri is computed exactly as in the notebook (N^2 / V_z^2 from b + B and v + vb).
Run: Sim9 (S = 1, Ri_inv = 0.95, delta = 0.8), March output flow_fields_..._Ri_inv_0.95_delta_0.8.nc.
"""
import pathlib
import sys

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import SIM_DIR, save  # noqa: E402

# Run parameters (Sim9), cell 1
theta = 0.1
f = 1e-4
Sinf = 1.0
Ri_inv = 0.95
N2 = Sinf**2 * f**2 / np.tan(theta) ** 2
delta = 0.8
Vinf = 0.05
gamma = (1 + (1 - Ri_inv ** (-1)) * Sinf**2) * (1 + Sinf**2) ** (-1)
lmbd = N2 * gamma * np.tan(theta) / (f * np.cos(theta))
h = Vinf / lmbd
fstar = f * np.cos(theta) * (1 + Sinf**2) ** (0.5)

file = "flow_fields_Sinf_1.0_Ri_inv_0.95_delta_0.8.nc"
ocean_data = xr.open_dataset(SIM_DIR / file).chunk(time=1)
seconds = ocean_data.time / pd.Timedelta("1s")
ocean_data = ocean_data.assign_coords(time=('time', seconds.data))

# Richardson number (cells 6 and 8)
B_BG = ocean_data.B.isel(yC=0)
B = ocean_data.b.isel(yC=0) + B_BG
N2 = B.differentiate("zC", edge_order=2)
V_BG = ocean_data.vb.isel(yF=0)
V = ocean_data.v.isel(yF=0) + V_BG
VZ = V.differentiate("zC", edge_order=2)
Ri = N2 / VZ**2
Ri = Ri.sel(zC=slice(0, h))

# Fraction of BBL points with Ri <= 0.25 for the first 501 outputs (cells 10 and 14). The notebook
# counts these with a Python loop over Ri.where(Ri <= 0.25, other="nan"); summing the boolean mask
# gives the same count. 135 is the number of grid levels below h for Sim9.
count_list = (Ri.isel(time=slice(0, 501)) <= 0.25).sum(dim=[d for d in Ri.dims if d != "time"]).values
prcnt = np.array(count_list) / (135 * 1024)

fig = plt.figure(figsize=(10, 7), constrained_layout=True)
gs = gridspec.GridSpec(3, 2, figure=fig, height_ratios=[1.2, 1.2, 1])
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])
ax_big = fig.add_subplot(gs[2, :])

upper_limit = -1
isopycnal_numbers = 15
t1 = 340
t2 = 350
t3 = 360
t4 = 370
v_max = 2e-3
# TODO(decision): colormap. Base_Case cell 28 uses "bwr"; the submission commit (27f8414, cell 27)
# uses "RdYlBu", and the preprint Fig 13 colour bar (blue at +0.002, red at -0.002) looks like
# RdYlBu. Panel (b), however, matches only the Base_Case version (BBL fraction, y-limit 0.41).
c_map = None
if c_map is None:
    raise NotImplementedError('Fig 13 colormap: "bwr" (Base_Case notebook) or "RdYlBu" (submission commit, matches preprint)?')


def T(ti):
    return str(round(ocean_data.time[ti].values * fstar / (2 * np.pi), 2))


ocean_data.u.isel(time=t1, zC=slice(0, upper_limit)).plot(ax=ax1, label=" ", add_colorbar=False, cmap=c_map, vmin=-1 * v_max, vmax=v_max)
(ocean_data.B + ocean_data.ba).isel(time=t1, yC=0, zC=slice(0, upper_limit)).plot.contour(ax=ax1, colors="k", levels=isopycnal_numbers, alpha=0.75)
Ri.isel(time=t1, zC=slice(0, upper_limit)).plot.contour(ax=ax1, levels=[0.25], colors="k", linewidths=1, linestyles="--")
ax1.set(title=r"T=" + T(t1), xlabel=" ", ylabel="z [$m$]")
ocean_data.u.isel(time=t2, zC=slice(0, upper_limit)).plot(ax=ax2, label=" ", add_colorbar=False, cmap=c_map, vmin=-1 * v_max, vmax=v_max)
Ri.isel(time=t2, zC=slice(0, upper_limit)).plot.contour(ax=ax2, levels=[0.25], colors="k", linewidths=1, linestyles="--")
(ocean_data.B + ocean_data.ba).isel(time=t2, yC=0, zC=slice(0, upper_limit)).plot.contour(ax=ax2, colors="k", levels=isopycnal_numbers, alpha=0.75)
ax2.set(title=r"T=" + T(t2), xlabel=" ", ylabel=" ")
ocean_data.u.isel(time=t3, zC=slice(0, upper_limit)).plot(ax=ax3, label=" ", add_colorbar=False, cmap=c_map, vmin=-1 * v_max, vmax=v_max)
(ocean_data.B + ocean_data.ba).isel(time=t3, yC=0, zC=slice(0, upper_limit)).plot.contour(ax=ax3, colors="k", levels=isopycnal_numbers, alpha=0.75)
ax3.set(title=r"T=" + T(t3), xlabel="x [$m$]", ylabel="z [$m$]")
im = ocean_data.u.isel(time=t4, zC=slice(0, upper_limit)).plot(ax=ax4, label=" ", add_colorbar=False, cmap=c_map, vmin=-1 * v_max, vmax=v_max)
(ocean_data.B + ocean_data.ba).isel(time=t4, yC=0, zC=slice(0, upper_limit)).plot.contour(ax=ax4, colors="k", levels=isopycnal_numbers, alpha=0.75)
ax4.set(title=r"T=" + T(t4), xlabel="x [$m$]", ylabel=" ")
fig.colorbar(im, ax=[ax1, ax2, ax3, ax4], label=r"u-velocity [$m\text{ }s^{-1}$]",
             ticks=[5e-3, -5e-3, 4e-3, -4e-3, 3e-3, -3e-3, 2e-3, -2e-3, -1e-3, 1e-3, 0])
ax1.text(-200, 185, "(a.)")

ax_big.plot(ocean_data.time.values[0:501] * fstar / (2 * np.pi), prcnt, color="k")
ax_big.vlines(ocean_data.time.values[t1] * fstar / (2 * np.pi), -0.01, 0.42, color="grey", linestyle="-.", label=r"T=" + T(t1))
ax_big.vlines(ocean_data.time.values[t2] * fstar / (2 * np.pi), -0.01, 0.42, color="grey", linestyle="--", label=r"T=" + T(t2))
ax_big.legend()
ax_big.text(-1, 0.45, "(b.)")
ax_big.set_xlim([0, 25])
ax_big.set_ylim([-0.01, 0.41])
ax_big.set_xlabel("T")
ax_big.set_ylabel(r"Fraction w/ Ri$\leq 0.25$")

save(fig, __file__)
