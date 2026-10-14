"""Figure 12: Sim9 (a) power spectral density of the quasi-Lagrangian u velocity; (b) quasi-Lagrangian
u (and v) normalised by u after 5 modified inertial periods, symmetric-log scale.

Source: Result_Calculations/time_series_and_fourier_transform.ipynb, cells 1, 3-7, 15 and 17
(origin/Base_Case numbering; cell 16 in the submission commit 27f8414).
Run: Sim9 (S = 1, Ri_inv = 0.95, delta = 0.8), flow_fields_..._Ri_inv_0.95_delta_0.8.nc.

Follows the submission-commit version (chosen 2026-10-14): the spectrum covers T = 5-15 (after the
initial transient) and panel (b) shows u and v. Fixes relative to the notebook: the legend depth is
taken from the data (the notebook's "z=16.02 m" was stale; the series is at grid level 64,
z = 50.39 m for Sim9); the PSD, which xrft returns per Hz, is divided by 2 pi so that it is a density
per unit angular frequency, consistent with the omega axis; axis labels give the correct units.

Usage: python 12_spectrum_sim9.py [--zi N]   (grid level of the time series; default 64, z = 50.39 m
for Sim9; earlier notebook versions used 20, z = 16.02 m). A non-default level is added to the
output name, e.g. 12_spectrum_sim9_zi20.png.
"""
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import xrft as xf

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
ocean_data = xr.open_dataset(SIM_DIR / file)
seconds = ocean_data.time / pd.Timedelta("1s")
ocean_data = ocean_data.assign_coords(time=('time', seconds.data))


def z_change(t, z, x):  # advection by isopycnals
    z_osc = delta * lmbd / fstar * np.sin(fstar * t + np.pi / 2) * (h - z)
    return z_osc


# Quasi-Lagrangian velocities (cells 6 and 7)
ti = 100   # T = 5
tf = 300   # T = 15
zi = int(sys.argv[sys.argv.index("--zi") + 1]) if "--zi" in sys.argv else 64
xi = 256
v_lagrange = ocean_data.v.sel(xC=ocean_data.xC[xi] + z_change(ocean_data.time[ti:tf], ocean_data.zC[zi], ocean_data.xC[xi]),
                              zC=ocean_data.zC[zi], time=ocean_data.time[ti:tf], method="nearest")
v_lagrange = v_lagrange.drop_vars(['xC', 'yF', 'zC'])
u_lagrange = ocean_data.u.sel(xF=ocean_data.xF[xi] + z_change(ocean_data.time[ti:tf], ocean_data.zC[zi], ocean_data.xF[xi]),
                              zC=ocean_data.zC[zi], time=ocean_data.time[ti:tf], method="nearest")
u_lagrange = u_lagrange.drop_vars(['xF', 'yC', 'zC'])

transform_lagrange = xf.power_spectrum(u_lagrange, dim="time")   # density per Hz (cycles/s)
omega = transform_lagrange.freq_time * 2 * np.pi                   # angular frequency [rad/s]
psd = transform_lagrange / (2 * np.pi)                             # density per unit omega
z_label = f"z={float(ocean_data.zC[zi]):.2f} m"

fig, ax_dict = plt.subplot_mosaic([['(a.)', '(b.)']], layout="constrained", figsize=(9, 4))
ax = ax_dict['(b.)']
ax.plot(u_lagrange.time[:-1] * (fstar) / (2 * np.pi), u_lagrange.values[:-1] / u_lagrange.values[0], color="r", label="$u/u_{5}$ at " + z_label)
ax.plot(v_lagrange.time[:-1] * (fstar) / (2 * np.pi), v_lagrange.values[:-1] / u_lagrange.values[0], color="k", alpha=0.8, label="$v/u_{5}$ at " + z_label)
ax.set_xlim([u_lagrange.time[0] * (fstar) / (2 * np.pi), u_lagrange.time[-1] * (fstar) / (2 * np.pi)])
legend_loc = "lower left"
ax.set_xlabel("Modified Inertial Period")
ax.set_yscale("symlog", linthresh=1)
ax.legend(loc=legend_loc)
ax.set_box_aspect(0.75)
ax.set_title("(b.)", fontfamily='serif', loc='left', fontsize='medium')

ax = ax_dict['(a.)']
ax.plot(omega, psd, color="k")
ax.vlines(0.5 * fstar, ymin=0, ymax=float(psd.max()) * 1.05, color="b", alpha=0.5, linestyles="dotted", label="$f^*/2$")
ax.vlines(fstar, ymin=0, ymax=float(psd.max()) * 1.05, color="k", alpha=0.5, linestyles="dotted", label="$f^*$")
ax.set_xlim([1e-05, float(omega.max())])
ax.set_ylim([0, float(psd.max()) * 1.05])
ax.set_box_aspect(0.9)
ax.set_title("(a.)", fontfamily='serif', loc='left', fontsize='medium')
ax.legend()
ax.set_xlabel(r"Angular frequency $\omega$ [rad s$^{-1}$]")
ax.set_xscale("log")
ax.set_ylabel(r"PSD [m$^2$ s$^{-1}$ rad$^{-1}$]")

save(fig, __file__ if zi == 64 else str(pathlib.Path(__file__).with_name(f"12_spectrum_sim9_zi{zi}.py")))
