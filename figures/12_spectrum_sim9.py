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
Panel (a) averages the quasi-Lagrangian power spectrum over all starting x, and panel (b) follows the
starting x where |u| is largest at T = 5 (a crest). The notebook used the single point x = 500 m,
which in the 2026 runs sits on a node of the u pattern, so its series flips sign as the node drifts
across it.

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

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import SIM_DIR, save  # noqa: E402
import quasi_lagrangian as ql  # noqa: E402

P = ql.background(Sinf=1.0, Ri_inv=0.95, delta=0.8)   # Sim9
fstar = P["fstar"]

file = "flow_fields_Sinf_1.0_Ri_inv_0.95_delta_0.8.nc"
ocean_data = xr.open_dataset(SIM_DIR / file)
seconds = ocean_data.time / pd.Timedelta("1s")
ocean_data = ocean_data.assign_coords(time=('time', seconds.data))

# Quasi-Lagrangian velocities from every starting x (figures/quasi_lagrangian.py)
ti = 100   # T = 5
tf = 300   # T = 15
zi = int(sys.argv[sys.argv.index("--zi") + 1]) if "--zi" in sys.argv else 64
u_ql = ql.quasi_lagrangian(ocean_data, "u", zi, ti, tf, P)
v_ql = ql.quasi_lagrangian(ocean_data, "v", zi, ti, tf, P)

# (a) power spectrum averaged over all starting positions, so it does not depend on where the
# sample sits relative to the crests and nodes of the perturbation pattern
omega, psd = ql.spectrum(u_ql)

# (b) time series at a crest: the starting x where |u| is largest at T = 5
j_crest = int(np.argmax(np.abs(u_ql.isel(time=0).values)))
u_lagrange, v_lagrange = u_ql.isel(x0=j_crest), v_ql.isel(x0=j_crest)
z_label = f"z={float(ocean_data.zC[zi]):.2f} m, $x_0$={float(u_ql.x0[j_crest]):.0f} m"

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
