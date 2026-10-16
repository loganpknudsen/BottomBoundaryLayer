"""Figure 17: Sim9D (Sim9 parameters with quadratic bottom drag). Perturbation u with isopycnals at
(a) T = 5 and (b) T = 15; (c) quasi-Lagrangian power spectral density of u over T = 0-20.

Run: Sim9D, full_code_diagnostics_stblty_comp_drag.jl with S = 1, Ri_inv = 0.95, delta = 0.8
(c_d = 3.1e-3 from the log law, zero total buoyancy gradient at the bottom); output
flow_fields_Sinf_1.0_Ri_inv_0.95_delta_0.8_drag.nc in BBL_SIM_DIR.

Based on Result_Calculations/drag_figure.ipynb (origin/Base_Case, cell 15), redrawn: each snapshot
has its own linear colour scale (the amplitudes differ ~10x), limited to the 99th percentile of |u|
above Z_LIMIT_MIN, so the drag-driven turbulence next to the bed (which saturates) does not wash out
the PSI structure; the spectrum is computed as in
Fig 12 (averaged over all starting x, per unit angular frequency), at z = 50.39 m (level 64).
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

P = ql.background(Sinf=1.0, Ri_inv=0.95, delta=0.8)
FILE = "flow_fields_Sinf_1.0_Ri_inv_0.95_delta_0.8_drag.nc"
SNAPSHOTS = (100, 300)        # outputs every 0.05 inertial periods: T = 5, 15
SPECTRUM = (0, 400)           # T = 0-20
ZI = 64                       # z = 50.39 m
PERCENTILE = 99               # colour limit for each snapshot ...
Z_LIMIT_MIN = 10.0            # ... computed from |u| above this height [m]
N_ISOPYCNALS = 15

ds = xr.open_dataset(SIM_DIR / FILE)
ds = ds.assign_coords(time=("time", (ds.time / pd.Timedelta("1s")).values))
T = ds.time.values * P["fstar"] / (2 * np.pi)


def nice_limit(x):
    """Round up to 1, 2 or 5 x 10^n."""
    e = np.floor(np.log10(x))
    return next(float(m * 10**e) for m in (1, 2, 5, 10) if m * 10**e >= x)


# Snapshots stacked as wide panels (the domain is 2000 m x 200 m), sharing the x axis. One grid
# with a narrow colourbar column, so the spectrum below has the same width as the snapshots.
fig = plt.figure(figsize=(7.5, 8), layout="constrained")
fig.get_layout_engine().set(h_pad=0.01, hspace=0)   # tighter between (a) and (b)
gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 0.1, 1.1], width_ratios=[1, 0.025])  # row 2: spacer
axs = {"s0": fig.add_subplot(gs[0, 0])}
axs["s1"] = fig.add_subplot(gs[1, 0], sharex=axs["s0"])
axs["psd"] = fig.add_subplot(gs[3, 0])
caxs = {"s0": fig.add_subplot(gs[0, 1]), "s1": fig.add_subplot(gs[1, 1])}

for key, ti in zip(("s0", "s1"), SNAPSHOTS):
    ax = axs[key]
    u = ds.u.isel(time=ti).squeeze(drop=True)
    vmax = nice_limit(float(np.nanpercentile(np.abs(u.where(u.zC > Z_LIMIT_MIN).values), PERCENTILE)))
    print(f"T = {T[ti]:5.2f}: max |u| = {float(np.abs(u).max()):.2e} m/s, colour limit +/-{vmax:.0e} m/s")
    im = u.plot(ax=ax, cmap="RdBu_r", vmin=-vmax, vmax=vmax, add_colorbar=False, rasterized=True)
    (ds.B + ds.ba).isel(time=ti).squeeze(drop=True).plot.contour(
        ax=ax, colors="k", levels=N_ISOPYCNALS, linewidths=0.6, linestyles="solid", alpha=0.7)
    fig.colorbar(im, cax=caxs[key], label=r"$u$ [m s$^{-1}$]", extend="both")
    ax.set(title="", xlabel="x [m]" if key == "s1" else "", ylabel="z [m]")
    if key == "s0":
        ax.tick_params(labelbottom=False)
    ax.text(0.01, 0.95, f"T = {T[ti]:.0f}", transform=ax.transAxes, ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))

u_ql = ql.quasi_lagrangian(ds, "u", ZI, *SPECTRUM, P)
omega, psd = ql.spectrum(u_ql)
pos = omega > 0
ax = axs["psd"]
ax.plot(omega[pos], psd[pos], color="k", lw=1.2)
ax.axvline(0.5 * P["fstar"], color="tab:blue", ls=":", lw=1.5, label=r"$f^*/2$")
ax.axvline(P["fstar"], color="tab:orange", ls=":", lw=1.5, label=r"$f^*$")
ax.set_xscale("log")
ax.set_xlim(float(omega[pos].min()), float(omega[pos].max()))
ax.set_ylim(0, None)
ax.set_xlabel(r"Angular frequency $\omega$ [rad s$^{-1}$]")
ax.set_ylabel(r"PSD [m$^2$ s$^{-1}$ rad$^{-1}$]")
ax.set_title(f"Quasi-Lagrangian $u$, z = {float(ds.zC[ZI]):.1f} m, T = {T[SPECTRUM[0]]:.0f}-{T[SPECTRUM[1]]:.0f}",
             fontsize="medium")
ax.legend(frameon=False)

for key, label in (("s0", "(a.)"), ("s1", "(b.)"), ("psd", "(c.)")):
    axs[key].text(-0.1, 1.04, label, transform=axs[key].transAxes, fontsize="medium")

save(fig, __file__)
