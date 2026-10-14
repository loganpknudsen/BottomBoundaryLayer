"""Quasi-Lagrangian sampling (Eq. 3.3) and its power spectrum, shared by Figs 12 and 17."""
import numpy as np
import xarray as xr
import xrft as xf


def background(Sinf, Ri_inv, delta, theta=0.1, f=1e-4, Vinf=0.05):
    """Background parameters as defined in the simulations (full_code_diagnostics_*.jl)."""
    N2 = Sinf**2 * f**2 / np.tan(theta) ** 2
    gamma = (1 + (1 - Ri_inv ** (-1)) * Sinf**2) * (1 + Sinf**2) ** (-1)
    lmbd = N2 * gamma * np.tan(theta) / (f * np.cos(theta))
    return dict(delta=delta, lmbd=lmbd, h=Vinf / lmbd, fstar=f * np.cos(theta) * (1 + Sinf**2) ** 0.5)


def chi(t, z, p):
    """Background-wave displacement at height z (the QL shift of Eq. 3.3), t in seconds."""
    return p["delta"] * p["lmbd"] / p["fstar"] * np.sin(p["fstar"] * t + np.pi / 2) * (p["h"] - z)


def quasi_lagrangian(ds, field, zi, ti, tf, p):
    """(time, x0) DataArray: `field` at level zi following the QL position of every starting grid
    point x0, outputs ti:tf. On the uniform periodic grid the nearest grid point to x0 + chi(t) is
    x0's index + round(chi/dx) (wrapping in x). ds.time must be in seconds."""
    t = ds.time.values[ti:tf]
    dx = float(ds.xF[1] - ds.xF[0])
    n_shift = np.rint(chi(t, float(ds.zC[zi]), p) / dx).astype(int)
    a = ds[field].isel(time=slice(ti, tf), zC=zi).squeeze(drop=True).transpose("time", ...).values
    nx = a.shape[1]
    idx = (np.arange(nx)[None, :] + n_shift[:, None]) % nx
    return xr.DataArray(np.take_along_axis(a, idx, axis=1), dims=("time", "x0"),
                        coords={"time": t, "x0": ds.xF.values})


def spectrum(series):
    """Power spectral density of a (time, x0) QL series averaged over x0, per unit angular
    frequency: returns (omega [rad/s], psd [field units^2 s rad^-1])."""
    p_hz = xf.power_spectrum(series, dim="time").mean("x0")   # density per Hz (cycles/s)
    return p_hz.freq_time * 2 * np.pi, p_hz / (2 * np.pi)
