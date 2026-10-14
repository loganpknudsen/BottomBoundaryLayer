"""Shared pieces for the linear-theory figures (Figs 3-6, 8, 9).

Reads the grids written by linear_theory/stability_grid.py. Plotted fields pass through mask(),
which holds the masking criteria agreed so far:
  1. points above the convective limit (convectively unstable at some wave phase) are removed.
The analytic curves (gamma_l, gamma_u, convective limit) are drawn as lines. Two choices are kept
from stability_analysis_output_reader.ipynb as placeholders, to be settled with the masking:
  * the PSI boundary line: the 0.4999 contour of the smoothed frequency field (cells 15, 18);
  * the TKE-budget variant for Figs 8-9 (BUDGET_IC).
"""
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import LT_DIR  # noqa: E402

# "ic01" starts the budget from (Psi, zeta) = (0, 1) as in the preprint; "floquet" uses the
# growing Floquet mode. Placeholder: not yet decided.
BUDGET_IC = "ic01"

# Placeholder PSI boundary (notebook cells 15, 18): frequency set to 0 above 1.1 gamma_u,
# smoothed with this sigma (in grid cells), contoured at BOUNDARY_LEVEL.
SMOOTH_SIGMA = 0.275
BOUNDARY_LEVEL = 0.4999


def load(name):
    """'sgamma' (Figs 3-5, 8, 9) or 'freq' (Fig 6)."""
    return xr.open_dataset(LT_DIR / f"{name}.nc")


def gamma_lower_limit(S):
    return (3 - S**2) / (3 * (1 + S**2))


def gamma_upper_limit(S):
    return (1 + S**2) ** (-1)


def convective_limit(S, d):
    return np.sqrt(1 + S**2) / (np.sqrt(1 + S**2) + d * S**2)


def mask(field, strat, S, delta):
    """Return `field` with masked points set to NaN. strat (gamma), S and delta must broadcast to
    field's shape. Criteria:
      1. gamma above the convective limit for this delta.
    """
    keep = strat <= convective_limit(S, delta)
    return np.where(keep, field, np.nan)


def boundary_field(freq, strat, S):
    """Smoothed frequency used for the PSI boundary contour (placeholder method).
    freq, strat: 2-D arrays on the same grid; S: slope Burger number(s) broadcastable to them."""
    f = np.where(strat < 1.1 * gamma_upper_limit(S), freq, 0.0)
    return gaussian_filter(f, SMOOTH_SIGMA)


def sgamma_mesh(data, delta):
    """(S, gamma) 2-D arrays and the selected slice of `data` at one delta value."""
    d = data.sel(delta=delta)
    g = d.strat_values.values
    S = np.broadcast_to(d.slope_burger_number.values[:, None], g.shape)
    return S, g, d


def sgamma_boundary(d):
    """Boundary field for an S-gamma slice from sgamma_mesh."""
    S = d.slope_burger_number.values[:, None]
    return boundary_field(d.frequency.values, d.strat_values.values, S)


def budget(d, term):
    """Budget fraction for the chosen BUDGET_IC; term is GSP, ASP or BP (sum of components)."""
    parts = {"GSP": ["GSP"], "ASP": ["ASPx", "ASPy"], "BP": ["BPx", "BPz"]}[term]
    return sum(d[f"{p}_{BUDGET_IC}"] for p in parts).values


def analytic_lines(ax, S, delta=None):
    """gamma_u and gamma_l (black) and, if delta is given, the convective limit (grey)."""
    ax.plot(S, gamma_upper_limit(S), c="k")
    ax.plot(S, gamma_lower_limit(S), c="k")
    if delta is not None:
        ax.plot(S, convective_limit(S, delta), c="grey", linestyle="--", linewidth=2)
