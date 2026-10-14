"""Shared pieces for the linear-theory figures (Figs 3-6, 8, 9).

Reads the grids written by linear_theory/stability_grid.py. Plotted fields pass through mask(),
which holds the masking criteria agreed so far:
  1. points above the convective limit (convectively unstable at some wave phase) are removed.
The analytic curves (gamma_l, gamma_u, convective limit) are drawn as lines. The PSI boundary
line is the 0.5 contour of an unsmoothed indicator: 1 where the fastest mode is PSI and growing,
0 elsewhere; its accuracy is set by the grid spacing (see linear_theory REFINE). The TKE-budget
variant for Figs 8-9 (BUDGET_IC) is still a placeholder kept from the preprint.
"""
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import LT_DIR  # noqa: E402

# "ic01" starts the budget from (Psi, zeta) = (0, 1) as in the preprint; "floquet" uses the
# growing Floquet mode. Placeholder: not yet decided.
BUDGET_IC = "ic01"

# PSI-unstable: |Im sigma| = f*/2 (tested as >= PSI_FREQ) and growth above GROWTH_FLOOR. Stable
# points have |growth| < 4e-11 (round-off) on the grids checked; the weakest real growth found was
# 2e-6, so the floor sits in the gap between them.
PSI_FREQ = 0.4999
GROWTH_FLOOR = 1e-8
BOUNDARY_LEVEL = 0.5


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


def psi_unstable(freq, growth):
    """1.0 where the fastest mode is PSI and growing, else 0.0; contour at BOUNDARY_LEVEL."""
    return ((np.abs(freq) >= PSI_FREQ) & (growth > GROWTH_FLOOR)).astype(float)


def sgamma_mesh(data, delta):
    """(S, gamma) 2-D arrays and the selected slice of `data` at one delta value."""
    d = data.sel(delta=delta)
    g = d.strat_values.values
    S = np.broadcast_to(d.slope_burger_number.values[:, None], g.shape)
    return S, g, d


def sgamma_boundary(d):
    """PSI-unstable indicator for an S-gamma slice from sgamma_mesh."""
    return psi_unstable(d.frequency.values, d.growth_rate.values)


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
