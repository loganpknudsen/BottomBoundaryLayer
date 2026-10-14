"""Shared pieces for the linear-theory figures (Figs 3-6, 8, 9).

Reads the grids written by linear_theory/stability_grid.py. Plotted fields pass through mask(),
which holds the masking criteria agreed so far:
  1. points above the convective limit (convectively unstable at some wave phase) are removed;
  2. below gamma_l (omega_min^2/f*^2 > 0.25), points with sigma/f below GROWTH_CUTOFF_F are
     removed (the weak modes there depend on the m-scan resolution; above gamma_l they are
     converged and are kept);
  3. below gamma_l, points whose fastest mode is not PSI are removed. SI cannot occur there, so
     this removes only non-PSI (harmonic) growth; PSI and SI above gamma_l are kept.
The analytic curves (gamma_l, gamma_u, convective limit) are drawn as lines. The PSI boundary
line is the 0.5 contour of an unsmoothed indicator: 1 where the fastest mode is PSI, growing
(sigma/f >= GROWTH_CUTOFF_F below gamma_l, sigma/f* > GROWTH_FLOOR above it) and below the
convective limit, 0 elsewhere; its accuracy is set by the grid spacing (see linear_theory REFINE). The TKE-budget
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

# PSI: |Im sigma| = f*/2, tested as >= PSI_FREQ.
PSI_FREQ = 0.4999
# Growth-rate cutoff below gamma_l, in units of f (sigma/f = (sigma/f*) sqrt(1 + S^2), as
# plotted). Removes the weak, m-resolution-dependent modes there (artefacts of the 0.1 m scan reach
# ~0.02-0.03 f); below gamma_l the PSI boundary therefore marks sigma = GROWTH_CUTOFF_F * f.
GROWTH_CUTOFF_F = 0.025
# "Growing" above gamma_l: stable points have |sigma/f*| < 4e-11 (round-off) on the grids checked,
# the weakest real growth found was 2e-6.
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


def growth_over_f(growth, S):
    """Growth rate in units of f from the stored sigma/f*."""
    return growth * np.sqrt(1 + S**2)


def mask(field, strat, S, delta, growth, freq):
    """Return `field` with masked points set to NaN. strat (gamma), S, delta, growth (sigma/f*)
    and freq (|Im sigma|/f*) must broadcast to field's shape. Criteria:
      1. gamma above the convective limit for this delta;
      2. gamma below gamma_l and sigma/f below GROWTH_CUTOFF_F;
      3. gamma below gamma_l and the fastest mode not PSI.
    """
    below = strat < gamma_lower_limit(S)
    keep = ((strat <= convective_limit(S, delta))
            & ~(below & (growth_over_f(growth, S) < GROWTH_CUTOFF_F))
            & ~(below & (np.abs(freq) < PSI_FREQ)))
    return np.where(keep, field, np.nan)


def psi_unstable(freq, growth, S, strat, delta):
    """1.0 where the fastest mode is PSI, growing (sigma/f >= GROWTH_CUTOFF_F below gamma_l,
    sigma/f* > GROWTH_FLOOR above it) and below the convective limit, else 0.0; contour at
    BOUNDARY_LEVEL. growth is sigma/f*."""
    below = strat < gamma_lower_limit(S)
    growing = np.where(below, growth_over_f(growth, S) >= GROWTH_CUTOFF_F, growth > GROWTH_FLOOR)
    return ((np.abs(freq) >= PSI_FREQ) & growing & (strat <= convective_limit(S, delta))).astype(float)


def sgamma_mesh(data, delta):
    """(S, gamma) 2-D arrays and the selected slice of `data` at one delta value."""
    d = data.sel(delta=delta)
    g = d.strat_values.values
    S = np.broadcast_to(d.slope_burger_number.values[:, None], g.shape)
    return S, g, d


def sgamma_boundary(d):
    """PSI-unstable indicator for an S-gamma slice from sgamma_mesh."""
    return psi_unstable(d.frequency.values, d.growth_rate.values, d.slope_burger_number.values[:, None],
                        d.strat_values.values, float(d.delta))


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
