"""Floquet growth rate for each simulation, evaluated at its exact (S, gamma, delta).

Uses linear_theory/floquet.py: Eq. 2.26 with the Eq. 2.28 coefficients (time scaled
by f*, so sigma is in units of f*), the same system used for the parameter-space figures.
The growth rate is maximised over the perturbation slope mu0 directly rather than looked up
on a precomputed grid.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "linear_theory"))
import floquet as fl  # noqa: E402

from sims import SIMS, gamma  # noqa: E402

M_GRID = np.arange(-50, 5 + 0.1, 0.1)


def max_growth(S, gm, delta):
    return fl.fastest_mode(fl.Flow(S, gm, delta), M_GRID)


if __name__ == "__main__":
    print(f"{'sim':<6}{'S':>5}{'Ri_inv':>8}{'delta':>7}{'gamma':>8}{'sigma/f*':>10}{'Im/f*':>8}{'m_max':>8}")
    for name, s in SIMS.items():
        g = gamma(s["S"], s["Ri_inv"])
        sigma, freq, m = max_growth(s["S"], g, s["delta"])
        print(f"{name:<6}{s['S']:>5}{s['Ri_inv']:>8}{s['delta']:>7}{g:>8.4f}{sigma:>10.4f}{freq:>8.3f}{m:>8.2f}")
