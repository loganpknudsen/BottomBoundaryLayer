"""Diagnostic plots of TKE growth for the March runs, to design the fit-window logic.

For each simulation this plots
  (top)    ln k vs T = f_s t / 2pi
  (bottom) the one-period growth rate  sigma_1P(T) = [ln k(T+1) - ln k(T)] / (2 * 2pi),
           in units of f_s, with the Floquet prediction overlaid.

TKE from a Floquet mode is exp(2 sigma t) times a function periodic in the wave period,
so differencing ln k over exactly one wave period removes the intra-period modulation and
leaves the period-averaged growth rate; a fit must likewise span whole periods.

Usage: python diagnose_growth.py DATA_DIR [OUT_PREFIX]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from sims import SIMS, tke_filename

THETA = 0.1
F = 1e-4
# Floquet growth rates (sigma / f*) from theory_growth.py (Eq 2.28 system)
THEORY = {"Sim1": 0.1346, "Sim2": 0.1768, "Sim3": 0.1827, "Sim4": 0.1843, "Sim5": 0.2403,
          "Sim6": 0.0126, "Sim7": 0.1074, "Sim8": 0.1292, "Sim9": 0.1453, "Sim10": 0.2834,
          "Sim11": 0.0309, "Sim12": 0.0713, "Sim13": 0.1027, "Sim14": 0.2655, "Sim15": 0.1200,
          "Sim16": 0.0873, "Sim17": 0.0571, "Sim18": 0.0282}


def load_lnk(path, S):
    """Return T (wave periods) and ln k, with k the domain-averaged TKE."""
    k = xr.open_dataset(path, decode_timedelta=True).k.squeeze()
    fs = F * np.cos(THETA) * np.sqrt(1 + S**2)
    t = (k.time / pd.Timedelta("1s")).values
    return t * fs / (2 * np.pi), np.log(k.values)


def one_period_growth(T, lnk):
    """sigma_1P(T) in units of f_s, using linear interpolation to land exactly one period later."""
    ok = T + 1 <= T[-1]
    return T[ok], (np.interp(T[ok] + 1, T, lnk) - lnk[ok]) / (4 * np.pi)


def main(data_dir, prefix="growth_diagnostics"):
    data_dir = Path(data_dir)
    names = [n for n in SIMS if (data_dir / tke_filename(SIMS[n])).exists()]
    missing = [n for n in SIMS if n not in names]
    if missing:
        print("missing files for:", ", ".join(missing))

    ncol = 3
    nrow = -(-len(names) // ncol)
    fig, axs = plt.subplots(2 * nrow, ncol, figsize=(4.2 * ncol, 4.2 * nrow), squeeze=False)
    rows = []
    for i, name in enumerate(names):
        s = SIMS[name]
        T, lnk = load_lnk(data_dir / tke_filename(s), s["S"])
        Tg, sig = one_period_growth(T, lnk)
        rows.append(pd.DataFrame({"sim": name, "T": Tg, "sigma_1P": sig}))

        r, c = 2 * (i // ncol), i % ncol
        tag = " (SI)" if s["si_unstable"] else ""
        axs[r, c].plot(T, lnk - lnk[0], "k", lw=0.8)
        axs[r, c].set_title(f"{name}{tag}: S={s['S']}, Ri$^{{-1}}$={s['Ri_inv']}, $\\delta$={s['delta']}", fontsize=9)
        axs[r, c].set_ylabel(r"$\ln k - \ln k_0$")
        axs[r + 1, c].plot(Tg, sig, "k", lw=0.8)
        axs[r + 1, c].axhline(THEORY[name], color="r", ls="--", lw=1, label="Floquet")
        axs[r + 1, c].axhline(0, color="0.6", lw=0.5)
        axs[r + 1, c].set_ylim(-0.1, max(0.35, 1.2 * THEORY[name]))
        axs[r + 1, c].set_ylabel(r"$\sigma_{1P}/f^*$")
        axs[r + 1, c].set_xlabel(r"$T = f^* t / 2\pi$")
        for ax in axs[r:r + 2, c]:
            ax.set_xlim(0, T[-1])
    for j in range(len(names), nrow * ncol):
        for ax in axs[2 * (j // ncol):2 * (j // ncol) + 2, j % ncol]:
            ax.set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{prefix}.png", dpi=110)
    pd.concat(rows).to_csv(f"{prefix}.csv", index=False)
    print(f"wrote {prefix}.png and {prefix}.csv")


if __name__ == "__main__":
    main(*sys.argv[1:])
