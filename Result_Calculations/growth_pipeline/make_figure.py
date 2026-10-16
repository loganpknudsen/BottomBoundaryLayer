"""Simulation vs. Floquet growth rate (manuscript Fig. 10), from fit_growth.py output.

Choices (agreed 2026-10-14):
  * fit window of n = 4 amplitude e-folds ending at the peak growth rate (fit_growth.py);
  * SI-unstable runs (Sim5, Sim10, Sim14) are included, drawn in grey (the Floquet growth rate
    is that of the fastest mode, SI for these runs); r^2 is printed with and without them;
  * Sim6 is kept although its fastest Floquet mode is not PSI;
  * runs with no exponential phase (Sim6 in the no-flux runs) are plotted at sigma = 0 with an
    open marker.

With --no-si, the SI-unstable runs are dropped.

Usage: python make_figure.py growth_fits.csv [OUT_PREFIX] [--no-si]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy.stats

from sims import SIMS

N_EFOLDS = 4
MARKERS = {0.5: "x", 1.0: "+", 1.5: "*"}
SI_COLOR = "0.6"


def main(fits_csv, prefix="simulation_vs_stability", include_si=True):
    df = pd.read_csv(fits_csv)
    df = df[df.n == N_EFOLDS].copy()
    if not include_si:
        df = df[~df.si_unstable]
    df["color"] = np.where(df.si_unstable, SI_COLOR, "k")
    df["no_growth"] = df.status != "ok"
    df.loc[df.no_growth, "sigma"] = 0.0
    df[["sim", "S", "Ri_inv", "delta", "theory", "sigma", "r2", "T_start", "T_end", "no_growth"]] \
        .to_csv(f"{prefix}.csv", index=False)

    r2 = scipy.stats.pearsonr(df.sigma, df.theory).statistic ** 2
    no_si = df[~df.si_unstable]
    r2_no_si = scipy.stats.pearsonr(no_si.sigma, no_si.theory).statistic ** 2
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    for S, m in MARKERS.items():
        sub = df[np.isclose(df.S, S) & ~df.no_growth]
        ax.scatter(sub.sigma, sub.theory, marker=m, c=sub.color, label=rf"$S_\infty={S:g}$")
    ng = df[df.no_growth]
    ax.scatter(ng.sigma, ng.theory, marker="o", facecolors="none", edgecolors=ng.color,
               label="no exponential phase")
    lim = 1.1 * max(df.sigma.max(), df.theory.max())
    ax.plot([0, lim], [0, lim], "k--", alpha=0.8)
    ax.set_xlim(-0.005, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.set_xlabel(r"Growth rate from simulations, $\sigma/f^*$")
    ax.set_ylabel(r"Growth rate from stability analysis, $\sigma/f^*$")
    ax.set_title(rf"$r^2={r2:.2f}$")
    handles, labels = ax.get_legend_handles_labels()
    if df.si_unstable.any():
        handles.append(Line2D([], [], ls="none", marker="o", color=SI_COLOR))
        labels.append("SI unstable")
    leg = ax.legend(handles, labels, fontsize=8)
    for h in leg.legend_handles[:len(MARKERS) + 1]:
        h.set_edgecolor("k")
    fig.tight_layout()
    fig.savefig(f"{prefix}.pdf")
    fig.savefig(f"{prefix}.png", dpi=200)
    print(df[["sim", "S", "Ri_inv", "delta", "sigma", "theory", "no_growth"]].round(4).to_string(index=False))
    print(f"r^2 = {r2:.3f} over {len(df)} runs (excl. SI: {r2_no_si:.3f} over {len(no_si)}); wrote {prefix}.pdf/.png/.csv")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--no-si"]
    main(*args, include_si="--no-si" not in sys.argv[1:])
