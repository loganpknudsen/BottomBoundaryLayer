"""Simulation vs. Floquet growth rate (manuscript Fig. 10), from fit_growth.py output.

Choices (agreed 2026-10-14):
  * window tolerance p = 10%;
  * only PSI runs are shown: SI-unstable runs and Sim6 (not PSI in theory) are dropped;
  * a PSI run with no exponential phase (Sim18) is plotted at sigma = 0 with an open marker.

Usage: python make_figure.py growth_fits.csv [OUT_PREFIX]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats

from sims import NOT_PSI, SIMS

P = 0.10
MARKERS = {0.5: "x", 1.0: "+", 1.5: "*"}


def main(fits_csv, prefix="simulation_vs_stability"):
    df = pd.read_csv(fits_csv)
    df = df[np.isclose(df.p, P) & ~df.sim.isin(NOT_PSI)].copy()
    df["no_growth"] = df.status != "ok"
    df.loc[df.no_growth, "sigma"] = 0.0
    df[["sim", "S", "Ri_inv", "delta", "theory", "sigma", "r2", "T_start", "T_end", "no_growth"]] \
        .to_csv(f"{prefix}.csv", index=False)

    r2 = scipy.stats.pearsonr(df.sigma, df.theory).statistic ** 2
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    for S, m in MARKERS.items():
        sub = df[np.isclose(df.S, S) & ~df.no_growth]
        ax.scatter(sub.sigma, sub.theory, marker=m, c="k", label=rf"$S_\infty={S:g}$")
    ng = df[df.no_growth]
    ax.scatter(ng.sigma, ng.theory, marker="o", facecolors="none", edgecolors="k",
               label="no exponential phase")
    lim = 1.1 * max(df.sigma.max(), df.theory.max())
    ax.plot([0, lim], [0, lim], "k--", alpha=0.8)
    ax.set_xlim(-0.005, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.set_xlabel(r"Growth rate from simulations, $\sigma/f^*$")
    ax.set_ylabel(r"Growth rate from stability analysis, $\sigma/f^*$")
    ax.set_title(rf"$r^2={r2:.2f}$")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{prefix}.pdf")
    fig.savefig(f"{prefix}.png", dpi=200)
    print(df[["sim", "S", "Ri_inv", "delta", "sigma", "theory", "no_growth"]].round(4).to_string(index=False))
    print(f"r^2 = {r2:.3f} over {len(df)} runs; wrote {prefix}.pdf/.png/.csv")


if __name__ == "__main__":
    main(*sys.argv[1:])
