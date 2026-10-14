"""Simulation vs. Floquet growth rate (manuscript Fig. 10), from fit_growth.py output.

Choices (agreed 2026-10-14):
  * window tolerance p = 10%;
  * only PSI runs are shown: SI-unstable runs and Sim6 (not PSI in theory) are dropped;
  * a PSI run with no exponential phase (Sim18) is plotted at sigma = 0 with an open marker.

With --all, the excluded runs are added back for testing, coloured by why they would be
excluded (red: SI-unstable, blue: not PSI in theory), and r^2 is reported both ways.

Usage: python make_figure.py growth_fits.csv [OUT_PREFIX] [--all]
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


def main(fits_csv, prefix="simulation_vs_stability", include_all=False):
    df = pd.read_csv(fits_csv)
    df = df[np.isclose(df.p, P)].copy()
    if not include_all:
        df = df[~df.sim.isin(NOT_PSI)]
    df["color"] = np.where(df.si_unstable, "tab:red", np.where(df.sim.isin(NOT_PSI), "tab:blue", "k"))
    df["no_growth"] = df.status != "ok"
    df.loc[df.no_growth, "sigma"] = 0.0
    df[["sim", "S", "Ri_inv", "delta", "theory", "sigma", "r2", "T_start", "T_end", "no_growth"]] \
        .to_csv(f"{prefix}.csv", index=False)

    r2 = scipy.stats.pearsonr(df.sigma, df.theory).statistic ** 2
    psi = df[~df.sim.isin(NOT_PSI)]
    r2_psi = scipy.stats.pearsonr(psi.sigma, psi.theory).statistic ** 2
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    for S, m in MARKERS.items():
        sub = df[np.isclose(df.S, S) & ~df.no_growth]
        ax.scatter(sub.sigma, sub.theory, marker=m, c=sub.color, label=rf"$S_\infty={S:g}$")
    ng = df[df.no_growth]
    ax.scatter(ng.sigma, ng.theory, marker="o", facecolors="none", edgecolors=ng.color,
               label="no exponential phase")
    if include_all:
        for _, r in df[df.sim.isin(NOT_PSI)].iterrows():
            ax.annotate(r.sim, (r.sigma, r.theory), xytext=(4, -8), textcoords="offset points",
                        fontsize=7, color=r.color)
    lim = 1.1 * max(df.sigma.max(), df.theory.max())
    ax.plot([0, lim], [0, lim], "k--", alpha=0.8)
    ax.set_xlim(-0.005, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.set_xlabel(r"Growth rate from simulations, $\sigma/f^*$")
    ax.set_ylabel(r"Growth rate from stability analysis, $\sigma/f^*$")
    ax.set_title(rf"$r^2={r2:.2f}$" + (rf" (PSI only: {r2_psi:.2f})" if include_all else ""))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{prefix}.pdf")
    fig.savefig(f"{prefix}.png", dpi=200)
    print(df[["sim", "S", "Ri_inv", "delta", "sigma", "theory", "no_growth"]].round(4).to_string(index=False))
    print(f"r^2 = {r2:.3f} over {len(df)} runs (PSI only: {r2_psi:.3f} over {len(psi)}); wrote {prefix}.pdf/.png/.csv")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--all"]
    main(*args, include_all="--all" in sys.argv[1:])
