"""Fit the exponential-phase growth rate of each run, with an automatic window.

Window rule (agreed 2026-10-14, replacing the earlier p rule):
  * sigma_1P(T) is the one-period growth rate from diagnose_growth.py; T_peak = argmax sigma_1P,
    sigma_peak its value there.
  * The window spans n e-folds of perturbation amplitude at the simulation's own peak rate,
    L = n / (2 pi sigma_peak) periods (sigma in units of f_s, T in periods), rounded to a whole
    number of periods (at least 1) so the intra-period modulation cancels, with n = N_EFOLDS.
  * The window ends at T_peak. sigma_1P(T_peak) describes [T_peak, T_peak + 1], so the peak
    period itself, which can carry a short spike as the growth phase ends (Sim5), is excluded.
  * A run whose window would start before the first output is flagged, not fitted. This is how a
    run without growth (sigma_peak ~ 0, so L is very large; Sim6) is caught.
  * sigma_sim is half the least-squares slope of ln k over the window, in units of f_s.

The earlier p rule (window where sigma_1P >= (1 - p) max sigma_1P) failed for Sim5, whose spike
at the end of the growth phase sets the maximum above the shelf. Every run is also fitted for
each n in N_VALUES so the sensitivity can be judged; n = 3 is too short for the fastest runs
(2 periods, spike-dominated) and n = 5 starts on the initial rise in Sim14.

Usage: python fit_growth.py DATA_DIR [OUT_PREFIX]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats

from diagnose_growth import F, THEORY, THETA, load_lnk, one_period_growth
from sims import SIMS, tke_filename

N_EFOLDS = 4
N_VALUES = (3, 4, 5)


def find_window(T, Tg, sig, n):
    """Return (T_start, L) for the fit, or (None, reason) if the rule is not met."""
    i = int(np.argmax(sig))
    if sig[i] <= 0:
        return None, "no growth"
    L = max(1, int(round(n / (2 * np.pi * sig[i]))))
    if Tg[i] - L < T[0]:
        return None, f"window {L} periods starts before T = 0"
    return Tg[i] - L, L


def fit(T, lnk, T0, L, S):
    sel = (T >= T0) & (T <= T0 + L)
    lr = scipy.stats.linregress(T[sel], lnk[sel])
    # ln k = 2 sigma t, and T = f_s t / 2pi, so sigma / f_s = slope / (4 pi)
    return lr.slope / (4 * np.pi), lr.rvalue**2


def main(data_dir, prefix="growth_fits"):
    data_dir = Path(data_dir)
    rows, curves = [], {}
    for name, s in SIMS.items():
        path = data_dir / tke_filename(s)
        if not path.exists():
            print(f"{name}: missing {path.name}")
            continue
        T, lnk = load_lnk(path, s["S"])
        Tg, sig = one_period_growth(T, lnk)
        curves[name] = (Tg, sig)
        for n in N_VALUES:
            T0, L = find_window(T, Tg, sig, n)
            row = {"sim": name, "S": s["S"], "Ri_inv": s["Ri_inv"], "delta": s["delta"],
                   "si_unstable": s["si_unstable"], "n": n, "theory": THEORY[name]}
            if T0 is None:
                row.update(status=L)
            else:
                sigma, r2 = fit(T, lnk, T0, L, s["S"])
                row.update(status="ok", T_start=T0, T_end=T0 + L, periods=L, sigma=sigma, r2=r2,
                           ratio=sigma / THEORY[name])
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(f"{prefix}.csv", index=False)

    wide = df.pivot(index="sim", columns="n", values="sigma").reindex(list(curves))
    wide.columns = [f"n={n}" for n in wide.columns]
    wide["spread"] = wide.max(axis=1) - wide.min(axis=1)
    wide["theory"] = [THEORY[n] for n in wide.index]
    wide["SI"] = [SIMS[n]["si_unstable"] for n in wide.index]
    print(wide.round(4).to_string())
    flagged = df[df.status != "ok"]
    if len(flagged):
        print("\nnot fitted:")
        print(flagged[["sim", "n", "status"]].to_string(index=False))
    windows = df[df.status == "ok"].pivot(index="sim", columns="n", values=["T_start", "T_end"])
    print("\nwindows (T_start, T_end):")
    print(windows.round(2).to_string())

    # Diagnostic: sigma_1P with each n's window shaded
    names = list(curves)
    ncol = 3
    nrow = -(-len(names) // ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.6 * nrow), squeeze=False)
    colors = dict(zip(N_VALUES, ("tab:blue", "tab:orange", "tab:green")))
    for ax, name in zip(axs.flat, names):
        Tg, sig = curves[name]
        ax.plot(Tg, sig, "k", lw=0.8)
        ax.axhline(THEORY[name], color="r", ls="--", lw=1)
        ax.axhline(0, color="0.6", lw=0.5)
        sub = df[(df.sim == name) & (df.status == "ok")]
        for k, (_, r) in enumerate(sub.iterrows()):
            ax.axvspan(r.T_start, r.T_end, ymin=k / 3, ymax=(k + 1) / 3, color=colors[r.n], alpha=0.3,
                       label=f"n={r.n:.0f} ({r.periods:.0f} P): {r.sigma:.3f}")
        tag = " (SI)" if SIMS[name]["si_unstable"] else ""
        ax.set_title(f"{name}{tag}", fontsize=9)
        ax.set_ylim(-0.1, max(0.35, 1.2 * THEORY[name]))
        ax.set_xlim(0, Tg[-1] + 1)
        if len(sub):
            ax.legend(fontsize=6, loc="upper right")
    for ax in axs.flat[len(names):]:
        ax.set_visible(False)
    fig.supxlabel(r"$T = f^* t/2\pi$")
    fig.supylabel(r"$\sigma_{1P}/f^*$")
    fig.tight_layout()
    fig.savefig(f"{prefix}_windows.png", dpi=110)
    print(f"\nwrote {prefix}.csv and {prefix}_windows.png")


if __name__ == "__main__":
    main(*sys.argv[1:])
