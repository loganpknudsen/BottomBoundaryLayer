"""Trial fit-window rule: a window of n amplitude e-folds ending at the peak growth rate.

Alternative to the p rule in fit_growth.py, being tested (2026-10-14) because a short spike in
sigma_1P at the end of Sim5's growth phase sets max(sigma_1P) above the shelf, so the p = 5% and
10% windows fail there.

  * sigma_1P(T) is the one-period growth rate from diagnose_growth.py; T_peak = argmax sigma_1P,
    sigma_peak its value there.
  * The window length is n e-folds of perturbation amplitude at the simulation's own peak rate:
    L = n / (2 pi sigma_peak) periods (sigma in units of f_s, T in periods), rounded to a whole
    number of periods (at least 1) so the intra-period modulation cancels.
  * sigma_1P(T) describes [T, T+1], so the window can end at T_peak (excluding the peak period,
    which may already be nonlinear) or at T_peak + 1 (including it). Both are fitted.
  * A run whose window would start before the first output is flagged, not fitted. This is how
    a run without growth (sigma_peak ~ 0, so L is very large) is caught; no extra threshold.
  * sigma_sim is half the least-squares slope of ln k over the window, as in fit_growth.py.

Every run is fitted for each n in N_VALUES and both end choices, next to the p = 10% fit.

Usage: python fit_growth_efold.py DATA_DIR [OUT_PREFIX]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from diagnose_growth import THEORY, load_lnk, one_period_growth
from fit_growth import find_window, fit
from sims import SIMS, tke_filename

N_VALUES = (3, 4, 5)
ENDS = {"peak": 0, "peak+1": 1}   # window end = T_peak + offset
P_REF = 0.10


def efold_window(T, Tg, sig, n, offset):
    """Return (T_start, L, T_end) for the fit, or (None, reason, None) if the rule is not met."""
    i = int(np.argmax(sig))
    if sig[i] <= 0:
        return None, "no growth", None
    L = max(1, int(round(n / (2 * np.pi * sig[i]))))
    T_end = Tg[i] + offset
    if T_end - L < T[0]:
        return None, f"window {L} periods starts before T = 0", None
    return T_end - L, L, T_end


def main(data_dir, prefix="growth_fits_efold"):
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
        base = {"sim": name, "S": s["S"], "Ri_inv": s["Ri_inv"], "delta": s["delta"],
                "si_unstable": s["si_unstable"], "theory": THEORY[name],
                "T_peak": Tg[np.argmax(sig)], "sigma_peak": sig.max()}

        T0, L = find_window(Tg, sig, P_REF)
        row = dict(base, method=f"p={P_REF:.0%}")
        if T0 is None:
            row.update(status=L)
        else:
            sigma, r2 = fit(T, lnk, T0, L, s["S"])
            row.update(status="ok", T_start=T0, T_end=T0 + L, periods=L, sigma=sigma, r2=r2,
                       ratio=sigma / THEORY[name])
        rows.append(row)

        for end, offset in ENDS.items():
            for n in N_VALUES:
                T0, L, T1 = efold_window(T, Tg, sig, n, offset)
                row = dict(base, method=f"n={n}, end={end}", n=n, end=end)
                if T0 is None:
                    row.update(status=L)
                else:
                    sigma, r2 = fit(T, lnk, T0, L, s["S"])
                    row.update(status="ok", T_start=T0, T_end=T1, periods=L, sigma=sigma, r2=r2,
                               ratio=sigma / THEORY[name])
                rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(f"{prefix}.csv", index=False)

    methods = list(dict.fromkeys(df.method))
    wide = df.pivot(index="sim", columns="method", values="sigma").reindex(index=list(curves),
                                                                           columns=methods)
    wide["theory"] = [THEORY[n] for n in wide.index]
    wide["SI"] = [SIMS[n]["si_unstable"] for n in wide.index]
    print("sigma / f*:")
    print(wide.round(3).to_string())
    print("\nperiods in window:")
    print(df.pivot(index="sim", columns="method", values="periods")
          .reindex(index=list(curves), columns=methods).to_string())
    flagged = df[df.status != "ok"]
    if len(flagged):
        print("\nnot fitted:")
        print(flagged[["sim", "method", "status"]].to_string(index=False))

    # Diagnostic: sigma_1P with each n's window shaded, one figure per end choice
    names = list(curves)
    ncol = 3
    nrow = -(-len(names) // ncol)
    colors = dict(zip(N_VALUES, ("tab:blue", "tab:orange", "tab:green")))
    for end in ENDS:
        fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.6 * nrow), squeeze=False)
        for ax, name in zip(axs.flat, names):
            Tg, sig = curves[name]
            ax.plot(Tg, sig, "k", lw=0.8)
            ax.axhline(THEORY[name], color="r", ls="--", lw=1)
            ax.axhline(0, color="0.6", lw=0.5)
            sub = df[(df.sim == name) & (df.end == end) & (df.status == "ok")]
            for k, (_, r) in enumerate(sub.iterrows()):
                ax.axvspan(r.T_start, r.T_end, ymin=k / 3, ymax=(k + 1) / 3, color=colors[r.n],
                           alpha=0.3, label=f"n={r.n:.0f} ({r.periods:.0f} P): {r.sigma:.3f}")
            tag = " (SI)" if SIMS[name]["si_unstable"] else ""
            ax.set_title(f"{name}{tag}", fontsize=9)
            ax.set_ylim(-0.1, max(0.35, 1.2 * THEORY[name]))
            ax.set_xlim(0, Tg[-1] + 1)
            if len(sub):
                ax.legend(fontsize=6, loc="upper right")
        for ax in axs.flat[len(names):]:
            ax.set_visible(False)
        fig.suptitle(f"window of n amplitude e-folds ending at T_peak{ENDS[end] and ' + 1' or ''}")
        fig.supxlabel(r"$T = f^* t/2\pi$")
        fig.supylabel(r"$\sigma_{1P}/f^*$")
        fig.tight_layout()
        fig.savefig(f"{prefix}_windows_end-{end.replace('+1', 'plus1')}.png", dpi=110)
    print(f"\nwrote {prefix}.csv and {prefix}_windows_end-*.png")


if __name__ == "__main__":
    main(*sys.argv[1:])
