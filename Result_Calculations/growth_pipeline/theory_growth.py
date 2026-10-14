"""Floquet growth rate for each simulation, evaluated at its exact (S, gamma, delta).

Uses the small-angle PSI system from stability_analysis_full_run_small_angle.py (time
scaled by f*, so sigma is in units of f*). Instead of looking up the nearest point on
a precomputed grid, the Floquet exponent is maximised over the perturbation slope m
directly: a coarse sweep followed by a bounded scalar refinement.
"""
import numpy as np
import scipy.integrate as sc
import scipy.optimize as so

from sims import SIMS, gamma

THETA = 0.1
F = 1e-4
TAU = 2 * np.pi


def floquet(m, S, gm, delta):
    """Return the Floquet exponent (sigma, frequency) with the largest real part, in units of f*."""
    S2 = S**2
    N2 = S2 * F**2 / THETA**2
    fstar = F * np.sqrt(1 + S2)
    lmbd = N2 * THETA * gm / F
    A0 = -S2 / (THETA * (1 + S2)) * (m + (1 - gm) / THETA)
    A1 = m + S2 / (THETA * (1 + S2))

    def rhs(t, v):
        psi, D = v
        mt = m + delta * lmbd / fstar * np.sin(t + np.pi / 2)
        return [(A0 - A1 * mt) * D, psi / (1 + mt**2)]

    cols = [sc.solve_ivp(rhs, (0, TAU), ic, rtol=1e-10, atol=1e-12).y[:, -1] for ic in ([1, 0], [0, 1])]
    eigs = np.log(np.linalg.eigvals(np.array(cols).T).astype(complex)) / TAU
    e = eigs[np.argmax(eigs.real)]
    return e.real, abs(e.imag)


def max_growth(S, gm, delta):
    ms = np.arange(-40, 5, 0.1)
    sig = np.array([floquet(m, S, gm, delta)[0] for m in ms])
    i = np.argmax(sig)
    res = so.minimize_scalar(lambda m: -floquet(m, S, gm, delta)[0],
                             bounds=(ms[max(i - 1, 0)], ms[min(i + 1, len(ms) - 1)]),
                             method="bounded", options={"xatol": 1e-6})
    sigma, freq = floquet(res.x, S, gm, delta)
    return sigma, freq, res.x


if __name__ == "__main__":
    print(f"{'sim':<6}{'S':>5}{'Ri_inv':>8}{'delta':>7}{'gamma':>8}{'sigma/f*':>10}{'Im/f*':>8}{'m_max':>8}")
    for name, s in SIMS.items():
        g = gamma(s["S"], s["Ri_inv"])
        sigma, freq, m = max_growth(s["S"], g, s["delta"])
        print(f"{name:<6}{s['S']:>5}{s['Ri_inv']:>8}{s['delta']:>7}{g:>8.4f}{sigma:>10.4f}{freq:>8.3f}{m:>8.2f}")
