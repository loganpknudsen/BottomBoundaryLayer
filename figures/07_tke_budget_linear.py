"""Figure 7: TKE budget of the linear (Floquet) solution for S = 1, delta = 0.8, gamma = gamma_u.
(a) tendency terms divided by k over 4 inertial periods; (b) time-integrated terms.

Source: Floquet_Analysis_Basis_Diagram_Hilditch.ipynb (origin/Base_Case), cells 3-7 and 36.
No data needed. The equations are kept exactly as in the notebook (small-angle Eq 2.28 system,
the optimal slope found by scanning m in [-200, 200]); this is a faithful extraction, except that
the along-slope ASP term no longer carries the notebook's extra cos(theta) (see Eq 2.31).

With --residual, also plots the budget residual in grey (panel a: (dk/dt - GSP - ASP - BP)/k;
panel b: k(t) - k(0) - integral of GSP + ASP + BP) and saves as 07_tke_budget_linear_residual.
"""
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as sc

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import save  # noqa: E402


def m_t(t, m, lmbd, delta, fstar, phi):
    return m + delta * lmbd / (fstar) * np.sin(t + phi)


def PSI_system(v, t, m, S2, theta, gm, lmbd, delta, fstar, phi):
    psi, D = v
    A0 = -S2 / (theta * (1 + S2)) * (m + (1 - gm) / theta)
    A1 = m + (S2) * (theta * (1 + S2)) ** (-1)
    dpsidt = (A0 - A1 * m_t(t, m, lmbd, delta, fstar, phi)) * D
    dDdt = (1 + m_t(t, m, lmbd, delta, fstar, phi) ** 2) ** (-1) * psi
    return [dpsidt, dDdt]


def parameters(S, gamma_set):
    theta = 0.1
    f = 1e-4
    N2 = S**2 * f**2 / theta**2
    phi = np.pi / 2
    S2 = N2 * theta**2 / f**2
    beta = (1 + S2) ** (0.5)
    fstar = f * beta
    if gamma_set == "gamma_u":
        gm = 1 / (1 + S2)
    elif gamma_set == "gamma_m":
        gm = 0.5 * (1 / (1 + S2) + (3 - S2) / (3 * (1 + S2)))
    elif gamma_set == "gamma_l":
        gm = (3 - S2) / (3 * (1 + S2))
    delta = 0.8
    lmbd = N2 * theta * gm / f
    return theta, f, N2, phi, S2, beta, fstar, gm, delta, lmbd


def floquet_m_max(S, gamma_set):
    theta, f, N2, phi, S2, beta, fstar, gm, delta, lmbd = parameters(S, gamma_set)
    m = np.arange(-200, 200 + 0.05, 0.05)
    growth_rates = []
    tau = 2 * np.pi
    dt = 100
    t = np.linspace(0, tau + 1 / dt, dt)
    for i in m:
        sol1 = sc.odeint(PSI_system, [1, 0], t, args=(i, S2, theta, gm, lmbd, delta, fstar, phi))
        sol2 = sc.odeint(PSI_system, [0, 1], t, args=(i, S2, theta, gm, lmbd, delta, fstar, phi))
        M = np.array([[sol1[-1, 0], sol1[-1, 1]], [sol2[-1, 0], sol2[-1, 1]]])
        eigs = np.log(np.linalg.eig(M)[0] + 0 * 1j) / (tau)
        growth_rates.append(eigs.real)
    growth_rates = np.array(growth_rates)
    return m[np.argmax(np.abs(growth_rates[:, 0]))]


def TKE_and_velocities(S, gamma_set, m_max, tf):
    theta, f, N2, phi, S2, beta, fstar, gm, delta, lmbd = parameters(S, gamma_set)
    n = fstar / lmbd
    Ri = N2 * (1 - gm) / lmbd**2
    alpha = (1 - gm) / gm
    dt = 0.01
    t = np.arange(0, 2 * np.pi * tf + dt, dt)
    output = sc.odeint(PSI_system, [0, 1], t, args=(m_max, S2, theta, gm, lmbd, delta, fstar, phi))
    w = ((m_t(t, m_max, lmbd, delta, fstar, phi)) ** 2 + 1) ** (-1) * output[:, 0]
    u = m_t(t, m_max, lmbd, delta, fstar, phi) * w * n
    v = -f * (gm * S2 / theta - theta + m_max) * output[:, 1] / lmbd
    b = -theta * (m_max + (1 - gm) / theta) * output[:, 1] / (1 - gm)
    ke = 0.5 * (u**2 + v**2 + n**2 * w**2)
    dkdt = np.gradient(ke, t[-1] - t[-2])
    GSP = -v * w
    # Eq 2.31 along-slope ASP term: delta (f/f*) sin(t + phi) v w (the notebook had an extra cos(theta),
    # left over from the simulations' f* = f cos(theta) sqrt(1 + S^2) convention)
    WSP = -u * w * delta * np.cos(t + np.pi / 2) + v * w * (f * delta / fstar * np.sin(t + np.pi / 2))
    BFLUX = alpha * beta ** (-1) * u * b + w * b * Ri
    intGSP = sc.cumulative_trapezoid(GSP[:], t[:])
    intWSP = sc.cumulative_trapezoid(WSP[:], t[:])
    intBFLUX = sc.cumulative_trapezoid(BFLUX[:], t[:])
    return t, [dkdt, GSP, WSP, BFLUX], [ke, intGSP, intWSP, intBFLUX]


SHOW_RESIDUAL = "--residual" in sys.argv[1:]

m_max_S1_u = floquet_m_max(1, "gamma_u")
time_S1_u, TKE_tends_S1_u, TKE_int_S1_u = TKE_and_velocities(1, "gamma_u", m_max_S1_u, 6)

fig, axes = plt.subplots(2, 1, sharex=False, figsize=(9, 9))

ti = 1256
tf = ti * 3
T = (time_S1_u[ti:tf] - time_S1_u[ti]) / (2 * np.pi)
k = TKE_int_S1_u[0][ti:tf]
axes[0].plot(T, TKE_tends_S1_u[0][ti:tf] / k, color="k", label=r"$\dfrac{1}{\hat{k}}\dfrac{\partial \hat{k}}{\partial \hat{t}}$")
axes[0].plot(T, TKE_tends_S1_u[1][ti:tf] / k, color="g", label=r"$GSP/\hat{k}$", alpha=0.9)
axes[0].plot(T, TKE_tends_S1_u[2][ti:tf] / k, color="r", label=r"$ASP/\hat{k}$")
# TODO(decision): the preprint Fig 7 labels this term BFLUX; the Base_Case notebook (used here)
# renamed it BP, matching the text's Eq 2.31. Same data either way.
axes[0].plot(T, TKE_tends_S1_u[3][ti:tf] / k, color="blue", linestyle="-", label=r"$BP/\hat{k}$")
if SHOW_RESIDUAL:
    resid = TKE_tends_S1_u[0] - TKE_tends_S1_u[1] - TKE_tends_S1_u[2] - TKE_tends_S1_u[3]
    axes[0].plot(T, resid[ti:tf] / k, color="grey", label=r"Residual$/\hat{k}$")
axes[0].hlines(0, 4, 0, color="k", linestyle="--")
axes[0].set_ylabel(r"TKE Tendency Term", fontsize=12)
axes[0].legend(loc="lower left", ncol=4)
axes[0].text(-0.375, 0.58, "(a.)", fontsize=11)
axes[0].set_xlim([0, 4])
axes[0].set_ylim([-0.8, 0.7])

axes[1].plot((time_S1_u) / (2 * np.pi), TKE_int_S1_u[0], color="k", label=r"$\hat{k}$")
axes[1].plot((time_S1_u[:-1]) / (2 * np.pi), TKE_int_S1_u[1], color="g", label=r"$\int_0^{\hat{t}}GSP\,d\hat{t}$", alpha=0.9)
axes[1].plot((time_S1_u[:-1]) / (2 * np.pi), TKE_int_S1_u[2], color="r", label=r"$\int_0^{\hat{t}}ASP\,d\hat{t}$")
axes[1].plot((time_S1_u[:-1]) / (2 * np.pi), TKE_int_S1_u[3], color="blue", linestyle="-", label=r"$\int_0^{\hat{t}}BP\,d\hat{t}$")
if SHOW_RESIDUAL:
    k_all = TKE_int_S1_u[0]
    resid_int = (k_all[1:] - k_all[0]) - (TKE_int_S1_u[1] + TKE_int_S1_u[2] + TKE_int_S1_u[3])
    axes[1].plot((time_S1_u[:-1]) / (2 * np.pi), resid_int, color="grey", label="Residual")
    print(f"max |residual| / max |dk/dt| (panel a range): "
          f"{np.max(np.abs(resid[ti:tf])) / np.max(np.abs(TKE_tends_S1_u[0][ti:tf])):.3f}")
    print(f"integrated residual / (k(T) - k(0)) at T = 6: {resid_int[-1] / (k_all[-1] - k_all[0]):.4f}")
axes[1].set_xlabel(r"Modified Inertial Period", fontsize=12)
axes[1].set_ylabel(r"$\hat{k}$", fontsize=14)
axes[1].legend(loc="lower left", ncol=4)
axes[1].text(-0.59, 0.5e5, "(b.)", fontsize=11)
axes[1].set_xlim([0, 6])
axes[1].set_yscale("symlog", linthresh=1)

save(fig, __file__ if not SHOW_RESIDUAL else str(pathlib.Path(__file__).with_name("07_tke_budget_linear_residual.py")))
