"""Figure 2: hodographs of the near-inertial oscillation, (a) u-v and (b) u-b, for S = 0.5, 1, 1.5.

Source: oscillations_mangitudes.ipynb (origin/Base_Case), cells 1, 3-5 and 10.
No data needed: the oscillation is analytic, evaluated at gamma = gamma_u, delta = 0.8,
phi = -pi/2, with velocities normalised by delta*Lambda and buoyancy by delta*Lambda*N.
"""
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paths import save  # noqa: E402


def oscillation(t, S, gamma_set):
    theta = 0.1
    f = 1e-4
    N = S * f / np.tan(theta)
    if gamma_set == "gamma_u":
        gm = 1 / (1 + S**2)
    elif gamma_set == "gamma_m":
        gm = 0.5 * (1 / (1 + S**2) + (3 - S**2) / (3 * (1 + S**2)))
    elif gamma_set == "gamma_l":
        gm = (3 - S**2) / (3 * (1 + S**2))
    phi = -np.pi / 2
    beta = (1 + S**2) ** (0.5)
    fstar = f * beta
    u = np.cos(t + phi)
    v = -f * np.sin(t + phi) / (fstar)
    b = (-N * theta * np.sin(t + phi)) / (fstar)
    return [u, v, b]


t = np.linspace(0, 2 * np.pi, 100)
u_S05_u, v_S05_u, b_S05_u = oscillation(t, 0.5, "gamma_u")
u_S1_u, v_S1_u, b_S1_u = oscillation(t, 1, "gamma_u")
u_S15_u, v_S15_u, b_S15_u = oscillation(t, 1.5, "gamma_u")

fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(8, 16 / 3))

axes[0].plot(u_S05_u, v_S05_u, color="k", zorder=0, label=r"$S_\infty =0.5$")
axes[0].scatter(u_S05_u[0], v_S05_u[0], zorder=1, color="k", label="Initial Condtions")
axes[0].plot(u_S1_u, v_S1_u, color="b", zorder=0, label=r"$S_\infty =1$")
axes[0].scatter(u_S1_u[0], v_S1_u[0], zorder=1, color="k")
axes[0].plot(u_S15_u, v_S15_u, color="g", zorder=0, label=r"$S_\infty =1.5$")
axes[0].scatter(u_S15_u[0], v_S15_u[0], zorder=1, color="k")
axes[0].set_aspect("equal")
axes[0].text(-1.1, 1.05, "(a.)", fontsize=12)
axes[0].set_yticks([-1, -0.5, 0, 0.5, 1])
axes[0].set_ylim([-1.01, 1.01])
axes[0].set_xlabel(r'$\breve{u}/(\delta\Lambda)$', fontsize=12)
axes[0].set_ylabel(r'$\breve{v}/(\delta\Lambda)$', fontsize=12)

axes[1].plot(u_S05_u, b_S05_u, zorder=0, color="k")
axes[1].scatter(u_S05_u[0], b_S05_u[0], zorder=1, color="k")
axes[1].plot(u_S1_u, b_S1_u, zorder=0, color="b")
axes[1].scatter(u_S1_u[0], b_S1_u[0], zorder=1, color="k")
axes[1].plot(u_S15_u, b_S15_u, zorder=0, color="g")
axes[1].scatter(u_S15_u[0], b_S15_u[0], zorder=1, color="k")
axes[1].set_xlabel(r'$\breve{u}/(\delta\Lambda)$', fontsize=12)
axes[1].set_ylabel(r'$\breve{b}/(\delta\Lambda N_\infty)$', fontsize=12)
axes[1].text(-1.1, 1.05, "(b.)", fontsize=12)
axes[1].set_aspect('equal')

fig.legend(bbox_to_anchor=(0.55, 0.05), loc="lower center", ncol=4)
plt.tight_layout()
save(fig, __file__)
