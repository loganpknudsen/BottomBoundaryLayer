"""Figure 5: perturbation slope mu_o / mu_B over S_inf-gamma for (a) delta = 0.5, (b) delta = 1.0.

Source: stability_analysis_output_reader.ipynb (origin/Base_Case) cell 28. Data: sgamma.nc.
mu_B = -(1 - gamma)/theta is the isopycnal-aligned slope. Masked with stability_common.mask
(currently: above the convective limit, weak growth below
gamma_l, or non-PSI below gamma_l). At stable points the "optimal" slope is the
argmax of a zero growth rate and is not meaningful.
"""
import matplotlib.pyplot as plt
import numpy as np

import stability_common as sc
from paths import save

THETA = 0.1
data = sc.load("sgamma")
S = data.slope_burger_number.values
vmin, vmax = 0.25, 1.75
levels = np.linspace(vmin, vmax + 1 / 20, 100)

fig, axs = plt.subplots(1, 2, figsize=(8, 3), sharey=True, layout="constrained")
for ax, delta, label in zip(axs, (0.5, 1.0), ("(a.)", "(b.)")):
    Sm, g, d = sc.sgamma_mesh(data, delta)
    mu_B = -(1 - g) / THETA
    cf = ax.contourf(Sm, g, sc.mask(d.slope_angle.values / mu_B, g, Sm, delta, d.growth_rate.values, d.frequency.values), extend="both", levels=levels, vmax=vmax, vmin=vmin,
                     cmap="RdYlBu")
    ax.contour(Sm, g, sc.sgamma_boundary(d), levels=[sc.BOUNDARY_LEVEL], colors="k", linestyles="--", linewidths=1.4)
    sc.analytic_lines(ax, S, delta)
    ax.text(0.075, 1.05, label)
    ax.set_xlim([0.15, 1.9])
    ax.set_xticks(np.arange(0.2, 2, 0.2))
    ax.set_xlabel(r"$S_\infty$")
axs[0].set_ylim([0.05, 1])
axs[0].set_ylabel(r"$\gamma$")
cbar = fig.colorbar(cf, ax=axs[1], ticks=np.linspace(0.25, 1.75, 7))
cbar.set_label(r"$\mu_o/\mu_B$")
save(fig, __file__)
