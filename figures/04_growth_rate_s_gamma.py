"""Figure 4: growth rate over S_inf-gamma for (a) delta = 0.5 and (b) delta = 1.0.

Source: stability_analysis_output_reader.ipynb (origin/Base_Case) cell 27. Data: sgamma.nc.
Masked with stability_common.mask (currently: points above the convective limit).
Growth rate is plotted as sigma/f = (sigma/f*) sqrt(1 + S^2), as in the notebook.
"""
import matplotlib.pyplot as plt
import numpy as np

import stability_common as sc
from paths import save

data = sc.load("sgamma")
S = data.slope_burger_number.values
vmax = 0.3
levels = np.linspace(0, vmax + 1 / 20, 100)

fig, axs = plt.subplots(1, 2, figsize=(8, 3), sharey=True, layout="constrained")
for ax, delta, label in zip(axs, (0.5, 1.0), ("(a.)", "(b.)")):
    Sm, g, d = sc.sgamma_mesh(data, delta)
    cf = ax.contourf(Sm, g, sc.mask(d.growth_rate.values * np.sqrt(1 + Sm**2), g, Sm, delta), levels=levels, extend="both",
                     vmax=vmax, cmap="RdYlBu_r")
    ax.contour(Sm, g, sc.sgamma_boundary(d), levels=[sc.BOUNDARY_LEVEL], colors="k", linestyles="--", linewidths=1.4)
    sc.analytic_lines(ax, S, delta)
    ax.text(0.075, 1.05, label)
    ax.set_xlim([0.15, 1.9])
    ax.set_xticks(np.arange(0.2, 2, 0.2))
    ax.set_xlabel(r"$S_\infty$")
axs[0].set_ylim([0.05, 1])
axs[0].set_ylabel(r"$\gamma$")
cbar = fig.colorbar(cf, ax=axs[1], ticks=np.arange(0, vmax + 0.06, 0.05))
cbar.set_label(r"Growth Rate ($\sigma/f$)")
save(fig, __file__)
