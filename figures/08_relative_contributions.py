"""Figure 8: relative contributions of (a) ASP, (b) GSP and (c) BP to the TKE change over one
period, for delta = 0.5.

Source: stability_analysis_output_reader.ipynb (origin/Base_Case) cell 87; the budget now comes
from sgamma.nc (computed with the growth-rate equations, see linear_theory/floquet.py) instead
of notebook cell 68. Masked with stability_common.mask (currently: above the convective
limit, weak growth below gamma_l, or non-PSI below gamma_l). Change from the notebook: the colourbar was labelled as a growth
rate; the plotted quantity is a dimensionless fraction of the TKE change.
"""
import matplotlib.pyplot as plt
import numpy as np

import stability_common as sc
from paths import save

DELTA = 0.5
data = sc.load("sgamma")
S = data.slope_burger_number.values
vmax = 1.5
levels = np.linspace(-vmax, vmax + 1 / 100, 100)
Sm, g, d = sc.sgamma_mesh(data, DELTA)
boundary = sc.sgamma_boundary(d)

fig, axs = plt.subplots(1, 3, figsize=(8, 3), sharey=True, layout="constrained")
for ax, term, label in zip(axs, ("ASP", "GSP", "BP"), ("(a.)", "(b.)", "(c.)")):
    cf = ax.contourf(Sm, g, sc.mask(sc.budget(d, term), g, Sm, DELTA, d.growth_rate.values, d.frequency.values), levels=levels, vmax=vmax, vmin=-vmax, extend="both", cmap="bwr")
    ax.contour(Sm, g, boundary, levels=[sc.BOUNDARY_LEVEL], colors="k", linestyles="--", linewidths=1.4)
    sc.analytic_lines(ax, S, DELTA)
    ax.text(0.05, 1.05, label)
    ax.set_xlim([0.15, 1.9])
    ax.set_ylim([0.05, 1])
    ax.set_xticks(np.arange(0.2, 2, 0.4))
    ax.set_xlabel(r"$S_\infty$")
axs[0].set_ylabel(r"$\gamma$")
cbar = fig.colorbar(cf, ax=axs[2], ticks=np.arange(-vmax, vmax + 0.06, 0.5))
cbar.set_label(r"Fraction of TKE change over one period")
save(fig, __file__)
