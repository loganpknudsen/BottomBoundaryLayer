"""Figure 9: contribution of GSP + BP to the TKE change over one period for (a) delta = 0.5 and
(b) delta = 1.0.

Source: stability_analysis_output_reader.ipynb (origin/Base_Case) cell 89; the budget now comes
from sgamma.nc instead of notebook cells 68-69. The colourbar is labelled (unlabelled in the notebook). Masked with stability_common.mask (currently:
above the convective limit, or below gamma_l with sigma/f < 1e-3 or a non-PSI fastest mode).
"""
import matplotlib.pyplot as plt
import numpy as np

import stability_common as sc
from paths import save

data = sc.load("sgamma")
S = data.slope_burger_number.values
vmax = 1.5
levels = np.linspace(-vmax, vmax + 1 / 100, 100)

fig, axs = plt.subplots(1, 2, figsize=(8, 3), sharey=True, layout="constrained")
for ax, delta, label in zip(axs, (0.5, 1.0), ("(a.)", "(b.)")):
    Sm, g, d = sc.sgamma_mesh(data, delta)
    cf = ax.contourf(Sm, g, sc.mask(sc.budget(d, "GSP") + sc.budget(d, "BP"), g, Sm, delta, d.growth_rate.values, d.frequency.values), levels=levels, vmax=vmax, vmin=-vmax,
                     extend="both", cmap="bwr")
    ax.contour(Sm, g, sc.sgamma_boundary(d), levels=[sc.BOUNDARY_LEVEL], colors="k", linestyles="--", linewidths=1.4)
    sc.analytic_lines(ax, S, delta)
    ax.text(0.1, 1.05, label)
    ax.set_xlim([0.15, 1.9])
    ax.set_xticks(np.arange(0.2, 2, 0.2))
    ax.set_xlabel(r"$S_\infty$")
axs[0].set_ylim([0.05, 1])
for ax in axs:
    sc.label_limits(ax)
axs[0].set_ylabel(r"$\gamma$")
cbar = fig.colorbar(cf, ax=axs, ticks=np.arange(-vmax, vmax + 0.06, 0.5))
cbar.set_label(r"Fraction of TKE change over one period")
save(fig, __file__)
