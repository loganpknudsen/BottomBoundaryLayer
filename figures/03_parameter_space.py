"""Figure 3: PSI parameter space in S_inf-gamma, with the PSI boundary for delta = 0.5 and 1.

Source: stability_analysis_output_reader.ipynb (origin/Base_Case) cells 18, 19, 21. Data:
sgamma.nc. Change from the notebook: the boundary colours now match the legend (the notebook
drew the delta = 0.5 boundary in red under a "delta = 1.0" label).
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import stability_common as sc
from paths import save

data = sc.load("sgamma")
S = data.slope_burger_number.values

fig, ax = plt.subplots()
for delta, color in ((1.0, "red"), (0.5, "blue")):
    Sm, g, d = sc.sgamma_mesh(data, delta)
    ax.contour(Sm, g, sc.sgamma_boundary(d), levels=[sc.BOUNDARY_LEVEL], colors=color, linestyles="--", linewidths=1)

ax.fill_between(S, sc.gamma_lower_limit(S), sc.gamma_upper_limit(S), color="b", alpha=0.2)
ax.fill_between(S, sc.gamma_upper_limit(S), 1.0, color="k", alpha=0.15)
ax.plot(S, sc.gamma_upper_limit(S), c="k")
ax.plot(S, sc.gamma_lower_limit(S), c="k")
ax.set_ylabel(r"$\gamma$", fontsize=14)
ax.set_xlabel(r"$S_\infty$", fontsize=14)
ax.text(1.27, 0.75, "SI", fontsize=14)
ax.text(1.275, 0.25, "PSI", fontsize=14)
ax.text(0.25, 0.25, "Stable", fontsize=14)
ax.legend(handles=[Line2D([], [], color="black", linestyle="-", label=r"$\delta=0.0$"),
                   Line2D([], [], color="blue", linestyle="--", label=r"$\delta=0.5$"),
                   Line2D([], [], color="red", linestyle="--", label=r"$\delta=1.0$")], loc="lower left")
ax.set_xlim([0.15, 1.9])
ax.set_ylim([0.02, 1])
sc.label_limits(ax, fontsize=11)
save(fig, __file__)
