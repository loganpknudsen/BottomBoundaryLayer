"""Figure 6: growth rate over omega_min^2/(f*)^2 - delta for S_inf = 0.5, 1, 1.5.

Source: stability_analysis_output_reader.ipynb (origin/Base_Case) cell 31, which does not run as
committed (it reads the S-gamma dataset and undefined boundary arrays); here it reads freq.nc and
draws the PSI boundary the same way as Figs 3-5. Masked with
stability_common.mask (currently: above the convective limit for each delta,
or non-PSI below gamma_l (omega_min^2/f*^2 > 0.25), or sigma/f < 1e-3). Colour levels are fixed
(0 to 0.3, as in Fig 4) rather than spanning the data range, so unmasked SI growth rates saturate
instead of compressing the colour scale.
"""
import matplotlib.pyplot as plt
import numpy as np

import stability_common as sc
from paths import save

data = sc.load("freq")
vmax = 0.3
levels = np.linspace(0, vmax, 100)

fig, axs = plt.subplots(1, 3, figsize=(8, 3), sharey=True, layout="constrained")
for k, (ax, label) in enumerate(zip(axs, ("(a.)", "(b.)", "(c.)"))):
    d = data.isel(slope_burger_number=k)
    S = float(d.slope_burger_number)
    x = d.omega_ratio.values
    y = np.broadcast_to(d.delta.values[:, None], x.shape)
    cf = ax.contourf(x, y, sc.mask(d.growth_rate.values * np.sqrt(1 + S**2), d.strat_values.values, S, y, d.growth_rate.values,
                             d.frequency.values),
                     levels=levels, extend="both", cmap="YlOrRd")
    ax.contour(x, y, sc.psi_unstable(d.frequency.values, d.growth_rate.values, S, d.strat_values.values, y),
               levels=[sc.BOUNDARY_LEVEL], colors="k", linewidths=1)
    ax.text(-0.075, 1.05, label)
    ax.set_xlabel(r"$\omega_{\mathrm{min}}^2/(f^*)^2$")
axs[0].set_xlim([-0.05, 0.45])
axs[0].set_ylim([0, 1])
axs[0].set_ylabel(r"$\delta$")
cbar = fig.colorbar(cf, ax=axs[2], ticks=np.arange(0, vmax + 0.01, 0.05))
cbar.set_label(r"Growth Rate ($\sigma /f$)")
save(fig, __file__)
