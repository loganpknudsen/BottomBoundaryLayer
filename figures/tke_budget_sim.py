"""Time-integrated TKE budget of a simulation (k, ASP = AGSP + WSP, GSP, BP, -epsilon), shared by
Figs 14 (Sim9) and 16 (Sim7); symmetric-log scale with linear threshold 1e-15.

Source: Result_Calculations/TKE_terms.ipynb, cells 3, 11 (time crop), 15-19, 28-34 and 37
(origin/Base_Case numbering; cell 33 in the submission commit 27f8414). Follows the submission
version, which matches the preprint figures: x-ticks every 2.5, y-ticks +-1e-6, +-1e-10, +-1e-14,
record cropped to outputs ts to tf + 999 (outputs every 0.005 periods, so ts = 1000 is T = 5; tf is
inferred from the preprint axes, not recorded in the code). Everything is measured from
t0 (T = 5): the curves are k - k(t0) and the integrals of each term from t0, labelled accordingly.

With --residual on the command line, also plots the budget residual
k - k(t0) - integral of (ASP + GSP + BP - epsilon) in grey (as in the Base_Case notebook, cell 37),
prints it as a fraction of k(T) - k(t0), and saves as <figure>_residual.
"""
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from paths import SIM_DIR, save

THETA = 0.1
F = 1e-4
TS = 1000   # first output of the integration (T = 5)


def plot(file, S, tf, script_file):
    show_residual = "--residual" in sys.argv[1:]
    ocean_data = xr.open_dataset(SIM_DIR / file).isel(time=slice(TS, tf + 1000))
    fstar = F * np.cos(THETA) * (1 + S**2) ** 0.5

    def with_coords(da):
        seconds = da.time / pd.Timedelta("1s")
        return da.assign_coords(intp=("time", (seconds * fstar / (2 * np.pi)).data), sec=("time", seconds.data))

    AGSP = with_coords(ocean_data.AGSP)
    WSP = with_coords(ocean_data.WSP)
    GSP = with_coords(ocean_data.GSP)
    BFLUX = with_coords(ocean_data.BFLUX)
    E = with_coords(ocean_data.E)
    integrated_AGSP = AGSP.isel(time=slice(0, -1)).cumulative_integrate(coord="sec")
    integrated_WSP = WSP.isel(time=slice(0, -1)).cumulative_integrate(coord="sec")
    integrated_GSP = GSP.isel(time=slice(0, -1)).cumulative_integrate(coord="sec")
    integrated_BFLUX = BFLUX.isel(time=slice(0, -1)).cumulative_integrate(coord="sec")
    integrated_E = (-1 * E.isel(time=slice(0, -1))).cumulative_integrate(coord="sec")
    KE_turb_ip = with_coords(ocean_data.k)

    start = 0
    end = -1
    fig = plt.figure(figsize=(7, 5))
    plt.plot(KE_turb_ip.intp[start:end], KE_turb_ip[start:end] - KE_turb_ip[start], label="$k - k(t_0)$", color="k")
    plt.plot(integrated_AGSP.intp[start:end], (integrated_AGSP[start:end] + integrated_WSP[start:end]) - (integrated_AGSP[start] + integrated_WSP[start]), label=r"$\int_{t_0}^tASP\,dt'$", color="red")
    plt.plot(integrated_GSP.intp[start:end], integrated_GSP[start:end] - integrated_GSP[start], label=r"$\int_{t_0}^tGSP\,dt'$", color="green")
    plt.plot(integrated_BFLUX.intp[start:end], integrated_BFLUX[start:end] - integrated_BFLUX[start], label=r"$\int_{t_0}^tBP\,dt'$", color="blue")
    plt.plot(integrated_E.intp[start:end], integrated_E[start:end] - integrated_E[start], label=r"$-\int_{t_0}^t\epsilon\,dt'$", color="pink")
    if show_residual:
        resid = KE_turb_ip - integrated_AGSP - integrated_GSP - integrated_WSP - integrated_BFLUX - integrated_E
        plt.plot(resid.intp[start:end], resid[start:end] - resid[start], label="Residual", color="grey", alpha=0.9)
        dk = float(KE_turb_ip[end] - KE_turb_ip[start])
        print(f"residual / (k(T) - k(t0)) at T = {float(resid.intp[end]):.2f}: {float(resid[end] - resid[start]) / dk:.4f}")
    plt.xlabel("Modified Inertial Period")
    plt.ylabel(r"k $(m^2 s^{-2})$")
    plt.yscale("symlog", linthresh=1e-15)
    plt.xticks([5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30])
    plt.yticks([1e-6, -1e-6, 1e-10, -1e-10, 1e-14, -1e-14])
    plt.legend(loc="lower left", prop={'size': 7.5}, ncol=2, handlelength=0.5)
    plt.xlim([KE_turb_ip.intp[start] - 0.05, KE_turb_ip.intp[end] + 0.01])

    out = pathlib.Path(script_file)
    save(fig, str(out) if not show_residual else str(out.with_name(out.stem + "_residual.py")))
