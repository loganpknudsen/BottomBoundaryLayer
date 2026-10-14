"""Refit TKE growth rates exactly as Result_Calculations/TKE_terms.ipynb does,
comparing the manuscript window (ts=1500, tf=3500) against the current one (ts=0, tf=-1).

Usage: python refit_growth.py FILE [FILE ...]
"""
import re
import sys

import numpy as np
import pandas as pd
import scipy as sc
import scipy.stats
import xarray as xr

theta = 0.1
f = 1e-4
WINDOWS = {"march (1500:3500)": (1500, 3500), "current (0:-1)": (0, -1)}

print(f"{'file':<62} {'window':<18} {'IP range':>13} {'growth':>9} {'r^2':>6}")
for path in sys.argv[1:]:
    name = path.split("/")[-1]
    S = float(re.search(r"Sinf_([0-9.]+?)_", name).group(1))
    fstar = f * np.cos(theta) * (1 + S**2) ** 0.5

    k = xr.open_dataset(path).k.squeeze()
    if k.ndim != 1:
        sys.exit(f"{name}: k has dims {k.dims}, expected a time series")
    seconds = k.time / pd.Timedelta("1s")
    ip = seconds * fstar / (2 * np.pi)

    for label, (ts, tf) in WINDOWS.items():
        if ts >= k.size:
            print(f"{name:<62} {label:<18} {'too short':>13}")
            continue
        y = np.log(k.values[ts:tf] / k.values[ts])
        lr = sc.stats.linregress(x=seconds[ts:tf], y=y)
        growth = np.round(lr.slope / (2 * f), 7)  # same normalization as the notebook
        span = f"{float(ip[ts]):.1f}-{float(ip[ts:tf][-1]):.1f}"
        print(f"{name:<62} {label:<18} {span:>13} {growth:>9.4f} {lr.rvalue**2:>6.3f}")
