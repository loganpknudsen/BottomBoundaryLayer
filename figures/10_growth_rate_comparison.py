"""Figure 10: simulated vs. Floquet growth rates.

Every run refits the growth rates from the TKE files in BBL_SIM_DIR
(Result_Calculations/growth_pipeline/fit_growth.py), then plots them with make_figure.py, so the
figure never comes from a stale fits file. Fit window: n = 4 amplitude e-folds ending at the peak
growth rate. SI-unstable runs are included, in grey; runs with no exponential phase are plotted
at zero with an open marker.

Besides the figure, figures/output/ gets 10_growth_rate_comparison_fits.csv (all fits) and
10_growth_rate_comparison_fits_windows.png (the fit windows, for checking).
"""
import os
import sys

# The fit can read TKE files that running simulations are still writing; without this the reader
# takes an HDF5 lock that can make the simulation's next write fail. Must be set before HDF5 loads.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

from paths import OUT_DIR, REPO, SIM_DIR  # noqa: E402

sys.path.insert(0, str(REPO / "Result_Calculations" / "growth_pipeline"))
import fit_growth  # noqa: E402
import make_figure  # noqa: E402

OUT_DIR.mkdir(parents=True, exist_ok=True)
stem = str(OUT_DIR / "10_growth_rate_comparison")
print(f"fitting growth rates from {SIM_DIR}")
fit_growth.main(str(SIM_DIR), stem + "_fits")
make_figure.main(stem + "_fits.csv", stem)
