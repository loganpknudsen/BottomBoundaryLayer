"""Figure 10: simulated vs. Floquet growth rates.

Runs Result_Calculations/growth_pipeline/make_figure.py on the fits from fit_growth.py
(BBL_FITS_CSV, default Result_Calculations/growth_pipeline/growth_fits.csv). Fit window: n = 4
amplitude e-folds ending at the peak growth rate. SI-unstable runs are included, in grey; runs
with no exponential phase are plotted at zero with an open marker.
"""
import sys

from paths import FITS_CSV, OUT_DIR, REPO

sys.path.insert(0, str(REPO / "Result_Calculations" / "growth_pipeline"))
import make_figure  # noqa: E402

OUT_DIR.mkdir(parents=True, exist_ok=True)
make_figure.main(str(FITS_CSV), str(OUT_DIR / "10_growth_rate_comparison"))
