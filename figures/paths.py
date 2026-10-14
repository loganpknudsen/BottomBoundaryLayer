"""Where the figure scripts read data and write figures.

Override any location with an environment variable, e.g.
    BBL_SIM_DIR=/path/to/runs python figures/11_u_linear_phase_sim9.py
"""
import os
from pathlib import Path

USER = os.environ.get("USER", "")
FIG_DIR = Path(__file__).resolve().parent
REPO = FIG_DIR.parent

# Simulation output (TKE_terms_*.nc and flow_fields_*.nc). Default: the March 2026 runs.
SIM_DIR = Path(os.environ.get("BBL_SIM_DIR", "/glade/derecho/scratch/knudsenl/data/new_data/final_draft"))

# Linear-theory grids from linear_theory/stability_grid.py (sgamma.nc, freq.nc).
LT_DIR = Path(os.environ.get("BBL_LT_DIR", f"/glade/derecho/scratch/{USER}/BottomBoundaryLayer/output/linear_theory"))

# Growth-rate fits from Result_Calculations/growth_pipeline/fit_growth.py.
FITS_CSV = Path(os.environ.get("BBL_FITS_CSV", REPO / "Result_Calculations/growth_pipeline/growth_fits.csv"))

OUT_DIR = Path(os.environ.get("BBL_FIG_OUT", FIG_DIR / "output"))


def save(fig, script_file, dpi=300):
    """Save `fig` as figures/output/<script name>.png and .pdf."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = Path(script_file).stem
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"{stem}.{ext}", dpi=dpi)
    print(f"wrote {OUT_DIR / stem}.png/.pdf")
