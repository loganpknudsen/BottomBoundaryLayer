"""Figure 14: time-integrated TKE budget terms for Sim9, from T = 5 to 20 (see tke_budget_sim.py).

Run: Sim9 (S = 1, Ri_inv = 0.95, delta = 0.8), TKE_terms_..._Ri_inv_0.95_delta_0.8.nc in BBL_SIM_DIR.
With --residual, adds the budget residual and saves as 14_tke_budget_sim9_residual.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import tke_budget_sim  # noqa: E402

tke_budget_sim.plot("TKE_terms_Sinf_1.0_Ri_inv_0.95_delta_0.8.nc", S=1.0, tf=3000, script_file=__file__)
