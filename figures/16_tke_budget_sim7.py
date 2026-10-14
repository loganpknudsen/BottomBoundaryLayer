"""Figure 16: time-integrated TKE budget terms for Sim7, from T = 5 to 27.5 (see tke_budget_sim.py).

Run: Sim7 (S = 1, Ri_inv = 0.75, delta = 0.8), TKE_terms_..._Ri_inv_0.75_delta_0.8.nc in BBL_SIM_DIR.
With --residual, adds the budget residual and saves as 16_tke_budget_sim7_residual.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import tke_budget_sim  # noqa: E402

tke_budget_sim.plot("TKE_terms_Sinf_1.0_Ri_inv_0.75_delta_0.8.nc", S=1.0, tf=4500, script_file=__file__)
