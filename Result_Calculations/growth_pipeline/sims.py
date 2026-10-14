"""Simulation parameters from Table 1 of the manuscript (and parameters_stblty.jl).

SI-unstable runs are those with Ri_HT^{-1} > 1, i.e. omega_min^2/(f*)^2 < 0 (starred in Table 1).
"""

def gamma(S, Ri_inv):
    return (1 + (1 - 1 / Ri_inv) * S**2) / (1 + S**2)


def _sim(S, Ri_inv, delta):
    return {"S": S, "Ri_inv": Ri_inv, "delta": delta, "si_unstable": Ri_inv > 1}


SIMS = {
    "Sim1": _sim(0.5, 0.65, 0.8), "Sim2": _sim(0.5, 0.75, 0.8), "Sim3": _sim(0.5, 0.85, 0.8),
    "Sim4": _sim(0.5, 0.95, 0.8), "Sim5": _sim(0.5, 1.05, 0.8),
    "Sim6": _sim(1.0, 0.65, 0.8), "Sim7": _sim(1.0, 0.75, 0.8), "Sim8": _sim(1.0, 0.85, 0.8),
    "Sim9": _sim(1.0, 0.95, 0.8), "Sim10": _sim(1.0, 1.05, 0.8),
    "Sim11": _sim(1.5, 0.75, 0.8), "Sim12": _sim(1.5, 0.85, 0.8), "Sim13": _sim(1.5, 0.95, 0.8),
    "Sim14": _sim(1.5, 1.05, 0.8),
    "Sim15": _sim(1.0, 0.8, 0.8), "Sim16": _sim(1.0, 0.8, 0.6), "Sim17": _sim(1.0, 0.8, 0.4),
    "Sim18": _sim(1.0, 0.8, 0.2),
}

# File naming used by full_code_diagnostics_stblty_comp.jl for the March runs, e.g.
# TKE_terms_Sinf_1.0_Ri_inv_0.95_delta_0.8.nc (Sim15-18 were written with Sinf_1).
def tke_filename(s):
    S = "1" if (s["Ri_inv"] == 0.8) else str(s["S"])
    return f"TKE_terms_Sinf_{S}_Ri_inv_{s['Ri_inv']}_delta_{s['delta']}.nc"


# Runs whose fastest Floquet mode is not PSI (Im sigma != f*/2), from theory_growth.py:
# the three SI-unstable runs (Im sigma = 0) and Sim6 (omega_min^2/(f*)^2 = 0.35, Im sigma = 0).
NOT_PSI = {"Sim5", "Sim6", "Sim10", "Sim14"}
