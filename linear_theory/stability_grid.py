"""Compute growth rate, frequency, optimal perturbation slope and the TKE budget over a
parameter grid, with the corrected Eq. 2.27 system (floquet.py).

Grids at REFINE = 1 (the same as the preprint figures used):
  sgamma  S in (0, 2) step 0.005, gamma in [0, 1] (5 x 40 points), delta in {0, 0.5, 1}
          -> Figs 3, 4, 5, 8, 9
  freq    S in {0.5, 1, 1.5}, omega_min^2/f*^2 from -0.05 to 0.45 (251 points),
          delta in [0, 1] (201 points)  -> Fig 6
REFINE = n divides the grid spacing by n in S and gamma (sgamma) or in omega and delta (freq),
costing about n^2 times as much.

Usage: python stability_grid.py GRID CHUNK NCHUNKS OUTDIR [NPROC] [REFINE]
  Splits the S values of GRID into NCHUNKS pieces and computes piece CHUNK (0-based),
  writing OUTDIR/GRID_CHUNK.nc. Combine the pieces with merge_grid.py.
"""
import os
import sys
from multiprocessing import Pool

import numpy as np
import xarray as xr

import floquet as fl

M_GRID = np.arange(-50, 5 + 0.1, 0.1)


def grid(name, refine=1):
    n = int(refine)
    if name == "sgamma":
        dS = 0.005 / n
        S = np.round(np.arange(dS, 2, dS), 6)
        delta = np.array([0.0, 0.5, 1.0])
        gamma = np.unique(np.concatenate([np.linspace(a, a + 0.2, 39 * n + 1) for a in np.arange(0, 1, 0.2)]))
        strat = np.broadcast_to(gamma, (S.size, gamma.size))
    elif name == "freq":
        S = np.array([0.5, 1.0, 1.5])
        delta = np.linspace(0, 1, 200 * n + 1)
        r = np.linspace(-0.05, 0.45, 250 * n + 1)          # omega_min^2 / f*^2
        strat = np.array([(1 - (1 + s**2) * r) / ((1 - r) * (1 + s**2)) for s in S])
    else:
        raise ValueError(name)
    return S, delta, strat


def point(args):
    S, gamma, delta = args
    flow = fl.Flow(S, gamma, delta)
    sig, freq, m = fl.fastest_mode(flow, M_GRID)
    out = {"growth_rate": sig, "frequency": freq, "slope_angle": m}
    for ic in ("floquet", "ic01"):
        for k, v in fl.budget(flow, m, ic).items():
            out[f"{k}_{ic}"] = v
    return out


def main(name, chunk, nchunks, outdir, nproc=None, refine=1):
    chunk, nchunks = int(chunk), int(nchunks)
    S, delta, strat = grid(name, refine)
    idx = np.array_split(np.arange(S.size), nchunks)[chunk]
    tasks = [(S[i], strat[i, j], d) for i in idx for d in delta for j in range(strat.shape[1])]
    with Pool(int(nproc) if nproc else os.cpu_count()) as pool:
        res = pool.map(point, tasks, chunksize=8)
    shape = (idx.size, delta.size, strat.shape[1])
    data = {k: (("slope_burger_number", "delta", "strat_index"), np.array([r[k] for r in res]).reshape(shape))
            for k in res[0]}
    data["strat_values"] = (("slope_burger_number", "delta", "strat_index"),
                            np.broadcast_to(strat[idx][:, None, :], shape).copy())
    ds = xr.Dataset(data, coords={"slope_burger_number": S[idx], "delta": delta,
                                  "strat_index": np.arange(strat.shape[1])})
    ds.attrs.update(equations="Eq 2.26 with corrected Eq 2.27 (linear_theory/floquet.py)",
                    units="growth_rate, frequency in units of f* = f sqrt(1+S^2)",
                    budget="Eq 2.31 terms integrated over one period / change in k_hat; "
                           "_floquet: growing Floquet mode, _ic01: (Psi, zeta) = (0, 1)",
                    refine=int(refine))
    S2 = ds.slope_burger_number**2
    with np.errstate(divide="ignore", invalid="ignore"):
        ds["omega_ratio"] = (1 - ds.strat_values * (1 + S2)) / ((1 + S2) * (1 - ds.strat_values))
    ds["omega_ratio"].attrs["long_name"] = "omega_min^2 / f*^2"
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{name}_{chunk:02d}.nc")
    ds.to_netcdf(path)
    print(f"wrote {path}: {len(tasks)} points")


if __name__ == "__main__":
    main(*sys.argv[1:])
