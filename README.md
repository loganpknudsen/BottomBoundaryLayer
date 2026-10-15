# Parametric subharmonic instability in the ocean bottom boundary layer

Code for Knudsen, Wenegrat, Hilditch & Thomas, *Parametric Subharmonic Instability in the Ocean
Bottom Boundary Layer*: the nonlinear Oceananigans simulations, the Floquet linear stability
analysis, the growth-rate comparison, and one script per manuscript figure.

The workflow below is the one used for the revision (branch `refit-growth-rates`). The original
notebooks it replaced are listed at the end.

## Layout

| Path | Contents |
|---|---|
| `full_code_diagnostics_noflux.jl` | Nonlinear simulation, no perturbation buoyancy flux at the bottom (∂b′/∂z = 0) |
| `full_code_diagnostics_stblty_comp_drag.jl` | The same with quadratic bottom drag (Sim9D, Fig 17) |
| `parameters_stblty.jl` | Parameter sets for every run (Table 1) |
| `submit_noflux.py` | Writes and submits the PBS jobs for the runs above |
| `linear_theory/` | Floquet analysis (Eq 2.26 with the Eq 2.28 coefficients) and the parameter-space grids |
| `Result_Calculations/growth_pipeline/` | Growth-rate fits of the simulations and their comparison with theory (Fig 10) |
| `figures/` | One script per figure, `make_all.py`, and shared helpers |

## 1. Simulations

Each run is a parameter set in `parameters_stblty.jl` (e.g. `S10R095delta08` is Sim9:
S∞ = 1, Ri_inv = 0.95, δ = 0.8). `Ri_inv` is the Hilditch & Thomas inverse Richardson number
Ri_HT⁻¹ of Eq 2.29 (so ω²min/f*² = 1 − Ri_inv), not the inverse of the boundary-layer Richardson
number N²(1 − γ)/Λ²; see the note at the top of `parameters_stblty.jl`. Submit from the repo root on Casper:

```
python submit_noflux.py                         # all runs in the growth-rate comparison
python submit_noflux.py S10R095delta08          # only these
python submit_noflux.py S05R105delta08          # SI-unstable runs are submitted only when named
python submit_noflux.py --drag S10R095delta08   # Sim9D, with bottom drag
python submit_noflux.py --dry-run ...           # write and print the job scripts only
```

Output goes to `OUTPUT_DIR`, set at the top of `submit_noflux.py`
(currently `/glade/derecho/scratch/wenegrat/BottomBoundaryLayer/output`; change it for another user):
`TKE_terms_*.nc` (domain-averaged TKE budget, every 0.005 inertial periods) and
`flow_fields_*.nc` (fields, every 0.05 periods); drag runs add `_drag`. Progress is written live
to `logs/noflux/<run>.log`.

The bottom boundary condition is ∂b′/∂z = 0, as stated in the manuscript. The March 2026 runs
used b′ = 0 (`full_code_diagnostics_stblty_comp.jl`, submitted with `run_all.sh`); growth rates
differ by at most 1% between the two.

## 2. Linear stability analysis

`linear_theory/floquet.py` solves the Floquet problem and computes the TKE budget (Eq 2.31) on the
same Eq 2.28 system as the growth rates; v is reconstructed consistently with the small-angle
approximation, so the budget closes exactly. `stability_grid.py` evaluates growth rate,
frequency, optimal slope and budget over a parameter grid:

- `sgamma`: S∞–γ, δ ∈ {0, 0.5, 1} (Figs 3, 4, 5, 8, 9)
- `freq`: ω²min/f*²–δ for S∞ ∈ {0.5, 1, 1.5} (Fig 6)

Submit both, chunked, and merge them:

```
DM=0.01 HOURS=8 bash linear_theory/submit_grids.sh sgamma freq
```

`DM` is the step of the scan over perturbation slope; use 0.01 (the preprint used 0.1, which
misses narrow instability ranges at small δ). `REFINE=n` refines the grid n-fold. The figures read
the merged grids from `/glade/derecho/scratch/$USER/BottomBoundaryLayer/output/linear_theory_dm0.01`.

## 3. Growth rates (Fig 10)

`Result_Calculations/growth_pipeline/`:

- `theory_growth.py`: Floquet growth rate of the fastest mode at each run's exact (S, γ, δ)
- `diagnose_growth.py`: the one-period growth rate σ₁P(T) = [ln k(T+1) − ln k(T)]/(4π), in units of f*
- `fit_growth.py`: fits ln k over a window of 4 amplitude e-folds at the peak of σ₁P,
  L = 4/(2π σ_peak) whole periods, ending at the peak (n = 3 and 5 are also fitted, for sensitivity).
  Runs whose window would start before T = 0 have no exponential phase and are flagged.
- `make_figure.py`: plots simulation against theory; SI-unstable runs are included, in grey (`--no-si` drops them)
- `sims.py`: run parameters and file names

`figures/10_growth_rate_comparison.py` runs the fit and the plot together, so the figure is always
refitted from the current output.

## 4. Figures

Every manuscript figure has its own script in `figures/`, named by figure number:

| Script | Figure | Data |
|---|---|---|
| `02_hodographs.py` | Near-inertial oscillation hodographs | analytic |
| `03_parameter_space.py` | PSI region in S∞–γ | `sgamma` grid |
| `04_growth_rate_s_gamma.py` | Growth rate over S∞–γ | `sgamma` grid |
| `05_slope_angle_s_gamma.py` | Optimal perturbation slope | `sgamma` grid |
| `06_growth_rate_omega_delta.py` | Growth rate over ω²min/f*²–δ | `freq` grid |
| `07_tke_budget_linear.py` | TKE budget of the linear solution | computed (`floquet.py`) |
| `08_relative_contributions.py` | ASP, GSP, BP fractions of the TKE change | `sgamma` grid |
| `09_gsp_plus_bp.py` | GSP + BP fraction | `sgamma` grid |
| `10_growth_rate_comparison.py` | Simulated vs. Floquet growth rates | all `TKE_terms_*.nc` |
| `11_u_linear_phase_sim9.py` | u in the linear phase, Sim9 | Sim9 `flow_fields` |
| `12_spectrum_sim9.py` | Quasi-Lagrangian spectrum and time series, Sim9 | Sim9 `flow_fields` |
| `13_u_finite_amplitude_sim9.py` | u at finite amplitude and Ri ≤ 0.25 fraction, Sim9 | Sim9 `flow_fields` |
| `14_tke_budget_sim9.py` | Integrated TKE budget, Sim9 | Sim9 `TKE_terms` |
| `15_u_linear_phase_sim7.py` | u in the linear phase, Sim7 | Sim7 `flow_fields` |
| `16_tke_budget_sim7.py` | Integrated TKE budget, Sim7 | Sim7 `TKE_terms` |
| `17_drag_sim9d.py` | Sim9D: u snapshots and quasi-Lagrangian spectrum | Sim9D `flow_fields` |

Shared modules:

- `paths.py`: data and output locations (see below) and `save()`
- `stability_common.py`: loading and masking of the linear-theory grids (Figs 3–6, 8, 9)
- `quasi_lagrangian.py`: quasi-Lagrangian sampling (Eq 3.3) and spectrum (Figs 12, 17)
- `tke_budget_sim.py`: the integrated-budget plot (Figs 14, 16)
- `colour_scales.py`: colour limits from the data (Figs 11, 15)

### Generating the figures

On a Casper login node:

```
module load conda && conda activate npl
python figures/make_all.py              # every figure
python figures/make_all.py 10 14 17     # only these
python figures/make_all.py --residual   # also the budget-residual versions of Figs 7, 14, 16
python figures/12_spectrum_sim9.py      # or any single script
```

`make_all.py` runs each script in its own process, so one failure (for example a run that has not
finished) does not stop the others, and prints a summary of what failed. It turns off HDF5 file
locking, so it is safe to run while simulations are still writing their output.

Figures are written to `figures/output/` as `NN_name.png` and `.pdf`. Fig 10 also writes its fits
(`10_growth_rate_comparison_fits.csv`) and a plot of the fit windows.

Figs 12 and 17 need `xrft` in the Python environment.

### Data locations

Set in `figures/paths.py`; override any of them with an environment variable:

| Variable | Default | Used by |
|---|---|---|
| `BBL_SIM_DIR` | `/glade/derecho/scratch/$USER/BottomBoundaryLayer/output` | Figs 10–17 |
| `BBL_LT_DIR` | `/glade/derecho/scratch/$USER/BottomBoundaryLayer/output/linear_theory_dm0.01` | Figs 3–6, 8, 9 |
| `BBL_FIG_OUT` | `figures/output` | all |

For example, to plot the March (b′ = 0) runs instead:
`BBL_SIM_DIR=/glade/derecho/scratch/knudsenl/data/new_data/final_draft python figures/make_all.py 11`.

## Environment (NCAR Casper)

- Python: `module load conda && conda activate npl`
- Julia: `module load ncarenv/23.10 gcc ncarcompilers netcdf cuda julia/1.10.2`; packages from
  `Project.toml`/`Manifest.toml` (Oceananigans, Oceanostics)

## Earlier code

These produced the submitted version and are kept for reference; the scripts above replace them.

- `run_all.sh` / `run_all.py` / `aux_pbs_psi.sh`: the March runs (b′ = 0 at the bottom) with
  `full_code_diagnostics_stblty_comp.jl`. `full_code_diagnostics.jl` with `sub.sh` runs a single simulation.
- `stability_analysis_full_run*.py` with `stability.sh`: the original parameter-space stability
  analysis; `stability_analysis_output_reader*.ipynb` read its output and drew Figs 3–9.
- `Floquet_Analysis_*.ipynb`, `oscillations_mangitudes.ipynb`: linear-theory notebooks (Figs 2, 7).
- `stability_analysis_against_simulation.ipynb`: the original growth-rate comparison (Fig 10).
- `Result_Calculations/*.ipynb`: notebooks for the simulation figures
  (`velocity_visualization`, `time_series_and_fourier_transform`, `TKE_terms`, `drag_figure`,
  `BFLUX_GSP_comp`, `TKE_terms_velocity_comparison`).
- `Result_Calculations/refit_growth.py`: one-off check of the notebook's fixed fit windows, superseded by `growth_pipeline/`.
- `Old_Code/`: superseded versions.
