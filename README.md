Hello, if you are reading this you are viewing the repository which contains the code used to study parametric subharmonic instability in the oceanic bottom boundary layer.
The important files which made this possible and general organization are described below:\
1. To run the full nonlinear simulations, there are two options:
     - To run all the simulations in the parameter space, run the shell script  `run_all.sh` which will automatically run the simulations using `run_all.py`, `full_code_diagnostics_stblty_comp.jl`, and `aux_pbs_psi.sh` for the parameters defined in `parameters.jl`.
     - To run a single simulation, modify `full_code_diagnostics.jl` and run the shell script `sub.sh`.
2. The linear stability analysis over the parameter space is contained in the code `stability_analysis_full_run.py` where submitting the shell script `stability.sh` will execute the code.
3. 
