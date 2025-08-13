from os import system
import math
import numpy as np
#+++ Define simnames

############################## f1e4theta029N21e5gammau
#   FREQF-THETA-N2_gamma.jld2
all_sims = [ #'S05-R065-delta08',
            # 'S05-R075-delta08',
            # 'S05-R085-delta08',
            # 'S05-R095-delta08',
            # 'S05-R105-delta08',
            # 'S10-R065-delta08',
            # 'S10-R075-delta08',
            # 'S10-R085-delta08',
            # 'S10-R095-delta08',
            # 'S10-R105-delta08',
            # 'S15-R075-delta08',
            # 'S15-R085-delta08',
            # 'S15-R095-delta08',
            # 'S15-R105-delta08',
            # 'S05-R080-delta08',
            # 'S05-R080-delta06',
            # 'S05-R080-delta04',
            # 'S05-R080-delta02',
            # 'S10-R080-delta08',
            # 'S10-R080-delta06',
            'S10-R080-delta04',
            # 'S10-R080-delta02',
            # 'S15-R080-delta08',
            # 'S15-R080-delta06',
            # 'S15-R080-delta04',
            # 'S15-R080-delta02',
            # 'S01-gammau',
            # 'S01-gammam',
            # 'S025-gammal',
            # 'S025-gammau',
            # 'S025-gammam',
            # 'S025-gammal',
            # 'S05-gammau',
            # 'S05-gammam',
            # 'S05-gammal',
            # 'S075-gammau',
            # 'S075-gammam',
            # 'S075-gammal',
            # 'S10-gammau',
            # 'S10-gammam',
            # 'S10-gammal',
            # 'S125-gammau',
            # 'S125-gammam',
            # 'S125-gammal',
            # 'S15-gammau',
            # 'S15-gammam',
            # 'S15-gammal',
            # 'S175-gammau',
            # 'S175-gammam',
            # 'S175-gamma005',
            # 'S20-gammau',
            # 'S20-gammam',
            # 'S20-gamma005',
             ]


#+++ Options
remove_checkpoints = False
only_one_job = False
dry_run = False

IPeriods = "30.0"
f = 1e-4 
    
verbose = 2
aux_filename = "aux_pbs_psi.sh"
julia_file = "full_code_diagnostics_stblty_comp.jl"
savepath = "/glade/derecho/scratch/knudsenl/data/new_data/"


#---

pbs_script = \
"""                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
#!/bin/bash -l
### Job Name
#PBS -N {simname_full}  
### Project Code Allocation
#PBS -A UMCP0023
### Resources :ngpus=1
#PBS -l select=1:ncpus=16:mem=80GB:ngpus=1
### Run Time
#PBS -l walltime=8:00:00
### To the casper queue
#PBS -q casper
### output
#PBS -o logs/{simname_full}.out 
### error
#PBS -e logs/{simname_full}.err
### type of GPU
#PBS -l gpu_type=v100
### email 
#PBS -M knudsen@umd.edu
#PBS -m abe                                                                                                                                                                                                                                                        

module --force purge
module --ignore-cache load ncarenv/23.10 gcc ncarcompilers netcdf
module --ignore-cache load cuda
module --ignore-cache load julia/1.10.2

### file to run                    

julia --pkgimages=no --project=. {julia_file} --path {savepath} --Sinf {Sinf}  --suffix\
    {simname_full} -T {IPeriods}
"""

def parseNaming(name):
    
    params = name.split('_')[-1].split('-')
    
    Sinf = float(params[0].replace('S',''))*10**(-1)
    Ri_inv = float(params[1].replace('R',''))*10**(-2)
    dlt = float(params[2].replace('delta',''))*10**(-1)
    # if params[1] == "gammau":
    #     gm = (1+Sinf**2)**(-1)
    # elif params[1] == "gammal":
    #     gm = (3-Sinf**2)*(3*(1+Sinf**2))**(-1)
    # elif params[1] == "gammam":
    #     gm = ((1+Sinf**2)**(-1)+(3-Sinf**2)*(3*(1+Sinf**2))**(-1))/2
    # elif params[1] == "gamma005":
    #     gm = 0.05
    return Sinf, Ri_inv, dlt
    
for sim in all_sims:

    #+++ Define simulation name
    #simname_full = f"IntWave-{dims}-{resScale}-" + sim["wavelength"] + "-" + sim["flowspd"]
    simname_full = sim
    
    Sinf, Ri_inv, delta = parseNaming(simname_full)
    simname_full = simname_full.replace("-","") 
    # print(f)
    # freq = f'{freqf*f:10}'
    pbs_script_filled = pbs_script.format(simname_full=simname_full, savepath=savepath, julia_file=julia_file, IPeriods=IPeriods,
                                          Sinf=Sinf)

    cmd1 = f"qsub {aux_filename}"
    if verbose>1: print(pbs_script_filled)
    if verbose>0: print(cmd1)
    #---

    #+++ Run command
    if not dry_run:
        with open(aux_filename, "w") as f:
            f.write(pbs_script_filled)
        system(cmd1)
    #---

    print()
