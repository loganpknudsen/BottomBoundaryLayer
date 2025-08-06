from os import system
import math
import numpy as np
#+++ Define simnames

############################## f1e4theta029N21e5gammau
#   FREQF-THETA-N2_gamma.jld2
all_sims = ['S01-gammau'
            #'f1e4-theta029-N21e5-delta05-Vinf005-gammau',
            #  'f1e4-theta029-N21e5-delta05-Vinf005-gammam',
            #  'f1e4-theta029-N21e5-delta05-Vinf005-gammal',
            #  'f1e4-theta10-N21e5-delta05-Vinf02-gammau',
            #  'f1e4-theta60-N21e6-delta05-Vinf005-gammau',
            #  'f1e4-theta05-N21e7-delta05-Vinf0001-gammau',
            #  'f1e4-theta0009-N21e5-delta05-Vinf0002-gammau',
            #  'f1e4-theta250-N21e5-delta05-Vinf02-gammau',
            #  'f1e4-theta250-N21e5-delta05-Vinf02-gammam',
            #  'f1e4-theta250-N21e5-delta05-Vinf02-gammal',
            #  'f1e4-theta20-N21e5-delta05-Vinf02-gammau',
            #  'f1e4-theta20-N21e5-delta05-Vinf02-gammam',
            #  'f1e4-theta20-N21e5-delta05-Vinf02-gammal',
            #  'f1e4-theta20-N21e5-delta025-Vinf02-gammau',
            #  'f1e4-theta20-N21e5-delta075-Vinf02-gammau',
            #  'f1e4-theta30-N21e6-delta05-Vinf005-gammau',
            #  'f1e4-theta30-N21e6-delta05-Vinf005-gammam',
            #  'f1e4-theta30-N21e6-delta05-Vinf005-gammal',
            #  'f1e4-theta30-N21e6-delta025-Vinf005-gammau',
            #  'f1e4-theta30-N21e6-delta075-Vinf005-gammau',
            #  'f1e4-theta01812-N21e5-delta05-Vinf01-gammau',
            #  'f1e4-theta01812-N21e5-delta05-Vinf01-gammam',
            #  'f1e4-theta01812-N21e5-delta05-Vinf01-gammal',
            #  'f1e4-theta09059-N21e5-delta05-Vinf01-gammau',
            #  'f1e4-theta09059-N21e5-delta05-Vinf01-gammam',
            #  'f1e4-theta09059-N21e5-delta05-Vinf01-gammal',
            #  'f1e4-theta09059-N21e5-delta025-Vinf01-gammau',
            #  'f1e4-theta09059-N21e5-delta075-Vinf01-gammau',
            #  'f1e4-theta181130-N21e5-delta05-Vinf01-gammau',
            #  'f1e4-theta181130-N21e5-delta05-Vinf01-gammam',
            #  'f1e4-theta181130-N21e5-delta05-Vinf01-gammal',
            #  'f1e4-theta181130-N21e5-delta025-Vinf01-gammau',
            #  'f1e4-theta181130-N21e5-delta075-Vinf01-gammau',
            #  'f1e4-theta27160-N21e5-delta05-Vinf01-gammau',
            #  'f1e4-theta27160-N21e5-delta05-Vinf01-gammam',
            #  'f1e4-theta27160-N21e5-delta05-Vinf01-gammal',
            #  'f1e4-theta27160-N21e5-delta025-Vinf01-gammau',
            #  'f1e4-theta27160-N21e5-delta075-Vinf01-gammau',
            #  'f1e4-theta36190-N21e5-delta05-Vinf01-gammau',
            #  'f1e4-theta36190-N21e5-delta05-Vinf01-gammam',
            #  'f1e4-theta36190-N21e5-delta05-Vinf01-gamma005',
            #  'f1e4-theta04530-N21e5-delta05-Vinf01-gammau',
            #  'f1e4-theta04530-N21e5-delta05-Vinf01-gammam',
            #  'f1e4-theta04530-N21e5-delta05-Vinf01-gammal',
            #  'f1e4-theta13587-N21e5-delta05-Vinf01-gammau',
            #  'f1e4-theta13587-N21e5-delta05-Vinf01-gammam',
            #  'f1e4-theta13587-N21e5-delta05-Vinf01-gammal',
            #  'f1e4-theta22637-N21e5-delta05-Vinf01-gammau',
            #  'f1e4-theta22637-N21e5-delta05-Vinf01-gammam',
            #  'f1e4-theta22637-N21e5-delta05-Vinf01-gammal',
            #  'f1e4-theta31676-N21e5-delta05-Vinf01-gammau',
            #  'f1e4-theta31676-N21e5-delta05-Vinf01-gammam',
            #  'f1e4-theta31676-N21e5-delta05-Vinf01-gamma005',
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
    print(Sinf)
    if params[1] == "gammau":
        gm = (1+Sinf**2)**(-1)
    elif params[1] == "gammal":
        gm = (3-Sinf**2)*(3*(1+Sinf**2))**(-1)
    elif params[1] == "gammam":
        gm = ((1+Sinf**2)**(-1)+(3-Sinf**2)*(3*(1+Sinf**2))**(-1))/2
    elif params[1] == "gamma005":
        gm = 0.05
    return Sinf, gm
    
for sim in all_sims:

    #+++ Define simulation name
    #simname_full = f"IntWave-{dims}-{resScale}-" + sim["wavelength"] + "-" + sim["flowspd"]
    simname_full = sim
    
    Sinf, gm = parseNaming(simname_full)
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
