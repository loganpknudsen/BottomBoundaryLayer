from os import system
import math
import numpy as np
#+++ Define simnames

##############################
#   FREQF-THETA-N2_gamma.jld2
all_sims = [ '1f-0.29-1e5-gammal',
             'f-0.29-1e5-gammau']


#+++ Options
remove_checkpoints = False
only_one_job = False
dry_run = False

IPeriods = "30.0"
f = 1e-4 
    
verbose = 2
aux_filename = "aux_pbs_psi.sh"
julia_file = "full_code_diagnostics.jl"
savepath = "/glade/derecho/scratch/knudsenl/data/new_data/paper_data/"


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

###module purge                                                                                                                                  
###module load ncarenv/1.3 gnu/9.1.0 ncarcompilers/0.5.0                                                                                         
###module load netcdf/4.7.4 openmpi/4.1.0                                                                                                        
###module load cuda/11.4.0                                                                                                                       

module --force purge
module --ignore-cache load ncarenv/23.10 gcc ncarcompilers netcdf
module --ignore-cache load cuda
module --ignore-cache load julia/1.10.2

### file to run                    

julia --pkgimages=no --project=../Project.toml {julia_file} --path {savepath} --strat {strat} --theta {theta} --freqf {freqf} --suffix\
    {simname_full} -T {IPeriods} | tee logs/{simname_full}.out
"""

def parseNaming(name):
    
    params = name.split('_')[-1].split('-')
    
    freqf = float(params[0].replace('f',''))
    theta = float(params[1])
    strat = 1 * 10**(-1 * float(params[3].replace('1e','')))
    if params[3] == "gammau":
        gamma = (1+0.5*(strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8))**(-1)
    elif params[3] == "gammal":
        gamma = (3-strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8)*(3*(1+strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8)-2*strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8)**(-1)
    return freqf, theta, strat, gamma
    
for sim in all_sims:

    #+++ Define simulation name
    #simname_full = f"IntWave-{dims}-{resScale}-" + sim["wavelength"] + "-" + sim["flowspd"]
    simname_full = sim 
    
    freqf, theta, strat, gamma= parseNaming(simname_full)
    # print(f)
    # freq = f'{freqf*f:10}'
    pbs_script_filled = pbs_script.format(simname_full=simname_full, savepath=savepath, julia_file=julia_file, IPeriods=IPeriods,
                                          freqf=freqf, theta=theta, strat=strat, gamma=gamma)

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
