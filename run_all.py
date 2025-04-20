from os import system
import math
import numpy as np
#+++ Define simnames

############################## f1e4theta029N21e5gammau
#   FREQF-THETA-N2_gamma.jld2
all_sims = [ #'f1e4-theta029-N21e5-delta05-Vinf005-gammau',
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
             'f1e4-theta27160-N21e5-delta05-Vinf01-gammal',
            #  'f1e4-theta27160-N21e5-delta025-Vinf01-gammau',
            #  'f1e4-theta27160-N21e5-delta075-Vinf01-gammau',
            #  'f1e4-theta36190-N21e5-delta05-Vinf01-gammau',
             'f1e4-theta36190-N21e5-delta05-Vinf01-gammam',
             'f1e4-theta36190-N21e5-delta05-Vinf01-gamma005',
             ]


#+++ Options
remove_checkpoints = False
only_one_job = False
dry_run = False

IPeriods = "30.0"
f = 1e-4 
    
verbose = 2
aux_filename = "aux_pbs_psi.sh"
julia_file = "full_code_diagnostics.jl"
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

julia --pkgimages=no --project=. {julia_file} --path {savepath} --strat {strat} --theta {theta} --freqf {freqf} --delta {delta} --suffix\
    {simname_full} -T {IPeriods}
"""

def parseNaming(name):
    
    params = name.split('_')[-1].split('-')
    
    freqf = float(params[0].replace('f',''))**(-1)
    theta = float(10**(-1*len(params[1].replace("theta",""))+1))*float(params[1].replace("theta",""))
    strat = 1 * 10**(-1 * float(params[2].replace('N21e','')))
    delta = float(params[3].replace('delta',''))*10**(-1)
    if params[5] == "gammau":
        gamma = (1+(1-delta)*(strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8))**(-1)
    elif params[5] == "gammal":
        gamma = (3-strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8)*(3*(1+strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8)-4*delta*strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8)**(-1)
    elif params[5] == "gammam":
        gamma = ((1+(1-delta)*(strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8))**(-1)+(3-strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8)*(3*(1+strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8)-4*delta*strat*np.tan(theta*np.pi/180)**2*freqf**(-2)*1e8)**(-1))/2
    elif params[5] == "gamma005":
        gamma = 0.05
    return freqf, theta, strat, delta, gamma
    
for sim in all_sims:

    #+++ Define simulation name
    #simname_full = f"IntWave-{dims}-{resScale}-" + sim["wavelength"] + "-" + sim["flowspd"]
    simname_full = sim
    
    freqf, theta, strat, delta, gamma= parseNaming(simname_full)
    simname_full = simname_full.replace("-","") 
    # print(f)
    # freq = f'{freqf*f:10}'
    pbs_script_filled = pbs_script.format(simname_full=simname_full, savepath=savepath, julia_file=julia_file, IPeriods=IPeriods,
                                          freqf=freqf, theta=theta, strat=strat, gamma=gamma, delta=delta)

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
