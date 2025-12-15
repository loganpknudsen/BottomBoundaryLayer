                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
#!/bin/bash -l
### Job Name
#PBS -N S10R080delta02  
### Project Code Allocation
#PBS -A UMCP0023
### Resources :ngpus=1
#PBS -l select=1:ncpus=128:mpiprocs=4:ngpus=4:gpu_type=a100
### Run Time
#PBS -l walltime=8:00:00
### To the casper queue
#PBS -q casper
### output
#PBS -o logs/S10R080delta02.out 
### error
#PBS -e logs/S10R080delta02.err
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

julia --pkgimages=no --project=. full_code_diagnostics_stblty_comp.jl --path /glade/derecho/scratch/knudsenl/data/new_data/ --Sinf 1.0  --suffix    S10R080delta02 -T 30.0
