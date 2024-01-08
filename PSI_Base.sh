#!/bin/bash -l
### Job Name
#PBS -N PSI_Base
### Project Code Allocation
#PBS -A UMCP0023
### Resources
#PBS -l select=1:ncpus=1:ngpus=1
### Run Time
#PBS -l walltime=20:00:00
### To the casper queue
#PBS -q casper
### output
#PBS -o PSI_b_change.out
### error
#PBS -j oe
### type of GPU
#PBS -l gpu_type=v100
### email 
#PBS -M knudsen@umd.edu
#PBS -m abe

## Clear environent and Load all the modules needed
module purge
module load ncarenv/1.3 gnu/10.1.0 ncarcompilers/0.5.0
module load netcdf/4.8.1 openmpi/4.1.1 
module load julia
module load cuda/4.4.1

### file to run
julia --project PSI_Base.jl /glade/scratch/knudsenl/BottomBoundaryLayer/


