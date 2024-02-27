#!/bin/bash -l
### Job Name
#PBS -N BBL_w_o
### Project Code Allocation
#PBS -A UMCP0023
### Resources
#PBS -l select=1:ncpus=1:ngpus=1
### Run Time
#PBS -l walltime=20:00:00
### To the casper queue
#PBS -q casper
### output
#PBS -o BBL_O.out
### error
#PBS -j oe
### type of GPU
#PBS -l gpu_type=v100
### email 
#PBS -M knudsen@umd.edu
#PBS -m abe

## Clear environent and Load all the modules needed
module purge
module load ncarenv/23.10 ncarcompilers/1.0.0
module load netcdf/4.9.2 openmpi/4.1.6 
module load julia
module load cuda/12.2.1

### file to run
julia --project BBL_with_oscillations_code.jl /glade/derecho/scratch/knudsenl/BottomBoundaryLayer/


