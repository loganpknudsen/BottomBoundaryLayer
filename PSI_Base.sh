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
module --force purge
module load ncarenv/23.10 gcc/12.2.0
module load ncarcompilers/1.0.0
module load netcdf/4.9.2 openmpi/4.1.6 
module load julia
module load cuda/11.8.0

### file to run
julia --project PSI_Base.jl /glade/derecho/scratch/knudsenl/BottomBoundaryLayer


