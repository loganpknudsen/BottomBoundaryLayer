#!/bin/bash -l
### Job Name
#PBS -N PSI_Base
### Project Code Allocation
#PBS -A UMCP0023
### Resources
#PBS -l select=1:ncpus=1
### Run Time
#PBS -l walltime=10:00:00
### To the casper queue
#PBS -q casper
### output
#PBS -o PSI_Base.out
### error
#PBS -j oe
### email 
#PBS -M knudsen@umd.edu
#PBS -m abe

## Clear environent and Load all the modules needed
module purge
module load ncarenv/1.3 gnu/10.1.0 ncarcompilers/0.5.0
module load netcdf/4.8.1 openmpi/4.1.1 
module load julia

### file to run
julia --project PSI_Base.jl /glade/scratch/knudsenl/BottomBoundaryLayer/



