#!/bin/bash -l
### Job Name
#PBS -N tiltedBBL
### Project Code Allocation
#PBS -A UMCP0023
### Resources
#PBS -l select=1:ncpus=16:mem=80GB:ngpus=1
### Run Time
#PBS -l walltime=10:00:00
### To the casper queue
#PBS -q casper
### output
#PBS -o Sim3.out
### error
#PBS -e Sim3.err
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

julia --pkgimages=no --project=. full_code_diagnostics_bac