#!/bin/bash -l
### Job Name
#PBS -N testcuda
### Project Code Allocation
#PBS -A UMCP0023
### Resources
#PBS -l select=1:ncpus=16:mem=80GB:ngpus=1
### Run Time
#PBS -l walltime=00:10:00
### To the casper queue
#PBS -q gpudev
### output
#PBS -o Sim.out
### error
#PBS -e Sim.err
### type of GPU
#PBS -l gpu_type=v100
### email 
#PBS -M knudsen@umd.edu
#PBS -m abe


module --force purge
module --ignore-cache load ncarenv/23.10 gcc ncarcompilers netcdf
module --ignore-cache load cuda
module --ignore-cache load julia/1.9 
### file to run

julia --pkgimages=no --project=. tilted_bottom_boundary_layer.jl