#!/bin/bash -l
### Job Name
#PBS -N stability2
### Project Code Allocation
#PBS -A UMCP0023
### Resources
#PBS -l select=1:ncpus=1:mem=80GB
### Run Time
#PBS -l walltime=24:00:00
### To the casper queue
#PBS -q casper
### output
#PBS -o stability2.out
### error
#PBS -e stability2.err
### email 
#PBS -M knudsen@umd.edu
#PBS -m abe


module --force purge
module --ignore-cache load ncarenv-basic/23.10
module --ignore-cache load conda
conda activate dedalus3
### file to run

python3 -u stability_analysis_full_run_small_angle.py
