#!/bin/bash -l
### Job Name
#PBS -N tiltedBBL
### Project Code Allocation
#PBS -A UMCP0023
### Resources
#PBS -l select=1:ncpus=4:mem=80GB
### Run Time
#PBS -l walltime=24:00:00
### To the casper queue
#PBS -q casper
### output
#PBS -o Sim6.out
### error
#PBS -e Sim6.err
### email 
#PBS -M knudsen@umd.edu
#PBS -m abe


module --force purge
module --ignore-cache load ncarenv-basic/23.10
module --ignore-cache load conda
conda activate dedalus3
### file to run

python3 -u PSI_non_dim_full_form_vis