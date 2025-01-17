#!/bin/bash

# Instructing SLURM to locate and assign
#X number of nodes with Y number of
#cores in each node.
# X,Y are integers. Refer to table for various combinations
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --gres=gpu:pascal:1
#SBATCH -p ug-gpu-small
#SBATCH --qos=short
#SBATCH -t 01:00:00

#SBATCH --job-name=RB_gpu_simulation

#SBATCH --mail-type=ALL
#SBATCH --mail-user sfbj55@durham.ac.uk

# Run the program
./julia-1.11.2/bin/julia ./RB_with_wind_2.jl