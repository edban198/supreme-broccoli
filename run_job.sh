#!/bin/bash

# Instructing SLURM to locate and assign
#X number of nodes with Y number of
#cores in each node.
# X,Y are integers. Refer to table for various combinations
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --gres=gpu
#SBATCH -p ug-gpu-small
#SBATCH --qos=short

#SBATCH --job-name=RB_gpu_simulation

#SBATCH --mail-type=ALL
#SBATCH --mail-user sfbj55@durham.ac.uk

# Source the bash profile (required to use the module command)
source /etc/profile

module load julia/1.9.2
# Run the program
julia RB.jl

# delete tempory file
rm -f RB_gpu_simulation.jld2