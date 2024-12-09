#!/bin/bash

# Instructing SLURM to locate and assign
#X number of nodes with Y number of
#cores in each node.
# X,Y are integers. Refer to table for various combinations
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --gres=gpu
#SBATCH -p gpu-small
#SBATCH --qos=short
#SBATCH --job-name=RB_gpu_simulation
#SBATCH --mail-type=ALL
#SBATCH --mail-user sfbj55@durham.ac.uk

# Source the bash profile (required to use the module command)
source /etc/profile

module load julia/1.9.2

# Run your program (replace this with your program)
julia RB.jl

# add log file to log the progress messages
julia RB.jl > output.log 2>&1