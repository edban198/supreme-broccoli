#!/bin/bash

# Instructing SLURM to locate and assign
#X number of nodes with Y number of
#cores in each node.
# X,Y are integers. Refer to table for various combinations
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 02:30:00

#SBATCH --job-name=wave_propagation_simulation

#SBATCH --mail-type=ALL
#SBATCH --mail-user sfbj55@durham.ac.uk

# Run the program
cd

./julia-1.11.2/bin/julia ./CODE/supreme-broccoli/wave_propagation.jl