#!/bin/bash

# Instructing SLURM to locate and assign
#X number of nodes with Y number of
#cores in each node.
# X,Y are integers. Refer to table for various combinations
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 08:00:00

#SBATCH --job-name=white_noise_cpu_sim_small

#SBATCH --mail-type=ALL
#SBATCH --mail-user sfbj55@durham.ac.uk

# Run the program
cd

./julia-1.11.2/bin/julia ./CODE/supreme-broccoli/RB_with_wind_3.jl