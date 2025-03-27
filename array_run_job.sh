#!/bin/bash
#SBATCH --job-name=RB_gpu_simulation_array
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=14G
#SBATCH --gres=gpu:pascal:1
#SBATCH -p ug-gpu-small
#SBATCH --qos=short
#SBATCH -t 00-12:00:00
#SBATCH --array=1-14
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=sfbj55@durham.ac.uk

cd

# Gamma values to simulate
GAMMAS=(1 1.1 1.2 1.4 2 3 4 6 10 15 20 30 40 50)

# Get the gamma for this task
GAMMA=${GAMMAS[$SLURM_ARRAY_TASK_ID - 1]}

# Run the Julia script with this gamma
./julia-1.11.2/bin/julia ./CODE/supreme-broccoli/array_RB.jl $GAMMA
