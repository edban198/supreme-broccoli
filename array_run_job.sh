#!/bin/bash
#SBATCH --job-name=RB_cpu_simulation_array
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 18:00:00
#SBATCH --array=1-14
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=sfbj55@durham.ac.uk


cd

# Gamma values to simulate
Prs=(1 2 4 6 8 10 16 20 40 50 60 70 80 100)

# Get the gamma for this task
Prs=${Prs[$SLURM_ARRAY_TASK_ID - 1]}

# Run the Julia simulation
~/julia-1.11.2/bin/julia ./CODE/supreme-broccoli/array_RB.jl $Prs