#!/bin/bash
#SBATCH --job-name=RB_cpu_simulation_array
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 18:00:00
#SBATCH --array=1-14
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=sfbj55@durham.ac.uk

# Optional: ensure you’re in the correct directory
cd ~/CODE/supreme-broccoli

# Gamma values to simulate
GAMMAS=(1 1.1 1.2 1.4 2 3 4 6 10 15 20 30 40 50)

# Get the gamma for this task
GAMMA=${GAMMAS[$SLURM_ARRAY_TASK_ID - 1]}

# Run the Julia simulation
~/julia-1.11.2/bin/julia ./array_RB.jl $GAMMA

# Extract the "data for csv" line and append to master CSV file
grep "data for csv:" RB_sim_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err | \
    sed 's/.*data for csv: //' >> RB_results_${SLURM_JOB_ID}.csv