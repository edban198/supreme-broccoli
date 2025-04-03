#!/bin/bash
#SBATCH --job-name=RB_sim
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 2-00:00:00
#SBATCH --array=0-4%5     # Array index matches number of Pr values (or use 0-215%64 later)
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sfbj55@durham.ac.uk

# Define Prandtl numbers
Prs=(2 4 6 8 10)

# Get index from SLURM
index=$SLURM_ARRAY_TASK_ID

# Extract Pr value
Pr_val=${Prs[$index]}

# Fixed chi value
Chi_val=10

echo "Running simulation with: Pr = $Pr_val, chi = $Chi_val"

# Run the Julia script
~/julia-1.11.2/bin/julia ~/CODE/supreme-broccoli/array_RB.jl $Pr_val $Chi_val