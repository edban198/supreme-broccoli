#!/bin/bash
#SBATCH --job-name=RB_sim
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 1-00:00:00
#SBATCH --array=0-186    # Adjust depending on total combinations
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err

# Define arrays (NO COMMAS, NO QUOTES)
Prs=(1 1.5 2 3 4 5 6 8 10 12 15 20 30 40 60 80 100)
Chis=(1 1.2 1.5 2 3 4 5 7.5 10 20 50)

N_Prs=${#Prs[@]}
N_Chis=${#Chis[@]}

index=$SLURM_ARRAY_TASK_ID
Pr_index=$((index % N_Prs))
Chi_index=$((index / N_Prs))

Pr_val=${Prs[$Pr_index]}
Chi_val=${Chis[$Chi_index]}

echo "Running: Pr = $Pr_val, chi = $Chi_val"

# Run your Julia script WITHOUT quotes
~/julia-1.11.2/bin/julia ~/CODE/supreme-broccoli/array_RB.jl $Pr_val $Chi_val