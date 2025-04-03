#!/bin/bash
#SBATCH --job-name=RB_sim
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 2-00:00:00
#SBATCH --array=0-215%64
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sfbj55@durham.ac.uk

Prs=(1 1.5 2 3 4 5 6 8 10 12 15 20 30 40 50 60 80 100)
Chis=(1 1.25 1.5 1.75 2 3 4 5 8 10 20 50)

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