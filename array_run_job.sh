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

Prs=(1 7)
CHIS=(1 1.1 1.2 1.4 2 3 4 6 10 15 20 30 40 50)

N_Prs=${#Prs[@]}
N_Chis=${#Chis[@]}

index=$SLURM_ARRAY_TASK_ID
Pr_index=$((index % N_Prs))
Chi_index=$((index / N_Prs))

Pr_val=${Prs[$Pr_index]}
Chi_val=${Chis[$Chi_index]}

echo "Running: Pr = $Pr_val, chi = $Chi_val"

~/julia-1.11.2/bin/julia --threads=$SLURM_CPUS_PER_TASK ~/CODE/supreme-broccoli/array_RB.jl $Pr_val $Chi_val