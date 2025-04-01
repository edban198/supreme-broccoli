#!/bin/bash
#SBATCH --job-name=RB_cpu_simulation_array
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 18:00:00
#SBATCH --array=1-266
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=sfbj55@durham.ac.uk

cd

# Parameter values to simulate
Prs=(1 1.25 1.5 1.75 2 3 4 5 6 8 10 12 15 20 30 40 60 80 100)19
Chis=(1 1.1 1.2 1.3 1.5 1.75 2 3 4 5 7.5 10 20 50)14

N_Prs=${#Prs[@]}
N_Chis=${#Chis[@]}

index=$((SLURM_ARRAY_TASK_ID - 1))

Prs_index=$(( index % N_Prs ))
Chis_index=$(( index / N_Prs ))

Prs_vals=${Prs[$Prs_index]}
Chis_vals=${Chis[$Chis_index]}

echo "Running simulation for Prs: $Prs_vals, Chis: $Chis_vals"

# Run the Julia simulation
~/julia-1.11.2/bin/julia ./CODE/supreme-broccoli/array_RB.jl $Prs_vals $Chis_vals