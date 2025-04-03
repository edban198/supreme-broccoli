#!/bin/bash
#SBATCH --job-name=RB_sim
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 2-00:00:00
#SBATCH --array=0-809%64        # 800 jobs, up to 64 running at once
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sfbj55@durham.ac.uk

Prs=(2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9
     3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9
     4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9
     5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9
     6.0 6.1 6.2 6.3 6.4 6.5 6.6 6.7 6.8 6.9
     7.0 7.1 7.2 7.3 7.4 7.5 7.6 7.7 7.8 7.9
     8.0 8.1 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.9
     9.0 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9
     10.0)

Chis=(5 6 7 8 9 10 11 12 13 14)

# Compute array indices
nPr=${#Prs[@]}      # 80
nChi=${#Chis[@]}    # 10

index=$SLURM_ARRAY_TASK_ID
iPr=$(( index % nPr ))
iChi=$(( index / nPr ))

# Get parameter values
Pr_val=${Prs[$iPr]}
Chi_val=${Chis[$iChi]}

echo "Running simulation with: Pr = $Pr_val, Chi = $Chi_val"

# Run Julia script
~/julia-1.11.2/bin/julia ~/CODE/supreme-broccoli/array_RB.jl $Pr_val $Chi_val