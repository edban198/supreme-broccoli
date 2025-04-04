#!/bin/bash
#SBATCH --job-name=RB_sim
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 2-00:00:00
#SBATCH --array=0-728%64
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sfbj55@durham.ac.uk

Prs=(1.0 1.048 1.099 1.151 1.207 1.265 1.326 1.389 1.456 1.526 1.6 1.677 1.758 1.842 1.931 2.024 2.121 2.223 2.33 2.442 2.56 2.683 2.812 2.947 3.089 3.237 3.393 3.556 3.728 3.907 4.095 4.292 4.498 4.715 4.942 5.179 5.429 5.69 5.964 6.251 6.551 6.866 7.197 7.543 7.906 8.286 8.685 9.103 9.541 10.0)

Chis=(908.6 1046.1 1204.5 1386.8 1596.8 1838.5 2116.9 2437.4 2806.4 3231.3 3720.5 4283.7 4932.2 5679.0 6538.7 7528.7 8668.5 9980.8 11491.9 13231.7 15234.9 17541.4 20197.1 23254.9 26775.6 30829.3 35496.7 40870.7 47058.3 54182.7 62385.8 71830.7 82705.5 95226.7 109643.6 126243.2 145355.8 167362.0 192699.8 221873.7 255464.4 294140.5 338672.0 389945.4 448981.4 516955.1 595219.7 685333.2 789089.5 908554.0)

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