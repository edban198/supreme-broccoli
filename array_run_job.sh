#!/bin/bash
#SBATCH --job-name=RB_sim
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 2-00:00:00
#SBATCH --array=1-33
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sfbj55@durham.ac.uk

# Get line number based on SLURM_ARRAY_TASK_ID
line_number=$SLURM_ARRAY_TASK_ID

# Read Pr and R from the corresponding line in the CSV (skip header)
read R_val Pr_val <<< $(awk -F, -v line=$line_number 'NR == line + 1 { print $1, $2 }' ~/CODE/supreme-broccoli/points_in_region_I.csv)

echo "Running simulation with: Pr = $Pr_val, R = $R_val"

# Run Julia script (now expects Pr and R)
~/julia-1.11.2/bin/julia ~/CODE/supreme-broccoli/array_RB.jl $Pr_val $R_val