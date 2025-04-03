#!/bin/bash
#SBATCH --job-name=RB_sim
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 2-00:00:00
#SBATCH --array=0-63     # Array index matches number of Pr values (or use 0-215%64 later)
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sfbj55@durham.ac.uk

# Define Prandtl numbers
Prs=(2.0 2.127 2.254 2.381 2.508 2.635 2.762 2.889 3.016 3.143 3.27 3.397 3.524 3.651 3.778 3.905 4.032 4.159 4.286 4.413 4.54 4.667 4.794 4.921 5.048 5.175 5.302 5.429 5.556 5.683 5.81 5.937 6.064 6.19 6.317 6.444 6.571 6.698 6.825 6.952 7.079 7.206 7.333 7.46 7.587 7.714 7.841 7.968 8.095 8.222 8.349 8.476 8.603 8.73 8.857 8.984 9.111 9.238 9.365 9.492 9.619 9.746 9.873 10.0)

# Get index from SLURM
index=$SLURM_ARRAY_TASK_ID

# Extract Pr value
Pr_val=${Prs[$index]}

# Fixed chi value
Chi_val=10

echo "Running simulation with: Pr = $Pr_val, chi = $Chi_val"

# Run the Julia script
~/julia-1.11.2/bin/julia ~/CODE/supreme-broccoli/array_RB.jl $Pr_val $Chi_val