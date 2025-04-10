#!/bin/bash
#SBATCH --job-name=RB_sim
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 2-00:00:00
#SBATCH --array=1-128%16
#SBATCH -o RB_sim_%A_%a.out
#SBATCH -e RB_sim_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sfbj55@durham.ac.uk

# Define the simulation parameter arrays
Prs=(1.0 1.1659 1.3594 1.5849 1.8478 2.1544 2.5119 2.9286 3.4145 3.9811 4.6416 5.4117 6.3096 7.3564 8.577 10.0)
CHIS=(25 30 40 50)
# Define wind forcing parameters (change these as desired)
WINDS=(0.0 0.00001603125)

# Get the number of values in each array
N_Prs=${#Prs[@]}
N_Chis=${#CHIS[@]}
N_Winds=${#WINDS[@]}

# Total number of tasks equals all possible combinations of Pr, Chi and Wind.
total_tasks=$(( N_Prs * N_Chis * N_Winds ))
echo "Total number of tasks = $total_tasks"

# Convert SLURM_ARRAY_TASK_ID (starting at 1) to a zero-indexed variable.
index=$(( SLURM_ARRAY_TASK_ID - 1 ))

# Determine the index for each parameter.
Pr_index=$(( index % N_Prs ))
Chi_index=$(( (index / N_Prs) % N_Chis ))
Wind_index=$(( index / (N_Prs * N_Chis) ))

# Get the parameter values.
Pr_val=${Prs[$Pr_index]}
Chi_val=${CHIS[$Chi_index]}
Wind_val=${WINDS[$Wind_index]}

# Print the current task information.
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "Zero-indexed task: $index"
echo "Pr_index: $Pr_index, Chi_index: $Chi_index, Wind_index: $Wind_index"
echo "Running simulation with: Pr = $Pr_val, chi = $Chi_val, wind forcing = $Wind_val"

# Run the Julia simulation with the extra wind forcing parameter.
# (Make sure that your Julia script is modified to accept the third argument.)
~/julia-1.11.2/bin/julia --threads=$SLURM_CPUS_PER_TASK ~/CODE/supreme-broccoli/array_RB.jl $Pr_val $Chi_val $Wind_val