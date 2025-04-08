#!/bin/bash
#SBATCH --job-name=RB_sim_wind
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p cpu
#SBATCH --qos=short
#SBATCH -t 2-00:00:00
#SBATCH --array=1-320%64
#SBATCH -o RB_sim_wind_%A_%a.out
#SBATCH -e RB_sim_wind_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sfbj55@durham.ac.uk

# Define the simulation parameter arrays
Prs=(1.0 1.0771 1.1602 1.2496 1.346 1.4497 1.5615 1.6819 1.8116 1.9513 2.1017 2.2638 2.4384 2.6264 2.8289 3.047 3.2819 3.535 3.8075 4.1011 4.4173 4.7579 5.1248 5.52 5.9456 6.404 6.8978 7.4296 8.0025 8.6195 9.2841 10.0)
CHIS=(2 5 10 15 20)
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
~/julia-1.11.2/bin/julia --threads=$SLURM_CPUS_PER_TASK ~/CODE/supreme-broccoli/array_RB_wind.jl $Pr_val $Chi_val $Wind_val