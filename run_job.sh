#!/bin/bash
# This line is required to inform the Linux
#command line to parse the script using
#the bash shell

# Instructing SLURM to locate and assign
#X number of nodes with Y number of
#cores in each node.
# X,Y are integers. Refer to table for
#various combinations
#SBATCH -N X
#SBATCH -c Y

# Governs the run time limit and
# resource limit for the job. Please pick values
# from the partition and QOS tables below
#for various combinations
#SBATCH -p "partitionname"
#SBATCH --qos="qosname"
#SBATCH -t DD-HH:MM:SS

# Source the bash profile (required to use the module command)
source /etc/profile


# Run your program (replace this with your program)
./hello-world