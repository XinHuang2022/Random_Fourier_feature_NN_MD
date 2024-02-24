#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A naiss2023-22-1095

# The name of the script is myjob
#SBATCH -J Correlation_tau_ref

# The partition
#SBATCH -p main

# The number of cores requested
#SBATCH -n 128

# Number of hours wall-clock time will be given to this job
#SBATCH -t 04:20:00

# Load the Matlab module
ml add PDC/22.06
ml matlab/r2023a

# Run matlab taking your_matlab_program.m as input and show the output in the file
# your_matlab_program.out. The input file must be in the directory where you submit this script. 
# This is also where the output will be created.

matlab -nosplash -nodesktop -nodisplay < Tester_approx_V_md_corr.m > Tester_correlation_ref.out