#!/bin/bash
#SBATCH --job-name=mpi_job        # Job name
#SBATCH --output=mpi_job_output.log  # Output log file
#SBATCH --ntasks=32               # Number of tasks (processes)
#SBATCH --time=01:00:00           # Time limit
#SBATCH --partition=instruction    # Partition (adjust based on your system)
#SBATCH --mem=120GB               # Memory allocation per node
#SBATCH --nodes=1                 # Number of nodes (1 node in this case)

# Load MPI module
module load openmpi

# Run the MPI job
mpirun -np 32 ./mpi_job 1.0 2.0 1024
