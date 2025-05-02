#!/bin/bash
#SBATCH --job-name=burgers_mpi       # Name of the job
#SBATCH --output=burgers_mpi_%j.txt  # Output file name with job ID
#SBATCH --ntasks=4                   # Number of MPI tasks (4 processors)
#SBATCH --time=00:30:00              # Maximum runtime (e.g., 30 minutes)
#SBATCH --partition=instruction         # SLURM partition (adjust if needed)

# Load the necessary modules
module load mpi/2021.13   # Adjust based on the MPI version available
module load intel/2024.2  # Adjust based on the compiler version available

# Compile the C program
mpicc burg_attempt2.c -o burgers_mpi.out -lm

# Run the program with 4 processes
mpirun -np 4 ./burgers_mpi.out

