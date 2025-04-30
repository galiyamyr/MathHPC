#!/bin/bash
#SBATCH --job-name=timing_test
#SBATCH --output=timing_output_%j.txt
#SBATCH --error=timing_error_%j.txt
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --partition=instruction

module load mpi
mpicc burg3.c -o burgers_mpi -lm

# Run with different numbers of processors
for np in 1 2 4 8
do
    echo "Running with $np processes:"
    mpirun -np $np ./burgers_mpi
    echo "----------------------------"
done

