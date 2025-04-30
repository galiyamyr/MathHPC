#!/bin/bash
#SBATCH --job-name=burgers_timing
#SBATCH --output=timing_output_%j.txt
#SBATCH --ntasks=8
#SBATCH --time=00:02:00
#SBATCH --partition=instruction

module load mpi
mpicc burg4.c -o burgers_explicit.out -lm

for procs in 1 2 4 8
do
    echo -n "Time for $procs procs: "
    mpirun -np $procs ./burgers_explicit.out | grep "Elapsed time" || echo "Failed or no output"
done

