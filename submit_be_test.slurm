#!/bin/bash
#SBATCH --job-name=be_burgers
#SBATCH --output=be_burgers_output.txt
#SBATCH --time=00:10:00
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --partition=instruction

module purge
module load intel/22.3.1
module load mpi/2021.7.1

echo "Compiling..."
mpicc -O3 -o burgers_be burgers_be.c -lm

# Grid size = 100000, time steps = 2000
for procs in 1 2 4 8; do
  echo "Running with $procs processors"
  mpirun -np $procs ./burgers_be 100000 2000
done

