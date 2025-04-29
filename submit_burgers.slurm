#!/bin/bash
#SBATCH --job-name=burgers_mpi
#SBATCH --output=burgers_output.txt
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
#SBATCH --partition=instruction
#SBATCH --mail-type=END

module purge
module load hdf5-intel/1.14.2
module load mpi/2021.7.1

mpicc crank_burgers.c -o crank_burgers -lm
mpirun -np 4 ./crank_burgers

