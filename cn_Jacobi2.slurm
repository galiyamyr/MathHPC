#!/bin/bash
#SBATCH --job-name=burgers_cn_mpi
#SBATCH --output=cn_output_%j.txt
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --partition=instruction

module purge
module load gcc
module load openmpi

echo "Running with 1 process:"
srun --ntasks=1 ./burgers_cn_mpi

echo "----------------------------"

echo "Running with 2 processes:"
srun --ntasks=2 ./burgers_cn_mpi

echo "----------------------------"

echo "Running with 4 processes:"
srun --ntasks=4 ./burgers_cn_mpi

echo "----------------------------"

echo "Running with 8 processes:"
srun --ntasks=8 ./burgers_cn_mpi

