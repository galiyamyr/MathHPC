#!/bin/bash
#SBATCH --job-name=burg_scaling
#SBATCH --output=timing_output_%j.txt
#SBATCH --error=timing_error_%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:10:00
#SBATCH --partition=instruction


module load gcc
module load openmpi

echo "Running with 1 process:"
srun -n 1 ./burg_cn
echo "----------------------------"

echo "Running with 2 processes:"
srun -n 2 ./burg_cn
echo "----------------------------"

echo "Running with 4 processes:"
srun -n 4 ./burg_cn
echo "----------------------------"

echo "Running with 8 processes:"
srun -n 8 ./burg_cn
echo "----------------------------"

