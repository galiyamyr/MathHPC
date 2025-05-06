#!/bin/bash
#SBATCH --job-name=burgers_benchmark
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
#SBATCH --partition=instruction

# === Configuration ===
SCALING_MODE="weak"   # Change to "strong" for strong scaling
BASE_N=1000           # N per process (weak) or total N (strong)

# === Load Modules ===
module load intel/2024.2
module load mpi/2021.13

# === Compile C Code ===
mpicc -O3 -o burgers_solver test86.c -lm
if [ ! -f burgers_solver ]; then
    echo "Compilation failed!"
    exit 1
fi

# === Log Info ===
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Execution times for ${SCALING_MODE} scaling" > execution_times.txt
echo "Job ID: $SLURM_JOB_ID" >> execution_times.txt
echo "-------------------------------------" >> execution_times.txt

# === Loop Over Processor Counts ===
for procs in 1 2 4 8; do
    if [ "$SCALING_MODE" = "weak" ]; then
        N=$((procs * BASE_N))
    else
        N=$BASE_N
    fi
    echo "Running with $procs process(es), N=$N..." >> execution_times.txt
    { time mpirun -np $procs ./burgers_solver $N; } 2>> execution_times.txt
    echo "-------------------------------------" >> execution_times.txt
done


