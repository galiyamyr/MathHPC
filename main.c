#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, size;
    double a, b;
    int N;
    
    // Check command-line arguments
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <a> <b> <N>\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    a = atof(argv[1]);
    b = atof(argv[2]);
    N = atoi(argv[3]);

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Simulated computation for T and Tex (use actual logic here)
    double T = 0.0, Tex = 1.0; // Replace with actual calculations
    for (int i = rank; i < N; i += size) {
        // Example computation (replace with actual logic)
        T += exp(-a * i + b);
    }

    // Reduce all T values from different processes to process 0
    double global_T;
    MPI_Reduce(&T, &global_T, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Root process (rank 0) prints the result
    if (rank == 0) {
        // Calculate relative error (|T-Tex|/|Tex|)
        double error = fabs(global_T - Tex) / Tex;
        printf("NP = %d, N = %d, T = %.13e, |T-Tex|/|Tex| = %.13e\n", size, N, global_T, error);

        // Record the start time of the computation
        double start_time = MPI_Wtime();
        
        // Do actual computations (your task here)

        double end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;

        // Print elapsed time
        printf("Elapsed time = %.13e\n", elapsed_time);
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
