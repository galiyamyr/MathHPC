#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#define PI 3.141592653589793

// Constants for the simulation
const double nu = 0.01; // Viscosity
const double L = 1.0;    // Length of the domain
const double T = 1.0;    // Total time for the simulation

int main(int argc, char *argv[]) {
    int rank, size;
    int N = 100;            // Number of spatial grid points
    double dx = L / (N - 1); // Grid spacing
    double dt = 0.4 * dx * dx / nu; // Time step based on the CFL condition
    int nsteps = (int)(T / dt);  // Number of time steps

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Debugging start
    if (rank == 0) {
        printf("Starting the computation with %d processes...\n", size);
    }

    // Start timing the computation
    double start_time = MPI_Wtime(); // Start time

    // Allocate memory for the solution arrays
    double *u = (double *)malloc(N * sizeof(double)); // Solution at current time step
    double *unew = (double *)malloc(N * sizeof(double)); // Solution at next time step

    // Initialize the solution u(x, 0) = sin(πx)
    for (int i = 0; i < N; i++) {
        u[i] = sin(PI * i * dx);
    }

    // Perform the time-stepping loop
    for (int step = 0; step < nsteps; step++) {
        // Time-stepping: Explicit scheme for Burgers' equation
        for (int i = 1; i < N-1; i++) {
            unew[i] = u[i] - dt / dx * u[i] * (u[i] - u[i-1]) + nu * dt / (dx*dx) * (u[i+1] - 2*u[i] + u[i-1]);
        }

        // Copy the new values back to u
        for (int i = 1; i < N-1; i++) {
            u[i] = unew[i];
        }
    }

    // End timing the computation
    double end_time = MPI_Wtime();   // End time

    // Print the elapsed time for rank 0 (the root process)
    if (rank == 0) {
        printf("Elapsed time with %d processes: %f seconds\n", size, end_time - start_time);
    }

    // Save the results to a file (only rank 0 writes)
    if (rank == 0) {
        FILE *fout = fopen("u_final.txt", "w");
        for (int i = 0; i < N; i++) {
            fprintf(fout, "%f %f\n", i * dx, u[i]);
        }
        fclose(fout);
    }

    // Free allocated memory
    free(u);
    free(unew);

    // Finalize MPI
    MPI_Finalize();
    return 0;
}

