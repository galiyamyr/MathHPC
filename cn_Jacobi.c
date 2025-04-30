#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define PI 3.141592653589793

// Constants for the simulation
const double nu = 0.01;      // Viscosity
const double L = 1.0;        // Length of the domain
const double T = 1.0;        // Total time for the simulation
const int N = 5000;           // Number of spatial grid points

// Function to update the solution using the Jacobi method
void jacobi_update(double *u, double *unew, double *u_left, double *u_right, double dx, double dt, double nu, int local_N) {
    for (int i = 1; i < local_N - 1; i++) {
        unew[i] = u[i] - dt / (2 * dx) * u[i] * (u[i + 1] - u[i - 1]) + (nu * dt) / (dx * dx) * (u[i + 1] - 2 * u[i] + u[i - 1]);
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    int local_N;              // Local grid size per processor
    double dx, dt, start_time, end_time;
    int nsteps;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Setup grid and time step
    dx = L / (N - 1);               // Grid spacing
    dt = 0.4 * dx * dx / nu;        // Time step based on the CFL condition
    nsteps = (int)(T / dt);         // Number of time steps
    local_N = N / size;             // Local grid size for each processor

    // Allocate memory for the solution arrays
    double *u = (double *)malloc(local_N * sizeof(double));      // Current solution
    double *unew = (double *)malloc(local_N * sizeof(double));  // Next solution
    double *u_left = (double *)malloc(sizeof(double));          // Left boundary
    double *u_right = (double *)malloc(sizeof(double));         // Right boundary

    // Initialize the solution u(x, 0) = sin(πx)
    for (int i = 0; i < local_N; i++) {
        u[i] = sin(PI * (rank * local_N + i) * dx);
    }

    // Start timing the computation
    start_time = MPI_Wtime();

    // Perform the time-stepping loop
    for (int step = 0; step < nsteps; step++) {
        // Communicate boundary values with neighboring processors
        if (rank > 0) {
            MPI_Send(&u[1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);   // Send left boundary
            MPI_Recv(u_left, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // Receive left boundary
        }
        if (rank < size - 1) {
            MPI_Send(&u[local_N - 2], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);  // Send right boundary
            MPI_Recv(u_right, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  // Receive right boundary
        }

        // Update the solution using Jacobi method
        jacobi_update(u, unew, u_left, u_right, dx, dt, nu, local_N);

        // Swap the arrays (u and unew) for the next iteration
        double *temp = u;
        u = unew;
        unew = temp;
    }

    // End timing the computation
    end_time = MPI_Wtime();

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
    free(u_left);
    free(u_right);

    // Finalize MPI
    MPI_Finalize();
    return 0;
}

