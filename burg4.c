#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define PI 3.141592653589793

// Simulation constants
const double nu = 0.01; // Viscosity
const double L = 1.0;   // Length of the domain
const double T = 1.0;   // Total simulation time

int main(int argc, char *argv[]) {
    int rank, size;
    int N = 100; // Total grid points
    double dx = L / (N - 1);
    double dt = 0.4 * dx * dx / nu;
    int nsteps = (int)(T / dt);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (N % size != 0) {
        if (rank == 0) printf("N must be divisible by number of processes.\n");
        MPI_Finalize();
        return 1;
    }

    int local_N = N / size; // Points per process
    int local_start = rank * local_N;

    // Add 2 ghost cells
    double *u = (double *)malloc((local_N + 2) * sizeof(double));
    double *unew = (double *)malloc((local_N + 2) * sizeof(double));

    // Initial condition: sin(πx)
    for (int i = 1; i <= local_N; i++) {
        double x = (local_start + i - 1) * dx;
        u[i] = sin(PI * x);
    }

    u[0] = 0.0;                // Left ghost (placeholder)
    u[local_N + 1] = 0.0;      // Right ghost (placeholder)

    double start_time = MPI_Wtime();

    // Time stepping loop
    for (int step = 0; step < nsteps; step++) {
        // Halo exchange
        if (rank > 0) {
            MPI_Sendrecv(&u[1], 1, MPI_DOUBLE, rank - 1, 0,
                         &u[0], 1, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Sendrecv(&u[local_N], 1, MPI_DOUBLE, rank + 1, 0,
                         &u[local_N + 1], 1, MPI_DOUBLE, rank + 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Update internal points
        for (int i = 1; i <= local_N; i++) {
            unew[i] = u[i] - dt / dx * u[i] * (u[i] - u[i - 1])
                      + nu * dt / (dx * dx) * (u[i + 1] - 2 * u[i] + u[i - 1]);
        }

        // Swap pointers
        double *temp = u;
        u = unew;
        unew = temp;
    }

    double end_time = MPI_Wtime();

    if (rank == 0) {
        printf("Elapsed time with %d processes: %f seconds\n", size, end_time - start_time);
    }

    // Gather and write results
    double *global_u = NULL;
    if (rank == 0) {
        global_u = (double *)malloc(N * sizeof(double));
    }

    MPI_Gather(&u[1], local_N, MPI_DOUBLE, global_u, local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        FILE *fout = fopen("u_final.txt", "w");
        for (int i = 0; i < N; i++) {
            fprintf(fout, "%f %f\n", i * dx, global_u[i]);
        }
        fclose(fout);
        free(global_u);
    }

    free(u);
    free(unew);
    MPI_Finalize();
    return 0;
}

