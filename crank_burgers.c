#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define IDX(i) ((i) - start + 1)  // local index (with ghost cells)

// Function to initialize the initial condition
double u0(double x) {
    return sin(M_PI * x);  // Example initial condition
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const double L = 1.0;
    const int N = 100;        // total number of spatial grid points
    const double nu = 0.01;   // viscosity
    const double T = 1.0;     // final time
    const double dt = 0.001;
    const int Nt = T / dt;
    const double dx = L / (N - 1);

    int local_N = N / size;
    int remainder = N % size;

    if (rank < remainder)
        local_N++;

    int start = rank * (N / size) + (rank < remainder ? rank : remainder);
    int end = start + local_N - 1;

    // Allocate arrays with ghost cells
    double *u = calloc(local_N + 2, sizeof(double));
    double *unew = calloc(local_N + 2, sizeof(double));
    double *rhs = calloc(local_N + 2, sizeof(double));

    // Initialize
    for (int i = start; i <= end; ++i) {
        double x = i * dx;
        u[IDX(i)] = u0(x);
    }

    // Time loop
    for (int n = 0; n < Nt; ++n) {
        // Exchange boundary data
        if (rank > 0)
            MPI_Sendrecv(&u[1], 1, MPI_DOUBLE, rank - 1, 0,
                         &u[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        else
            u[0] = 0.0;

        if (rank < size - 1)
            MPI_Sendrecv(&u[local_N], 1, MPI_DOUBLE, rank + 1, 0,
                         &u[local_N + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        else
            u[local_N + 1] = 0.0;

        // Crank-Nicolson RHS
        for (int i = start + 1; i < end; ++i) {
            int li = IDX(i);
            double u_x = (u[li + 1] - u[li - 1]) / (2 * dx);
            double u_xx = (u[li + 1] - 2 * u[li] + u[li - 1]) / (dx * dx);
            rhs[li] = u[li] - 0.5 * dt * (u[li] * u_x) + 0.5 * dt * nu * u_xx;
        }

        // Solve tridiagonal system for Crank-Nicolson
        // For simplicity, using Jacobi iteration (can be replaced with Thomas or other solvers)
        for (int it = 0; it < 20; ++it) {
            for (int i = start + 1; i < end; ++i) {
                int li = IDX(i);
                double u_x = (u[li + 1] - u[li - 1]) / (2 * dx);
                double a = -nu * dt / (2 * dx * dx);
                double b = 1.0 + nu * dt / (dx * dx);
                double c = -nu * dt / (2 * dx * dx);
                unew[li] = (rhs[li] - a * u[li - 1] - c * u[li + 1]) / b;
            }
            // Swap pointers
            double *tmp = u;
            u = unew;
            unew = tmp;
        }
    }

    // Gather data to rank 0
    int *recvcounts = NULL, *displs = NULL;
    if (rank == 0) {
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
    }

    int sendcount = local_N;
    MPI_Gather(&sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; ++i)
            displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    double *global_u = NULL;
    if (rank == 0)
        global_u = malloc(N * sizeof(double));

    MPI_Gatherv(&u[1], local_N, MPI_DOUBLE, global_u, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Write output
    if (rank == 0) {
        FILE *fout = fopen("solution.txt", "w");
        for (int i = 0; i < N; ++i)
            fprintf(fout, "%f %f\n", i * dx, global_u[i]);
        fclose(fout);
        free(global_u);
        free(recvcounts);
        free(displs);
    }

    free(u);
    free(unew);
    free(rhs);

    MPI_Finalize();
    return 0;
}

