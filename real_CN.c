#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define PI 3.141592653589793

const double nu = 0.01;  // Viscosity
const double L = 1.0;    // Domain length
const double T = 1.0;    // Final time
const int N = 500;       // Total number of grid points

void initialize(double *u, int local_N, int rank, double dx) {
    for (int i = 0; i < local_N; i++) {
        int global_i = rank * (local_N - 2) + i;
        u[i] = sin(PI * global_i * dx);
    }
}

void solve_tridiagonal(int n, double *a, double *b, double *c, double *d, double *x) {
    for (int i = 1; i < n; i++) {
        double m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }

    x[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int local_N = N / size + 2; // 2 ghost cells
    double dx = L / (N - 1);
    double dt = 0.4 * dx * dx / nu;
    int nsteps = (int)(T / dt);

    double *u = malloc(local_N * sizeof(double));
    double *u_new = malloc(local_N * sizeof(double));
    double *rhs = malloc((local_N - 2) * sizeof(double));

    double *a = malloc((local_N - 2) * sizeof(double));
    double *b = malloc((local_N - 2) * sizeof(double));
    double *c = malloc((local_N - 2) * sizeof(double));

    initialize(u, local_N, rank, dx);

    double start_time = MPI_Wtime();

    for (int step = 0; step < nsteps; step++) {
        // Exchange ghost cells
        if (rank > 0)
            MPI_Sendrecv(&u[1], 1, MPI_DOUBLE, rank - 1, 0, &u[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rank < size - 1)
            MPI_Sendrecv(&u[local_N - 2], 1, MPI_DOUBLE, rank + 1, 0, &u[local_N - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Build RHS using explicit convection and Crank-Nicolson diffusion
        for (int i = 1; i < local_N - 1; i++) {
            double conv = -0.5 * dt / dx * u[i] * (u[i + 1] - u[i - 1]);
            double diff = 0.5 * nu * dt / (dx * dx) * (u[i + 1] - 2 * u[i] + u[i - 1]);
            rhs[i - 1] = u[i] + conv + diff;

            a[i - 1] = -0.5 * nu * dt / (dx * dx);
            b[i - 1] = 1 + nu * dt / (dx * dx);
            c[i - 1] = -0.5 * nu * dt / (dx * dx);
        }

        // Solve local tridiagonal system
        solve_tridiagonal(local_N - 2, a, b, c, rhs, &u_new[1]);

        // Apply boundary values (Dirichlet)
        u_new[0] = 0.0;
        u_new[local_N - 1] = 0.0;

        // Update
        for (int i = 0; i < local_N; i++)
            u[i] = u_new[i];
    }

    double end_time = MPI_Wtime();

    if (rank == 0) {
        printf("Elapsed time with %d processes: %f seconds\n", size, end_time - start_time);
    }

    // Gather all data to root
    double *global_u = NULL;
    if (rank == 0)
        global_u = malloc(N * sizeof(double));

    MPI_Gather(&u[1], local_N - 2, MPI_DOUBLE, global_u, local_N - 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        FILE *fout = fopen("u_cn_final.txt", "w");
        for (int i = 0; i < N; i++) {
            fprintf(fout, "%f %f\n", i * dx, global_u[i]);
        }
        fclose(fout);
        free(global_u);
    }

    free(u);
    free(u_new);
    free(rhs);
    free(a);
    free(b);
    free(c);

    MPI_Finalize();
    return 0;
}

