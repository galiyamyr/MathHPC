#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define N 100
#define L 1.0
#define T 1.0
#define nu 0.01
#define dt 0.001
#define PI 3.141592653589793

double u[N+1], u_new[N+1];
double a[N-1], b[N-1], c[N-1], d[N-1]; // d is now globally allocated

void initialize() {
    for (int i = 0; i <= N; i++) {
        double x = i * L / N;
        u[i] = sin(PI * x);
    }
}

void build_matrix() {
    double dx = L / N;
    double r = nu * dt / (2 * dx * dx);
    for (int i = 0; i < N - 1; i++) {
        a[i] = -r;
        b[i] = 1 + 2 * r;
        c[i] = -r;
    }
}

void solve_tridiagonal() {
    int n = N - 1;

    // Forward elimination
    for (int i = 1; i < n; i++) {
        double m = a[i] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }

    // Back substitution
    u_new[N - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        u_new[i + 1] = (d[i] - c[i] * u_new[i + 2]) / b[i];
    }

    u_new[0] = 0.0;
    u_new[N] = 0.0;
}

void write_output() {
    FILE *fp = fopen("burgers_output.txt", "w");
    double dx = L / N;
    for (int i = 0; i <= N; i++) {
        fprintf(fp, "%f %f\n", i * dx, u[i]);
    }
    fclose(fp);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int steps = T / dt;
    int base = (N - 1) / size;
    int remainder = (N - 1) % size;

    int local_N = base + (rank < remainder ? 1 : 0);
    int start = 1 + rank * base + (rank < remainder ? rank : remainder);
    int end = start + local_N - 1;

    if (rank == 0) {
        initialize();
    }

    MPI_Bcast(u, N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double dx = L / N;
    double r = nu * dt / (2 * dx * dx);

    for (int nstep = 0; nstep < steps; nstep++) {
        double* d_local = (double*)malloc(local_N * sizeof(double));

        // Compute RHS (explicit part) on each local segment
        for (int i = start; i <= end; i++) {
            double convection = u[i] * (u[i + 1] - u[i - 1]) / (2 * dx);
            double diffusion = r * (u[i + 1] - 2 * u[i] + u[i - 1]);
            d_local[i - start] = u[i] - dt * convection + diffusion;
        }

        // Set up gather information
        int* counts = NULL;
        int* displs = NULL;
        if (rank == 0) {
            counts = (int*)malloc(size * sizeof(int));
            displs = (int*)malloc(size * sizeof(int));
            for (int i = 0; i < size; i++) {
                counts[i] = base + (i < remainder ? 1 : 0);
                displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
            }
        }

        // Gather full RHS vector `d` on rank 0
        MPI_Gatherv(d_local, local_N, MPI_DOUBLE, d, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        free(d_local);

        if (rank == 0) {
            build_matrix();
            solve_tridiagonal();

            for (int i = 0; i <= N; i++) {
                u[i] = u_new[i];
            }

            free(counts);
            free(displs);
        }

        // Broadcast updated solution to all processes
        MPI_Bcast(u, N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        write_output();
    }

    MPI_Finalize();
    return 0;
}

