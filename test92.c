#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define nu 0.01
#define T 1.0

// Function for initial condition
double u0(double x) {
    return sin(M_PI * x);
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) fprintf(stderr, "Usage: %s N\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    int N = atoi(argv[1]); // Total number of points
    if (N % size != 0) {
        if (rank == 0) fprintf(stderr, "N must be divisible by number of processes!\n");
        MPI_Finalize();
        return 1;
    }

    double L = 1.0;
    double dx = L / (N - 1);
    double dt = 0.4 * dx * dx / nu;
    int Nt = T / dt;

    int local_N = N / size;
    int start = rank * local_N;
    int end = start + local_N;

    double *u = malloc((local_N + 2) * sizeof(double));
    double *unew = malloc((local_N + 2) * sizeof(double));

    for (int i = 0; i < local_N + 2; i++) {
        int global_idx = start + i - 1;
        double x = global_idx * dx;
        if (global_idx >= 0 && global_idx < N)
            u[i] = u0(x);
        else
            u[i] = 0.0;
    }

    for (int t = 0; t < Nt; t++) {
        MPI_Sendrecv(&u[1], 1, MPI_DOUBLE, (rank + size - 1) % size, 0,
                     &u[0], 1, MPI_DOUBLE, (rank + 1) % size, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&u[local_N], 1, MPI_DOUBLE, (rank + 1) % size, 1,
                     &u[local_N + 1], 1, MPI_DOUBLE, (rank + size - 1) % size, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 1; i <= local_N; i++) {
            unew[i] = u[i] - dt / (2 * dx) * u[i] * (u[i + 1] - u[i - 1])
                          + nu * dt / (dx * dx) * (u[i + 1] - 2 * u[i] + u[i - 1]);
        }

        double *temp = u;
        u = unew;
        unew = temp;
    }

    // Output final result from rank 0
    if (rank == 0) {
        FILE *f = fopen("output_final.txt", "w");
        for (int i = 0; i < local_N; i++) {
            double x = (i + start) * dx;
            fprintf(f, "%f %f\n", x, u[i + 1]);
        }

        for (int p = 1; p < size; p++) {
            for (int i = 0; i < local_N; i++) {
                double val;
                MPI_Recv(&val, 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                double x = (i + p * local_N) * dx;
                fprintf(f, "%f %f\n", x, val);
            }
        }
        fclose(f);
    } else {
        for (int i = 0; i < local_N; i++) {
            MPI_Send(&u[i + 1], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    free(u);
    free(unew);
    MPI_Finalize();
    return 0;
}

