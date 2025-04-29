#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#define NU 0.01  // viscosity

// Initial condition
double u0(double x) {
    return sin(M_PI * x);
}

// Right-hand side of Burgers' equation (could be zero)
double f(double x, double t) {
    return 0.0;
}

int main(int argc, char *argv[]) {
    int N = 1000;          // default number of grid points
    int nsteps = 1000;     // default number of time steps
    double T = 1.0;

    if (argc >= 2) N = atoi(argv[1]);
    if (argc >= 3) nsteps = atoi(argv[2]);

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double dx = 1.0 / (N - 1);
    double dt = T / nsteps;
    int local_N = N / size;
    if (rank == size - 1) local_N += N % size; // handle remainder

    double *u = malloc((local_N + 2) * sizeof(double)); // +2 for ghost cells
    double *unew = malloc((local_N + 2) * sizeof(double));

    // Initialize
    for (int i = 1; i <= local_N; i++) {
        double x = (rank * (N / size) + i - 1) * dx;
        u[i] = u0(x);
    }

    double t_start = MPI_Wtime();

    for (int step = 0; step < nsteps; step++) {
        // Communicate ghost cells
        if (rank > 0)
            MPI_Sendrecv(&u[1], 1, MPI_DOUBLE, rank - 1, 0,
                         &u[0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (rank < size - 1)
            MPI_Sendrecv(&u[local_N], 1, MPI_DOUBLE, rank + 1, 0,
                         &u[local_N + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 1; i <= local_N; i++) {
            double ux = (u[i + 1] - u[i - 1]) / (2 * dx);
            double uxx = (u[i + 1] - 2 * u[i] + u[i - 1]) / (dx * dx);
            unew[i] = u[i] + dt * (-u[i] * ux + NU * uxx);
        }

        double *temp = u;
        u = unew;
        unew = temp;
    }

    double t_end = MPI_Wtime();

    if (rank == 0) {
        printf("Execution time: %.3f sec\n", t_end - t_start);
    }

    free(u);
    free(unew);
    MPI_Finalize();
    return 0;
}

