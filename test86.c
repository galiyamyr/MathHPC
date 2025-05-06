#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define N 10000
#define L 1.0
#define T 1.0
#define nu 0.01
#define dt 0.001
#define PI 3.141592653589793

void initialize(double* u, int local_N, int start) {
    for (int i = 0; i < local_N; i++) {
        double x = (start + i) * L / N;
        u[i + 1] = sin(PI * x);  // +1 to skip ghost cell
    }
}

void exchange_ghost_cells(double* u, int local_N, int rank, int size, MPI_Comm comm) {
    MPI_Status status;
    if (rank > 0) {
        MPI_Send(&u[1], 1, MPI_DOUBLE, rank - 1, 0, comm);
        MPI_Recv(&u[0], 1, MPI_DOUBLE, rank - 1, 0, comm, &status);
    } else {
        u[0] = 0.0; // Dirichlet BC
    }
    if (rank < size - 1) {
        MPI_Send(&u[local_N], 1, MPI_DOUBLE, rank + 1, 0, comm);
        MPI_Recv(&u[local_N + 1], 1, MPI_DOUBLE, rank + 1, 0, comm, &status);
    } else {
        u[local_N + 1] = 0.0; // Dirichlet BC
    }
}

void apply_A(double* u_full, double* x, double* Ax, int local_N, double r, int rank, int size, MPI_Comm comm) {
    double* x_full = calloc(local_N + 2, sizeof(double));
    for (int i = 0; i < local_N; i++) x_full[i + 1] = x[i];
    exchange_ghost_cells(x_full, local_N, rank, size, comm);

    for (int i = 0; i < local_N; i++) {
        Ax[i] = (1 + 2 * r) * x_full[i + 1] - r * x_full[i] - r * x_full[i + 2];
    }

    free(x_full);
}

double dot(double* x, double* y, int n, MPI_Comm comm) {
    double local = 0.0, global = 0.0;
    for (int i = 0; i < n; i++) local += x[i] * y[i];
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, comm);
    return global;
}

void axpy(double* y, double alpha, double* x, int n) {
    for (int i = 0; i < n; i++) y[i] += alpha * x[i];
}

void axpby(double* y, double alpha, double* x, double beta, int n) {
    for (int i = 0; i < n; i++) y[i] = alpha * x[i] + beta * y[i];
}

void cg_solve(double* u, double* d, int local_N, double r, int rank, int size, MPI_Comm comm) {
    double *rvec = malloc(local_N * sizeof(double));
    double *p = malloc(local_N * sizeof(double));
    double *Ap = malloc(local_N * sizeof(double));

    apply_A(u, u + 1, Ap, local_N, r, rank, size, comm);

    for (int i = 0; i < local_N; i++) {
        rvec[i] = d[i] - Ap[i];
        p[i] = rvec[i];
    }

    double rdot = dot(rvec, rvec, local_N, comm);

    for (int iter = 0; iter < 1000; iter++) {
        apply_A(u, p, Ap, local_N, r, rank, size, comm);

        double alpha = rdot / dot(p, Ap, local_N, comm);
        axpy(u + 1, alpha, p, local_N);
        axpy(rvec, -alpha, Ap, local_N);

        double rdot_new = dot(rvec, rvec, local_N, comm);
        if (sqrt(rdot_new) < 1e-6) break;

        double beta = rdot_new / rdot;
        axpby(p, 1.0, rvec, beta, local_N);
        rdot = rdot_new;
    }

    free(rvec); free(p); free(Ap);
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

    double *u = calloc(local_N + 2, sizeof(double)); // includes ghost cells
    double *d = malloc(local_N * sizeof(double));

    initialize(u, local_N, start);
    double dx = L / N;
    double r = nu * dt / (dx * dx);

    for (int nstep = 0; nstep < steps; nstep++) {
        exchange_ghost_cells(u, local_N, rank, size, MPI_COMM_WORLD);
        for (int i = 0; i < local_N; i++) {
            double convection = u[i + 1] * (u[i + 2] - u[i]) / (2 * dx);
            double diffusion = r * (u[i + 2] - 2 * u[i + 1] + u[i]);
            d[i] = u[i + 1] - dt * convection + diffusion;
        }
        cg_solve(u, d, local_N, r, rank, size, MPI_COMM_WORLD);
    }

    double* full_u = NULL;
    int* counts = NULL;
    int* displs = NULL;
    if (rank == 0) {
        full_u = malloc((N + 1) * sizeof(double));
        counts = calloc(size, sizeof(int));
        displs = calloc(size, sizeof(int));

        for (int i = 0; i < size; i++) {
            counts[i] = base + (i < remainder ? 1 : 0);
            displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
        }
    }

    MPI_Gatherv(u + 1, local_N, MPI_DOUBLE, full_u + 1, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        full_u[0] = 0.0;
        full_u[N] = 0.0;
        FILE* fp = fopen("burgers_output.txt", "w");
        for (int i = 0; i <= N; i++) {
            fprintf(fp, "%f %f\n", i * dx, full_u[i]);
        }
        fclose(fp);
        free(full_u);
        free(counts);
        free(displs);
    }

    free(u);
    free(d);
    MPI_Finalize();
    return 0;
}

