#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define L_VAL 1.0       // Length of domain
#define T 1.0           // Final time
#define N 100           // Number of spatial points (global)
#define TIME_STEP 0.001 // Time step
#define nu 0.01         // Viscosity
#define B (nu * TIME_STEP / (2 * (L_VAL / (N - 1)) * (L_VAL / (N - 1))))  // Coefficient for diffusion term

// Declare LU decomposition and solving functions
void lu_decomposition(double *M, double *Lmat, double *Umat, int n);
void solve(double *Lmat, double *Umat, double *b, double *x, int n);
void solve_nonlinear(double *Lmat, double *Umat, double *b, double *sol, int M, double *v, double dx, double dt);

int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double dx = L_VAL / (N - 1);      // Spatial step
    int M = N - 2;                    // Number of interior points (excluding boundaries)
    int steps = (int)(T / TIME_STEP); // Total time steps

    // Allocate memory for arrays
    double *M_left = malloc(M * M * sizeof(double));  // Matrix for Crank-Nicolson
    double *Lmat = malloc(M * M * sizeof(double));    // Lower triangular matrix for LU decomposition
    double *Umat = malloc(M * M * sizeof(double));    // Upper triangular matrix for LU decomposition
    double *b = malloc(M * sizeof(double));           // RHS vector
    double *sol = malloc(M * sizeof(double));         // Solution vector

    // Allocate solution array and initialize with initial condition
    double *v = malloc(N * sizeof(double));  // Local solution vector for each process
    for (int i = 0; i < N; i++) {
        double x = i * dx;
        v[i] = sin(M_PI * x);  // Initial condition: u(x,0) = sin(pi x)
    }

    // Construct the left-hand matrix (Crank-Nicolson)
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            M_left[i * M + j] = 0.0;
            if (i == j)
                M_left[i * M + j] = 1 + 2 * B;
            else if (abs(i - j) == 1)
                M_left[i * M + j] = -B;
        }
    }

    lu_decomposition(M_left, Lmat, Umat, M);  // Perform LU decomposition for the matrix

    // Time-stepping loop
    for (int k = 1; k <= steps; k++) {
        solve_nonlinear(Lmat, Umat, b, sol, M, v, dx, TIME_STEP);

        // Apply boundary conditions (no-flux boundary condition)
        v[0] = 0.0; 
        v[N - 1] = 0.0;

        // Collect the local solutions to form the full solution (you may need MPI_Gather here depending on your setup)
    }

    // Output the final solution (you can use MPI to collect results and output to file)
    if (rank == 0) {
        FILE *fp = fopen("output.txt", "w");
        for (int i = 0; i < N; i++) {
            fprintf(fp, "%f ", v[i]);
        }
        fclose(fp);
    }

    // Clean up allocated memory
    free(M_left);
    free(Lmat);
    free(Umat);
    free(b);
    free(sol);
    free(v);

    MPI_Finalize();
    return 0;
}

// LU Decomposition: A = LU
void lu_decomposition(double *M, double *Lmat, double *Umat, int n) {
    for (int i = 0; i < n; i++) {
        // Upper triangular
        for (int j = i; j < n; j++) {
            Umat[i * n + j] = M[i * n + j];
            for (int k = 0; k < i; k++) {
                Umat[i * n + j] -= Lmat[i * n + k] * Umat[k * n + j];
            }
        }

        // Lower triangular
        for (int j = i; j < n; j++) {
            if (i == j)
                Lmat[i * n + j] = 1.0;
            else {
                Lmat[j * n + i] = M[j * n + i];
                for (int k = 0; k < i; k++) {
                    Lmat[j * n + i] -= Lmat[j * n + k] * Umat[k * n + i];
                }
                Lmat[j * n + i] /= Umat[i * n + i];
            }
        }
    }
}

// Solve Ax = b using LU decomposition
void solve(double *Lmat, double *Umat, double *b, double *x, int n) {
    double *y = malloc(n * sizeof(double));

    // Forward substitution Ly = b
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= Lmat[i * n + j] * y[j];
        }
    }

    // Backward substitution Ux = y
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= Umat[i * n + j] * x[j];
        }
        x[i] /= Umat[i * n + i];
    }

    free(y);
}

// Solve nonlinear equation (Convection + Diffusion)
void solve_nonlinear(double *Lmat, double *Umat, double *b, double *sol, int M, double *v, double dx, double dt) {
    double *u_new = malloc(M * sizeof(double));
    double tolerance = 1e-6;
    int max_iter = 10;

    // Initial guess is the previous time step
    for (int i = 0; i < M; i++) {
        u_new[i] = v[i + 1];  // Assuming v[0] and v[N-1] are boundary conditions
    }

    // Start Newton's iterations
    for (int iter = 0; iter < max_iter; iter++) {
        // Construct the right-hand side (RHS) including convection term
        for (int i = 1; i < M - 1; i++) {
            double u_x = (v[i + 1] - v[i - 1]) / (2 * dx);  // Central difference for du/dx
            b[i - 1] = v[i] + B * (v[i + 1] - 2 * v[i] + v[i - 1]) - dt * u_new[i] * u_x;
        }

        // Solve the system using LU decomposition
        solve(Lmat, Umat, b, sol, M);

        // Update the solution
        for (int i = 0; i < M; i++) {
            u_new[i] = sol[i];
        }

        // Check for convergence (if the solution is not changing much)
        double max_diff = 0.0;
        for (int i = 0; i < M; i++) {
            max_diff = fmax(max_diff, fabs(u_new[i] - v[i + 1]));
        }
        if (max_diff < tolerance) {
            break;
        }
    }

    // Copy the new solution back to v
    for (int i = 1; i < M - 1; i++) {
        v[i + 1] = u_new[i];
    }

    free(u_new);
}

