#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 5  // Size of the matrix
#define MAX_ITER 100
#define TOLERANCE 1e-6

// Function to print a matrix
void print_matrix(double mat[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }
}

// Function to generate a lower triangular matrix L with 1s on the diagonal
void generate_lower_triangular(double L[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            // Diagonal elements set to 1
            if (i == j) {
                L[i][j] = 1.0;
            } else {
                // Below diagonal elements set to random values between 0 and 1
                L[i][j] = (double)rand() / RAND_MAX;
            }
        }
    }
}

// Function to print matrix L for debugging
void print_L_matrix(double L[N][N]) {
    printf("\nMatrix L (Lower Triangular):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f ", L[i][j]);
        }
        printf("\n");
    }
}

// Function to compute A = L * L^T
void compute_A(double L[N][N], double A[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = 0.0;
            for (int k = 0; k < N; k++) {
                A[i][j] += L[i][k] * L[j][k];
            }
        }
    }
}

// Function to check if any value is NaN or Inf in the matrix
int check_valid_matrix(double mat[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (isnan(mat[i][j]) || isinf(mat[i][j])) {
                return 0; // Invalid value found
            }
        }
    }
    return 1; // All values are valid
}

// Function to compute the Hessenberg form of a matrix
void hessenberg_form(double A[N][N], double H[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            H[i][j] = A[i][j];
        }
    }

    // Hessenberg reduction (simplified)
    for (int k = 0; k < N - 2; k++) {
        for (int i = k + 1; i < N; i++) {
            double factor = H[i][k] / H[k][k];
            if (fabs(H[k][k]) < 1e-10) {
                printf("Warning: Small pivot element in Hessenberg form at (%d, %d)\n", k, k);
            }
            for (int j = k; j < N; j++) {
                H[i][j] -= factor * H[k][j];
            }
            for (int j = k; j < N; j++) {
                H[k][j] -= factor * H[i][j];
            }
        }
    }
}

// Function to perform QR decomposition and update matrix H (Simplified QR)
void qr_algorithm(double H[N][N]) {
    double R[N][N], Q[N][N], H_new[N][N];
    
    // QR Decomposition and matrix updates
    for (int iter = 0; iter < MAX_ITER; iter++) {
        // Compute QR decomposition (simplified for this example)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                R[i][j] = H[i][j];
                Q[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        // Multiply R and Q to get H_new (approximation of QR update)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                H_new[i][j] = 0.0;
                for (int k = 0; k < N; k++) {
                    H_new[i][j] += R[i][k] * Q[k][j];
                }
            }
        }

        // Update H matrix
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                H[i][j] = H_new[i][j];
            }
        }

        // Check for convergence (off-diagonal elements)
        int converged = 1;
        for (int i = 1; i < N; i++) {
            if (fabs(H[i][i-1]) > TOLERANCE) {
                converged = 0;
                break;
            }
        }
        if (converged) break;
    }
}

// Function to extract the eigenvalues from the diagonal of H
void extract_eigenvalues(double H[N][N], double eigenvalues[N]) {
    for (int i = 0; i < N; i++) {
        eigenvalues[i] = H[i][i];
    }
}

int main() {
    // Seed random number generator
    srand(time(NULL));

    double L[N][N], A[N][N], H[N][N];
    double eigenvalues[N];

    // Step 1: Generate lower triangular matrix L with 1s on the diagonal
    generate_lower_triangular(L);
    
    // Step 2: Print the L matrix for debugging
    print_L_matrix(L);

    // Step 3: Compute A = L * L^T
    compute_A(L, A);

    // Step 4: Check if matrix A is valid (no NaN or Inf)
    if (!check_valid_matrix(A)) {
        printf("Error: A contains NaN or Inf\n");
        exit(1); // Exit if invalid value is found
    }

    // Step 5: Print matrix A
    printf("Original Matrix A (LL^T):\n");
    print_matrix(A);

    // Step 6: Bring A to Hessenberg form
    hessenberg_form(A, H);
    printf("\nHessenberg Form of A:\n");
    print_matrix(H);

    // Step 7: Apply QR algorithm to diagonalize H
    qr_algorithm(H);
    printf("\nDiagonalized Matrix (Approximately Diagonal):\n");
    print_matrix(H);

    // Step 8: Extract eigenvalues from diagonal
    extract_eigenvalues(H, eigenvalues);
    printf("\nEigenvalues (Diagonal elements):\n");
    for (int i = 0; i < N; i++) {
        printf("%f ", eigenvalues[i]);
    }
    printf("\n");

    return 0;
}
