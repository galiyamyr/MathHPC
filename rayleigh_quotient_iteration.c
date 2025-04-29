#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"

// Function to perform Gaussian elimination to solve Ax = b
void gaussian_elimination(double **A, double *b, double *x, int N) {
    double **augmented = malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        augmented[i] = malloc((N + 1) * sizeof(double));
        for (int j = 0; j < N; j++) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][N] = b[i];
    }

    // Forward elimination
    for (int i = 0; i < N; i++) {
        // Pivoting
        for (int k = i + 1; k < N; k++) {
            if (fabs(augmented[k][i]) > fabs(augmented[i][i])) {
                double *temp = augmented[i];
                augmented[i] = augmented[k];
                augmented[k] = temp;
            }
        }
        
        // Make the diagonal element 1 and eliminate below
        for (int k = i + 1; k < N; k++) {
            double factor = augmented[k][i] / augmented[i][i];
            for (int j = i; j <= N; j++) {
                augmented[k][j] -= factor * augmented[i][j];
            }
        }
    }

    // Back substitution
    for (int i = N - 1; i >= 0; i--) {
        x[i] = augmented[i][N];
        for (int j = i + 1; j < N; j++) {
            x[i] -= augmented[i][j] * x[j];
        }
        x[i] /= augmented[i][i];
    }

    for (int i = 0; i < N; i++) {
        free(augmented[i]);
    }
    free(augmented);
}

void rayleigh_quotient_iteration(double **A, double *eigenvalue, double *eigenvector, int N) {
    double *v = malloc(N * sizeof(double));
    double *w = malloc(N * sizeof(double));
    double norm, lambda_old, lambda_new;
    int k = 0, mstop = 0;

    // Initialize v(0) with ones
    for (int i = 0; i < N; i++) {
        v[i] = 1.0;
    }
    
    // Normalize v(0)
    norm = 0.0;
    for (int i = 0; i < N; i++) {
        norm += v[i] * v[i];
    }
    norm = sqrt(norm);
    for (int i = 0; i < N; i++) {
        v[i] /= norm;
    }
    
    // Initial eigenvalue estimate
    lambda_old = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            lambda_old += v[i] * A[i][j] * v[j];
        }
    }
    
    while (!mstop) {
        k++;
        
        // Construct (A - lambda_old * I)
        double **A_shifted = malloc(N * sizeof(double*));
        for (int i = 0; i < N; i++) {
            A_shifted[i] = malloc(N * sizeof(double));
            for (int j = 0; j < N; j++) {
                A_shifted[i][j] = A[i][j];
                if (i == j) {
                    A_shifted[i][j] -= lambda_old;
                }
            }
        }
        
        // Solve (A - lambda_old * I) * w = v using Gaussian elimination
        gaussian_elimination(A_shifted, v, w, N);
        
        // Normalize w to get v(k)
        norm = 0.0;
        for (int i = 0; i < N; i++) {
            norm += w[i] * w[i];
        }
        norm = sqrt(norm);
        for (int i = 0; i < N; i++) {
            v[i] = w[i] / norm;
        }
        
        // Compute new eigenvalue estimate
        lambda_new = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                lambda_new += v[i] * A[i][j] * v[j];
            }
        }
        
        // Free A_shifted
        for (int i = 0; i < N; i++) {
            free(A_shifted[i]);
        }
        free(A_shifted);
        
        // Check stopping condition
        if (fabs(lambda_new - lambda_old) < EPSILON || k == MAX_ITER) {
            mstop = 1;
        }
        
        lambda_old = lambda_new;
    }
    
    // Store computed eigenvalue
    *eigenvalue = lambda_new;
    
    // Store computed eigenvector
    for (int i = 0; i < N; i++) {
        eigenvector[i] = v[i];
    }
    
    free(v);
    free(w);
}
