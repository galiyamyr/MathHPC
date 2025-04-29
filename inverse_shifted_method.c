#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
// Function to perform Gaussian elimination to solve Ax = b
void gauss(double **A, double *b, double *x, int N) {
    double **U = malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        U[i] = malloc(N * sizeof(double));
        for (int j = 0; j < N; j++) {
            U[i][j] = A[i][j];
        }
    }
    double *y = malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        y[i] = b[i];
    }

    for (int i = 0; i < N; i++) {
        int maxRow = i;
        for (int k = i + 1; k < N; k++) {
            if (fabs(U[k][i]) > fabs(U[maxRow][i])) {
                maxRow = k;
            }
        }
        for (int k = 0; k < N; k++) {
            double temp = U[i][k];
            U[i][k] = U[maxRow][k];
            U[maxRow][k] = temp;
        }
        double temp = y[i];
        y[i] = y[maxRow];
        y[maxRow] = temp;

        if (fabs(U[i][i]) < EPSILON) { // Check for singular matrix
            printf("Matrix is singular or close to singular.\n");
            return;
        }

        for (int k = i + 1; k < N; k++) {
            double factor = U[k][i] / U[i][i];
            for (int j = i; j < N; j++) {
                U[k][j] -= factor * U[i][j];
            }
            y[k] -= factor * y[i];
        }
    }

    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    for (int i = 0; i < N; i++) {
        free(U[i]);
    }
    free(U);
    free(y);
}

// Inverse shifted method to compute eigenvalue and eigenvector
void inverse_shifted_method(double **A, double shift, double *eigenvalue, double *eigenvector, int N) {
    double *v = malloc(N * sizeof(double));
    double *w = malloc(N * sizeof(double));
    double norm, lambda_old, lambda_new;
    int k = 0, mstop = 0;

    // Initialize eigenvector v to 1
    for (int i = 0; i < N; i++) {
        v[i] = 1.0;
    }
    
    // Normalize v
    norm = 0.0;
    for (int i = 0; i < N; i++) {
        norm += v[i] * v[i];
    }
    norm = sqrt(norm);
    for (int i = 0; i < N; i++) {
        v[i] /= norm;
    }
    
    // Initial guess for eigenvalue (Rayleigh quotient)
    lambda_old = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            lambda_old += v[i] * A[i][j] * v[j];
        }
    }
    
    while (!mstop) {
        k++;
        
        // Create A_shifted matrix (A - shift * I)
        double **A_shifted = malloc(N * sizeof(double*));
        for (int i = 0; i < N; i++) {
            A_shifted[i] = malloc(N * sizeof(double));
            for (int j = 0; j < N; j++) {
                A_shifted[i][j] = A[i][j];
                if (i == j) {
                    A_shifted[i][j] -= shift;
                }
            }
        }
        
        // Solve (A - shift * I) v = w using Gaussian elimination
        gauss(A_shifted, v, w, N);
        
        // Normalize w to get the new eigenvector v
        norm = 0.0;
        for (int i = 0; i < N; i++) {
            norm += w[i] * w[i];
        }
        norm = sqrt(norm);
        for (int i = 0; i < N; i++) {
            v[i] = w[i] / norm;
        }
        
        // Compute new eigenvalue using Rayleigh quotient
        lambda_new = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                lambda_new += v[i] * A[i][j] * v[j];
            }
        }
        
        // Free the shifted matrix
        for (int i = 0; i < N; i++) {
            free(A_shifted[i]);
        }
        free(A_shifted);
        
        // Check for convergence
        if (fabs(lambda_new - lambda_old) < EPSILON || k == MAX_ITER) {
            mstop = 1;
        }
        
        lambda_old = lambda_new;
    }
    
    // Return the computed eigenvalue and eigenvector
    *eigenvalue = lambda_new;
    for (int i = 0; i < N; i++) {
        eigenvector[i] = v[i];
    }
    
    free(v);
    free(w);
}
