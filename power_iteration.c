#include <math.h>
#include <stdlib.h>
#include "config.h"

void power_iteration(double **A, double *eigenvalue, double *eigenvector, int N) {
    double *b = malloc(N * sizeof(double));
    double *b_next = malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) b[i] = 1.0;

    for (int iter = 0; iter < MAX_ITER; iter++) {
        for (int i = 0; i < N; i++) {
            b_next[i] = 0;
            for (int j = 0; j < N; j++) {
                b_next[i] += A[i][j] * b[j];
            }
        }

        double norm = 0;
        for (int i = 0; i < N; i++) norm += b_next[i] * b_next[i];
        norm = sqrt(norm);

        for (int i = 0; i < N; i++) b[i] = b_next[i] / norm;

        *eigenvalue = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                *eigenvalue += b[i] * A[i][j] * b[j];
            }
        }
    }

    for (int i = 0; i < N; i++) eigenvector[i] = b[i];

    free(b);
    free(b_next);
}
