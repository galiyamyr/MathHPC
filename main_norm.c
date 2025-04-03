#include <stdio.h>
#include <stdlib.h>
#include "parallel_norm.h"

int main() {
    int N = 10000;  
    double *vec_2norm = (double*)malloc(N * sizeof(double));
    double *vec_maxnorm = (double*)malloc(N * sizeof(double));
    
    // Initialize vectors
    for (int i = 0; i < N; i++) {
        vec_2norm[i] = 1.0;
        vec_maxnorm[i] = i + 1;
    }
    
    // Compute norms
    double norm_2 = compute_2norm_fine(vec_2norm, N);
    double norm_max = compute_maxnorm_fine(vec_maxnorm, N);
    
    // Normalize vectors
    double *normalized_2norm = normalize_2norm(vec_2norm, N);
    double *normalized_maxnorm = normalize_maxnorm(vec_maxnorm, N);
    
    // Print results
    printf("2-norm: %f\n", norm_2);
    printf("First component of normalized (2-norm): %f\n", normalized_2norm[0]);
    printf("Last component of normalized (2-norm): %f\n", normalized_2norm[N-1]);
    
    printf("Max-norm: %f\n", norm_max);
    printf("First component of normalized (max-norm): %f\n", normalized_maxnorm[0]);
    printf("Last component of normalized (max-norm): %f\n", normalized_maxnorm[N-1]);
    
    // Free memory
    free(vec_2norm);
    free(vec_maxnorm);
    free(normalized_2norm);
    free(normalized_maxnorm);
    
    return 0;
}
