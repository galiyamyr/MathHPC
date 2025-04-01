#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Fine-grain parallelization for 2-norm computation
double compute_2norm_fine(double *vec, int N) {
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; i++) {
        sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}

// Coarse-grain parallelization for 2-norm computation
double compute_2norm_coarse(double *vec, int N) {
    double sum = 0.0;
    int num_threads = omp_get_max_threads();
    double *partial_sums = (double*)calloc(num_threads, sizeof(double));

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        int chunk_size = (N + num_threads - 1) / num_threads;
        int start = id * chunk_size;
        int end = (start + chunk_size > N) ? N : start + chunk_size;
        
        for (int i = start; i < end; i++) {
            partial_sums[id] += vec[i] * vec[i];
        }
    }

    for (int i = 0; i < num_threads; i++) {
        sum += partial_sums[i];
    }
    free(partial_sums);
    return sqrt(sum);
}

// Fine-grain parallelization for max-norm computation
double compute_maxnorm_fine(double *vec, int N) {
    double max_val = 0.0;
    #pragma omp parallel for reduction(max:max_val)
    for (int i = 0; i < N; i++) {
        if (fabs(vec[i]) > max_val) {
            max_val = fabs(vec[i]);
        }
    }
    return max_val;
}

// Coarse-grain parallelization for max-norm computation
double compute_maxnorm_coarse(double *vec, int N) {
    int num_threads = omp_get_max_threads();
    double *partial_max = (double*)calloc(num_threads, sizeof(double));

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        int chunk_size = (N + num_threads - 1) / num_threads;
        int start = id * chunk_size;
        int end = (start + chunk_size > N) ? N : start + chunk_size;
        
        for (int i = start; i < end; i++) {
            if (fabs(vec[i]) > partial_max[id]) {
                partial_max[id] = fabs(vec[i]);
            }
        }
    }

    double max_val = 0.0;
    for (int i = 0; i < num_threads; i++) {
        if (partial_max[i] > max_val) {
            max_val = partial_max[i];
        }
    }
    free(partial_max);
    return max_val;
}

// Normalize vector using 2-norm
double* normalize_2norm(double *vec, int N) {
    double norm = compute_2norm_fine(vec, N);
    double *normalized_vec = (double*)malloc(N * sizeof(double));
    
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        normalized_vec[i] = vec[i] / norm;
    }
    return normalized_vec;
}

// Normalize vector using max-norm
double* normalize_maxnorm(double *vec, int N) {
    double norm = compute_maxnorm_fine(vec, N);
    double *normalized_vec = (double*)malloc(N * sizeof(double));
    
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        normalized_vec[i] = vec[i] / norm;
    }
    return normalized_vec;
}

