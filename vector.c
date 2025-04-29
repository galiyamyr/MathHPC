#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

void usage(const char* prog_name) {
    fprintf(stderr, "Usage: %s <thread_count> <vector_size> <k>\n", prog_name);
    fprintf(stderr, "thread_count, vector_size, and k should be positive integers.\n");
    exit(1);
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        usage(argv[0]);
    }

    const int thread_count = strtol(argv[1], NULL, 10);
    const int N = strtol(argv[2], NULL, 10);
    const int k = strtol(argv[3], NULL, 10);

    if (thread_count <= 0 || N <= 0 || k <= 0) {
        usage(argv[0]);
    }

    srand(time(NULL));
    double* x = (double*)malloc(N * sizeof(double));
    double* u = (double*)malloc(N * sizeof(double));

    if (x == NULL || u == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Generate random values for x
    for (int i = 0; i < N; i++) {
        x[i] = (double)rand() / RAND_MAX;  // Random values between 0 and 1
    }

    double time1 = omp_get_wtime();

    #pragma omp parallel for num_threads(thread_count)
    for (int i = 0; i < N; i++) {
        double sum = 1.0;
        double power = x[i];
        for (int j = 1; j <= k; j++) {
            sum += power;
            power *= x[i];  // Efficiently compute x^j iteratively
        }
        u[i] = sum;
    }

    double time2 = omp_get_wtime();
    double clock_time = time2 - time1;

    printf("With %d threads, clock_time = %e sec\n", thread_count, clock_time);

    free(x);
    free(u);

    return 0;
}
