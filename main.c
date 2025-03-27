#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Declare function from simpsons.c
double simpsons_rule(double a, double b, int n, int num_threads);

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <a> <b> <n>\n", argv[0]);
        return 1;
    }

    double a = atof(argv[1]);
    double b = atof(argv[2]);
    int n = atoi(argv[3]);

    // Exact integral of sin(x) from 0 to π is 2
    double exact_value = 2.0;
    int thread_counts[] = {1, 2, 4, 8, 16};
    FILE *file = fopen("output.txt", "w");

    if (file == NULL) {
        printf("Error opening file for writing!\n");
        return 1;
    }

    for (int i = 0; i < 5; i++) {
        int num_threads = thread_counts[i];
        double start_time = omp_get_wtime();
        double result = simpsons_rule(a, b, n, num_threads);
        double end_time = omp_get_wtime();

        double error = fabs(result - exact_value);

        // Ensure only one thread writes to file at a time
        #pragma omp critical
        {
            printf("Threads: %d, Integral: %.10f, Error: %.10e, Time: %f sec\n", num_threads, result, error, end_time - start_time);
            fprintf(file, "Threads: %d, Integral: %.10f, Error: %.10e, Time: %f sec\n", num_threads, result, error, end_time - start_time);
        }
    }

    fclose(file);
    return 0;
}

