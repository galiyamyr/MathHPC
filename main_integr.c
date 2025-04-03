#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

void usage(const char* prog_name);
double AdaptiveInt(const double a, const double b, const double TOL, char* filename, const int num_threads);
double AdaptiveIntSerial(const double a, const double b, const double TOL, char* filename);
double Q(const double a, const double b);
double f(const double x);

int main(int argc, char* argv[]) {
    if (argc != 3) {
        usage(argv[0]);
    }

    int thread_count = strtol(argv[1], NULL, 10);
    double TOL = strtod(argv[2], NULL);
    if (thread_count <= 0 || TOL < 5.0e-16) {
        usage(argv[0]);
    }

    double a = -2.0;
    double b = 4.0;
    char filename[] = "quadrature.data";
    FILE* file = fopen(filename, "w");
    fclose(file);
    
    const double Iex = 0.4147421694070212;
    omp_set_nested(1);
    omp_set_dynamic(0);

    double start_time = omp_get_wtime();
    double I;
    if (thread_count == 1) {
        I = AdaptiveIntSerial(a, b, TOL, filename);
    } else {
        I = AdaptiveInt(a, b, TOL, filename, thread_count);
    }
    double end_time = omp_get_wtime();

    printf("\nThread Count = %d\n", thread_count);
    printf("TOL = %24.15e\n", TOL);
    printf("Integral Approximation = %24.15e\n", I);
    printf("Error = %24.15e\n", fabs(I - Iex));
    printf("Execution Time = %f seconds\n", end_time - start_time);
    
    return 0;
}

void usage(const char* prog_name) {
    fprintf(stderr, "Usage: %s <num_threads> <TOL>\n", prog_name);
    fprintf(stderr, "num_threads should be positive\n");
    fprintf(stderr, "TOL should be positive\n");
    exit(1);
}

double AdaptiveIntSerial(const double a, const double b, const double TOL, char* filename) {
    double Qab = Q(a, b);
    double c = 0.5 * (a + b);
    double Qac = Q(a, c);
    double Qcb = Q(c, b);
    double error_est = (1.0 / 15.0) * fabs(Qac + Qcb - Qab);
    double result;
    
    int my_rank = omp_get_thread_num();
    FILE* file = fopen(filename, "a");
    fprintf(file, "%3i %24.15e %24.15e\n", my_rank + 1, a, b);
    fclose(file);
    
    if (error_est < TOL) {
        result = Qac + Qcb;
    } else {
        result = AdaptiveIntSerial(a, c, 0.5 * TOL, filename);
        result += AdaptiveIntSerial(c, b, 0.5 * TOL, filename);
    }
    return result;
}

double AdaptiveInt(const double a, const double b, const double TOL, char* filename, const int num_threads) {
    double Qab = Q(a, b);
    double c = 0.5 * (a + b);
    double Qac = Q(a, c);
    double Qcb = Q(c, b);
    double error_est = (1.0 / 15.0) * fabs(Qac + Qcb - Qab);
    double result = 0.0;

    FILE* file = fopen(filename, "a");
    fprintf(file, "%3i %24.15e %24.15e\n", omp_get_thread_num() + 1, a, b);
    fclose(file);
    
    if (error_est < TOL) {
        return Qac + Qcb;
    }

    #pragma omp parallel num_threads(num_threads)
    {
        double local_result = 0.0;
        int my_rank = omp_get_thread_num();
        double chunk_size = (b - a) / num_threads;
        double a_local = a + my_rank * chunk_size;
        double b_local = a_local + chunk_size;

        local_result = AdaptiveInt(a_local, b_local, TOL / num_threads, filename, num_threads);
        
        #pragma omp atomic
        result += local_result;
    }
    return result;
}

double Q(const double a, const double b) {
    double c = 0.5 * (a + b);
    return (1.0 / 6.0) * (b - a) * (f(a) + 4.0 * f(c) + f(b));
}

double f(const double x) {
    const double beta = 10.0;
    return exp(-pow(beta * x, 2)) + sin(x);
}
