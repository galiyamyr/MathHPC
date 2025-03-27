#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Function to integrate
double f(double x) {
    return sin(x);  // Example function: sin(x)
}

// Simpson’s Rule with OpenMP
double simpsons_rule(double a, double b, int n, int num_threads) {
    if (n % 2 != 0) n++;  // Ensure n is even for Simpson's Rule

    double h = (b - a) / n;
    double integral = f(a) + f(b);

    // Parallel computation using reduction for sum
    #pragma omp parallel num_threads(num_threads)
    {
        double local_sum = 0.0;
        
        #pragma omp for 
        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            if (i % 2 == 0)
                local_sum += 2 * f(x);
            else
                local_sum += 4 * f(x);
        }

        // Reduction to ensure correct summation
        #pragma omp critical
        {
            integral += local_sum;
        }
    }

    integral *= h / 3.0;
    return integral;
}

