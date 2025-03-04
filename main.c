/**
 * @file main.c
 * @brief Computes least squares and interpolation polynomials for sin(2πx) and visualizes results.
 *
 * This program generates N data points (x, y) where y = sin(2πx), fits a 
 * least squares polynomial (degree M) and an interpolation polynomial (degree N-1),
 * and then calls a Python script to visualize the results.
 */

 #include <stdio.h>
 #include <math.h>
 #include <stdlib.h>
 #include "least_squares.h"
 
 #define N 20  // Number of data points
 #define M 5   // Polynomial order (Least Squares)
 #define INTERP_ORDER (N-1)  // Polynomial interpolation order
 
 /**
  * @brief Writes polynomial coefficients to a file.
  * 
  * @param filename Name of the file to write to.
  * @param coeffs Array of polynomial coefficients.
  * @param degree Degree of the polynomial.
  */
 void write_coeffs_to_file(const char *filename, double coeffs[], int degree) {
     FILE *file = fopen(filename, "w");
     if (file == NULL) {
         printf("Error opening file!\n");
         exit(1);
     }
     for (int i = 0; i <= degree; i++) {
         fprintf(file, "%lf\n", coeffs[i]);
     }
     fclose(file);
 }
 
 int main() {
     double x_vals[N], y_vals[N];
     double coeffs_m[M + 1], coeffs_n[INTERP_ORDER + 1];
 
     // Generate x values from 0 to 1 and compute y = sin(2πx)
     for (int i = 0; i < N; i++) {
         x_vals[i] = (double)i / (N - 1);
         y_vals[i] = sin(2 * M_PI * x_vals[i]);
     }
 
     // Compute Least Squares Polynomial (Degree M)
     least_squares(M, N, x_vals, y_vals, coeffs_m);
     write_coeffs_to_file("coeffs_m.txt", coeffs_m, M);
 
     // Compute Interpolating Polynomial (Degree N-1)
     least_squares(INTERP_ORDER, N, x_vals, y_vals, coeffs_n);
     write_coeffs_to_file("coeffs_n.txt", coeffs_n, INTERP_ORDER);
 
     // Call Python script to visualize results
     system("python3 plot_polynomials.py");
 
     return 0;
 }
 