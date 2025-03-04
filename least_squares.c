/**
 *  least_squares.c
 *  Computes polynomial coefficients using the least squares method.
 *
 * Constructs the Vandermonde matrix and solves the normal equations
 * (A^T A x = A^T b) to obtain polynomial coefficients that best fit the given data.
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #include "least_squares.h"
 #include "gauss_elim.h"
 
 /**
  * Computes least squares polynomial fitting.
  * 
  * m Degree of the polynomial.
  * n Number of data points.
  * x Array of x values.
  *  y Array of y values.
  * coeffs Output array for polynomial coefficients.
  */
 void least_squares(int m, int n, double x[], double y[], double coeffs[]) {
     double A[n][m + 1];    // Rectangular Vandermonde matrix (n × (m+1))
     double ATA[m + 1][m + 1];  // Square matrix A^T A ((m+1) × (m+1))
     double ATb[m + 1];     // Right-hand side vector A^T b ((m+1) × 1)
 
     // Construct the Vandermonde matrix A (n × (m+1))
     for (int i = 0; i < n; i++) {
         for (int j = 0; j <= m; j++) {
             A[i][j] = pow(x[i], j);
         }
     }
 
     // Compute A^T A ((m+1) × (m+1))
     for (int i = 0; i <= m; i++) {
         for (int j = 0; j <= m; j++) {
             ATA[i][j] = 0;
             for (int k = 0; k < n; k++) {
                 ATA[i][j] += A[k][i] * A[k][j];  // Sum of (A^T * A)
             }
         }
     }
 
     // Compute A^T b ((m+1) × 1)
     for (int i = 0; i <= m; i++) {
         ATb[i] = 0;
         for (int k = 0; k < n; k++) {
             ATb[i] += A[k][i] * y[k];  // Sum of (A^T * b)
         }
     }
 
     // Solve the system (A^T A) x = (A^T b) using Gaussian elimination
     gauss_elim(m + 1, ATA, ATb);
 
     // Store the computed coefficients
     for (int i = 0; i <= m; i++) {
         coeffs[i] = ATb[i];
     }
 }
 