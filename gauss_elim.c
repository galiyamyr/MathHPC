/**
 Implements Gaussian elimination for solving linear systems.
 */

 #include <stdio.h>
 #include <math.h>
 
 /**
  *
  * 
  *  n Size of the matrix.
  *  A Coefficient matrix.
  * b Right-hand side vector (modified in-place).
  */
 void gauss_elim(int n, double A[n][n], double b[n]) {
     for (int i = 0; i < n - 1; i++) {
         // Partial pivoting: Find the maximum element in the current column
         int max = i;
         for (int j = i + 1; j < n; j++) {
             if (fabs(A[j][i]) > fabs(A[max][i])) {
                 max = j;
             }
         }
         // Swap rows if maximum is not in the current row
         if (max != i) {
             for (int k = 0; k < n; k++) {
                 double temp = A[i][k];
                 A[i][k] = A[max][k];
                 A[max][k] = temp;
             }
             double temp = b[i];
             b[i] = b[max];
             b[max] = temp;
         }
 
         // Eliminate row below pivot
         for (int j = i + 1; j < n; j++) {
             double factor = A[j][i] / A[i][i];
             for (int k = i; k < n; k++) {
                 A[j][k] -= factor * A[i][k];
             }
             b[j] -= factor * b[i];
         }
     }
 
     // Back-substitution
     for (int i = n - 1; i >= 0; i--) {
         for (int j = i + 1; j < n; j++) {
             b[i] -= A[i][j] * b[j];
         }
         b[i] /= A[i][i];
     }
 }
 