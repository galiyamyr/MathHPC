#include <stdio.h>
#include "gauss_elim.h"
#include "least_squares.h"

void least_squares(int m, int n, double x_vals[], double y_vals[], double coeffs[]) {
    int poly_order = m + 1; // m is the polynomial order, so matrix size is (m+1)
    
    double V[n][poly_order];  // Vandermonde Matrix
    double Vt[poly_order][n]; // Transpose of Vandermonde Matrix

    for (int i = 0; i < n; i++) {
        V[i][0] = 1.0;
        for (int j = 1; j < poly_order; j++) {
            V[i][j] = V[i][j - 1] * x_vals[i];
        }
    }

    for (int i = 0; i < poly_order; i++) {
        for (int j = 0; j < n; j++) {
            Vt[i][j] = V[j][i];
        }
    }

    double VtV[poly_order][poly_order + 1];

    for (int i = 0; i < poly_order; i++) {
        for (int j = 0; j < poly_order; j++) {
            VtV[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                VtV[i][j] += Vt[i][k] * V[k][j];
            }
        }
        VtV[i][poly_order] = 0.0;
        for (int k = 0; k < n; k++) {
            VtV[i][poly_order] += Vt[i][k] * y_vals[k];
        }
    }

    gauss_elim(poly_order, VtV, coeffs);
}
