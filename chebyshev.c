#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    const int Nmax = 5;  // Maximum order
    const int NumPts = 21;  // Number of points
    int N;

    // Input polynomial order
    printf("\n Input order (0 <= N <= %d): ", Nmax);
    scanf("%d", &N);
    if (N < 0 || N > Nmax) {
        printf("Error: N = %d out of range.\n", N);
        printf("Order has to be: 0 <= N <= %d\n\n", Nmax);
        exit(1);
    }

    printf("\n");

    // Input polynomial coefficients
    double b[Nmax + 1];
    for (int i = 0; i <= N; i++) {
        printf("Set coefficient b[%d]: ", i);
        scanf("%lf", &b[i]);
    }

    printf("\n");

    // Generate x values
    double x[NumPts];
    for (int i = 0; i < NumPts; i++) {
        x[i] = -1.0 + i * (2.0 / (NumPts - 1));
    }

    // Compute Chebyshev polynomial expansion
    double y[NumPts];
    void ComputeChebyshev(const int N, const int NumPts, const double b[], const double x[], double y[]);
    ComputeChebyshev(N, NumPts, b, x, y);

    void WritePoly(const int NumPts, const double x[], double y[]);
    WritePoly(NumPts, x, y);
    system("python PlotPoly.py");
    return 0;
}

void ComputeChebyshev(const int N, const int NumPts, const double b[], const double x[], double y[]) {
    for (int i = 0; i < NumPts; i++) {
        double T[N + 1];  // 
        T[0] = 1;
        if (N > 0) {
            T[1] = x[i];
        }
        for (int j = 2; j <= N; j++) {  
            T[j] = 2 * x[i] * T[j - 1] - T[j - 2];
        }

        y[i] = 0.0;
        for (int j = 0; j <= N; j++) {  
            y[i] += b[j] * T[j];
        }
    }
}


void WritePoly(const int NumPts, const double x[], double y[]) {
    FILE *file = fopen("poly.data", "w");
    if (file == NULL) {
        printf("Error: cannot open the file.\n");
        exit(1);
    }

    for (int i = 0; i < NumPts; i++) {
        fprintf(file, "%lf %lf\n", x[i], y[i]);
    }

    fclose(file);
}


   
   
