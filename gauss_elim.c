#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void gauss_elim(int N, double matrix[N][N], double b[N], double x[N]){
    for (int i=0; i<N-1; i++){
        for (int j=i+1; j<N; j++){
            double factor=matrix[j][i]/matrix[i][i];
            for (int k=0; k<N; k++){
                matrix[j][k]-=matrix[i][k]*factor;
            }
            b[j]-=b[i]*factor;
        }
    }
    for (int i = N - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < N; j++) {  // Fix: j should increment
            x[i] -= matrix[i][j] * x[j];
        }
        x[i] /= matrix[i][i]; // Fix: Final division only once
    }
}

