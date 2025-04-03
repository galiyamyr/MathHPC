#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double vector_norm(double* vector, int N){
    double sum=0;
    for (int i=0; i<N; i++){
        sum+=vector[i]*vector[i];
    }
    return sqrt(sum);
}
void normalize(double* vector, int N){
    double norm=vector_norm(vector, N);
    for (int i=0; i<N; i++){
        vector[i]=vector[i]/norm;
    }
}
void matrix_vector_mult(int n, double matrix[][n], double *vec, double *result) {
    for (int i = 0; i < n; i++) {
        result[i] = 0;
        for (int j = 0; j < n; j++) {
            result[i] += matrix[i][j] * vec[j];
        }
    }
}
double power(int N, double A[N][N], double vec[N], int maxiter, double tol){
    double w[N];
    int k=0;
    double lambda_old=0;
    double lambda_new=0;
    int mstop=0;
    normalize(vec, N);
    matrix_vector_mult(N, A, vec, w);
    for (int i=0; i<N; i++){
        lambda_old+=vec[i]*w[i];
    }
while(mstop!=1){
    k++;
    matrix_vector_mult(N, A, vec, w);
    normalize(w,N);
    lambda_new=0;
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            lambda_new+=w[i]*A[i][j]*w[j];
        }
    }
    if (fabs(lambda_new-lambda_old)<tol || k>=maxiter){
        mstop=1;
    }
    lambda_old=lambda_new;
    for(int i=0; i<N; i++){
        vec[i]=w[i];
    }
}
return lambda_new;

}
int main() {
    int N = 3;
    double A[3][3] = {
        {4, 1, 0},
        {1, 3, 1},
        {0, 1, 2}
    };
    double vec[3] = {1, 1, 1};  // Initial guess

    int maxiter = 100;
    double tol = 1e-6;

    double dominant_eigenvalue = power(N, A, vec, maxiter, tol);

    // Print results
    printf("Dominant eigenvalue: %lf\n", dominant_eigenvalue);
    printf("Dominant eigenvector: [");
    for (int i = 0; i < N; i++) {
        printf("%lf ", vec[i]);
    }
    printf("]\n");

    return 0;
}