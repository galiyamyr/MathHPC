#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double vec_norm(int N, double vec[N]){
    double sum=0;
    for (int i=0; i<N; i++){
        sum+=vec[i]*vec[i];
    }
    return sqrt(sum);
}
void normalize(int N, double vec[N]){
    double norm=vec_norm(N, vec);
    for (int i=0; i<N; i++){
        vec[i]/=norm;
    }

}
void matrix_vec_mult(int N, double A[N][N], double vec[N], double result[N]){
    for (int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            result[i]+=A[i][j]*vec[j];

        }
    }

}
double power2(int N, double A[N][N], double vec[N], int maxiter, double tol){
    double lambda_old=0;
    double lambda_new=0;
    double w[N];
    int mstop=0;
    int k=0;
    normalize(N, vec);
    matrix_vec_mult(N, A, vec, w);
        for (int i=0; i<N; i++){
        lambda_old+=vec[i]*w[i];
    }
    while(mstop!=1){
        k+=1;
        matrix_vec_mult(N, A, vec, w);
        normalize(N, w);
        lambda_new=0;
        for(int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                lambda_new+=w[i]*A[i][j]*w[j];
            }
                    
        }
        if (fabs(lambda_new-lambda_old)<tol || k>=maxiter){
            mstop=1;
        }
        lambda_old=lambda_new;
        for(int i=0; i<N; i++){
            vec[i]=w[i];        }

    }
    return lambda_new;
    

}