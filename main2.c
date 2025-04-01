#include<stdio.h>
int main(){
    double A[3][3]={
        {1,2,3},{4,5,6},{7,8,9}
    };
    double b[3]={1,2,3};
    double x[3];
    printf("The solution is:");
    for (int i=0; i<3; i++){
        printf("x[%d]=%lf", i, x[i]);
    }
    return 0;

}