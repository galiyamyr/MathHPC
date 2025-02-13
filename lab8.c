#include<stdio.h>
#include<math.h>
unsigned long long factorial(int n){
    if (n<0) return 0;
    if (n==0||n==1) return 1;
    return n*factorial(n-1);
    }

double my_exp(double x){
    int terms=10;
    double sum=1.0;
    int x0=round(x);
    double z=x-x0;
    double exp_x0=pow(2.718281828459, x0);
    for (int i=1; i<terms; i++){
        sum+=pow(z,i)/factorial(i);
    }
    return exp_x0*sum;
}
int main(){
    FILE *file=fopen("exp_data.txt", "w");
    if (file==NULL){
        printf("Error opening file\n");
        return 1;
    }
    for (double x=0.0; x<=1.0; x+=0.02){
        double exp_val=my_exp(x);
        printf("exp(%lf)=%lf\n", x, exp_val);
        fprintf(file, "%lf=%lf\n", x, exp_val);
    }
    fclose(file);
    printf("Data saved to exp_data.txt\n");
    return 0;
}