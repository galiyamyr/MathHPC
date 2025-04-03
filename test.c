#include<stdio.h>
unsigned long long factorial(int n){
if (n<0) return 0;
if (n==0||n==1) return 1;
return n*factorial(n-1);
}
double exponential(double x){
    double sum=1;
    double term=1;
    for (int i=1; i<=20; i++){
        term*=x/i;
        sum+=term;
    }
    return sum;
}
int main(){
    int n;
    double x;
    printf("Enter value for n:");
    scanf("%d",&n);
    printf("Factorial of %d=%llu\n", n, factorial(n));
    printf("Enter the value of x:");
    scanf("%lf", &x);
    printf("the exponent in power %lf is %lf", x, exponential(x));
}