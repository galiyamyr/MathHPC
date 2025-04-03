#include <stdio.h>

// Function to compute factorial (recursive)
unsigned long long factorial(int n) {
    if (n < 0) return 0; 
    if (n == 0 || n == 1) return 1;
    return n * factorial(n - 1);
}

// Function to compute exponential using Taylor series
double exponential(double x) {
    double sum = 1.0;  // First term of the series (1)
    double term = 1.0;
    int n = 1;
    
    // Loop to calculate the series sum
    for (int i = 1; i <= 20; i++) {  // 20 terms for better accuracy
        term *= x / i;  // x^i / i!
        sum += term;
    }
    
    return sum;
}

// Function to compute natural logarithm using Newton-Raphson method
double logarithm(double y) {
    if (y <= 0) {
        printf("Logarithm is undefined for non-positive numbers.\n");
        return 0;
    }

    double x = y - 1;  // Initial guess
    for (int i = 0; i < 20; i++) {  // Iterate 20 times for accuracy
        x = x - (exponential(x) - y) / exponential(x);
    }

    return x;
}

int main() {
    int n;
    double x, y;
    
    // Factorial input
    printf("Enter an integer for factorial calculation: ");
    scanf("%d", &n);
    printf("Factorial of %d = %llu\n", n, factorial(n));

    // Exponential input
    printf("Enter a real number to compute exponential: ");
    scanf("%lf", &x);
    printf("exp(%.2lf) ≈ %.5lf\n", x, exponential(x));

    // Logarithm input
    printf("Enter a positive real number to compute logarithm: ");
    scanf("%lf", &y);
    printf("log(%.2lf) ≈ %.5lf\n", y, logarithm(y));

    return 0;
}
