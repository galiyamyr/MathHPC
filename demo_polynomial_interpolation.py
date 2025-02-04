import numpy as np

def gauss_elimination(A, b):
    """Solves Ax = b using Gaussian elimination."""
    n = len(b)
    A = A.astype(float)
    b = b.astype(float)
    
    # Forward elimination
    for k in range(n-1):
        for i in range(k+1, n):
            factor = A[i, k] / A[k, k]
            b[i] = b[i] - factor * b[k]
            for j in range(k, n):
                A[i, j] = A[i, j] - factor * A[k, j]
    
    # Back substitution
    x = np.zeros(n, dtype=float)
    x[n-1] = b[n-1] / A[n-1, n-1]
    for i in range(n-2, -1, -1):
        sum_val = b[i]
        for j in range(i+1, n):
            sum_val = sum_val - A[i, j] * x[j]
        x[i] = sum_val / A[i, i]
    
    return x

# Define input points
x_points = np.array([-0.1, -0.02, 0.02, 0.1], dtype=float)
y_points = np.cos(x_points)

def matrixA(x_points):
    """Constructs the  matrix."""
    n = len(x_points)
    A = np.zeros((n, n), dtype=float)
    for i in range(n):   
        for j in range(n):
            A[i, j] = x_points[i] ** j
    return A

def cubic_polynomial_string(coefficients):
    terms = [f"{coefficients[i]} x^{i}" if i > 0 else f"{coefficients[i]}" for i in range(len(coefficients))]
    return " + ".join(terms)

def cubic_polynomial(coefficients, x):
    return sum(coefficients[i] * (x**i) for i in range(len(coefficients)))

# Compute matrix
A = matrixA(x_points)
coefficients = gauss_elimination(A, y_points)
print("Polynomial Coefficients:", coefficients)
print("Cubic Polynomial:", cubic_polynomial_string(coefficients))

# Compute maximum error at given points
errors = [abs(np.cos(x) - cubic_polynomial(coefficients, x)) for x in x_points]
max_error = max(errors)
print("Maximum Error:", max_error)
