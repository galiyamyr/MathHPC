import numpy as np
import matplotlib.pyplot as plt

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

def matrixA(x_points, degree=None):
    """Constructs the matrix for polynomial interpolation or least squares approximation."""
    if degree is None:
        degree = len(x_points) - 1
    n = len(x_points)
    A = np.zeros((n, degree + 1), dtype=float)
    for i in range(n):   
        for j in range(degree + 1):
            A[i, j] = x_points[i] ** j
    return A


def evaluate_polynomial(coefficients, x):
    """Evaluates the polynomial at a given x."""
    return sum(coefficients[i] * (x**i) for i in range(len(coefficients)))

def LeastSquareApprox(x, f, n):
    """Finds the least squares approximation of data (x, f) by a polynomial of degree <= n."""
    A = matrixA(x, n)  # Construct matrix using existing function
    ATA = A.T @ A  # Compute A^T * A
    ATf = A.T @ f  # Compute A^T * f
    coeffs = gauss_elimination(ATA, ATf)  # Solve the normal equations using Gaussian elimination
    return coeffs

def evaluate_polynomial(coefficients, x):
    """Evaluates the polynomial at a given x."""
    return sum(coefficients[i] * (x**i) for i in range(len(coefficients)))
def polynomial_string(coefficients):
    """Returns the string representation of a polynomial."""
    terms = [f"{coefficients[i]} x^{i}" if i > 0 else f"{coefficients[i]}" for i in range(len(coefficients))]
    return " + ".join(terms)

def main():
    """Main program for least squares approximation."""
    # Application: Least squares approximation of cos(x) on [-pi, pi]
    x_values = np.linspace(-np.pi, np.pi, 51)  # 51 nodes in [-pi, pi]
    f_values = np.cos(x_values)  # Compute cos(x) at these nodes

    # Find the polynomial of degree <= 5
    coefficients = LeastSquareApprox(x_values, f_values, 5)

    # Plot comparison of exact function and least squares approximation
    x_plot = np.linspace(-np.pi, np.pi, 200)
    y_exact = np.cos(x_plot)
    y_approx = [evaluate_polynomial(coefficients, x) for x in x_plot]

    plt.figure(figsize=(8, 6))
    plt.plot(x_plot, y_exact, label="cos(x)", linestyle="--", color="blue")
    plt.plot(x_plot, y_approx, label="Least Squares Approximation", linestyle="-", color="red")
    plt.scatter(x_values, f_values, color="black", marker="o", label="Data Points")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("Least Squares Approximation of cos(x)")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()
