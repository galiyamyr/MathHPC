import numpy as np
import matplotlib.pyplot as plt

# Function to evaluate a polynomial given coefficients
def evaluate_polynomial(coeffs, x):
    y = np.zeros_like(x)
    for i, c in enumerate(coeffs):
        y += c * x**i
    return y

# Load coefficients
coeffs_m = np.loadtxt("coeffs_m.txt")
coeffs_n = np.loadtxt("coeffs_n.txt")

# Generate x values for smooth plotting
x_vals = np.linspace(0, 1, 200)
y_true = np.sin(2 * np.pi * x_vals)

# Compute polynomial values
y_m = evaluate_polynomial(coeffs_m, x_vals)
y_n = evaluate_polynomial(coeffs_n, x_vals)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(x_vals, y_true, 'k-', label=r'$f(x) = \sin(2\pi x)$', linewidth=2)
plt.plot(x_vals, y_m, 'r--', label=f'$P_5(x)$ (Least Squares)', linewidth=2)
plt.plot(x_vals, y_n, 'b-.', label=f'$P_{{19}}(x)$ (Interpolation)', linewidth=2)

# Highlight original points
x_data = np.linspace(0, 1, 20)
y_data = np.sin(2 * np.pi * x_data)
plt.scatter(x_data, y_data, color='black', marker='o', label='Data Points')

plt.xlabel("x")
plt.ylabel("y")
plt.title("Least Squares Approximation vs Interpolation")
plt.legend()
plt.grid()
plt.show()
plt.savefig("polynomial_plot.png")  # Save the plot
print("Plot saved as 'polynomial_plot.png'. Open it manually to view.")
