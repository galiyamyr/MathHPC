def factorial(n):
    if n == 0 or n == 1:
        return 1
    else:
        return n * factorial(n - 1)

def exponential(z, terms=10):
    exp_value = 0
    for n in range(terms):
        exp_value += (z ** n) / factorial(n)
    return exp_value

def calculate_exp(x):
    x0 = int(x)
    z = x - x0
    exp_z = exponential(z)
    return exp_z * (2.718281828459045 ** x0)

def f(s, x):
    return calculate_exp(s) - x

def f_prime(s):
    return calculate_exp(s)

def newtons_method(x, initial_guess, tolerance=1e-7, max_iterations=100):
    s = initial_guess
    for iteration in range(1, max_iterations + 1):
        s_new = s - f(s, x) / f_prime(s)
        print(f"Iteration {iteration}: s = {s_new}")
        if abs(s_new - s) < tolerance:
            print(f"Converged after {iteration} iterations.")
            return s_new
        s = s_new
    print(f"Did not converge after {max_iterations} iterations.")
    return s

x = 2.5  
initial_guess = x
result = newtons_method(x, initial_guess)
print(f"\nThe solution to f(s) = e^s - {x} is s ≈ {result}")
