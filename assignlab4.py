import math

def factorial(n):
    if n == 0 or n == 1:
        return 1
    else:
        return n * factorial(n - 1)

def factorial_error(n):
    recursive_result = factorial(n)
    builtin_result = math.factorial(n)

    relative_error = abs(recursive_result - builtin_result) / builtin_result 

    print(f"Relative error for factorial n={n}: {relative_error:}")

n = 10 
factorial_error(n)

import math


def exponential(z, terms=10):
    exp_value = sum((z ** n) / math.factorial(n) for n in range(terms))
    return exp_value


def calculate_exp(x):
    x0 = int(x)  
    z = x - x0   
    exp_z = exponential(z)  
    return exp_z * (math.e ** x0)  


def exp_error(x):
    custom_result = calculate_exp(x)
    builtin_result = math.exp(x)

    relative_error = abs(custom_result - builtin_result) / builtin_result 

    print(f"Relative error for e^{x}: {relative_error:}")


x = 2.5  
exp_error(x)
def f(s, x):
    return calculate_exp(s) - x

def f_prime(s):
    return calculate_exp(s)

def newtons_method(x, initial_guess, tolerance=1e-7, max_iterations=100):
    s = initial_guess
    for iteration in range(1, max_iterations + 1):
        s_new = s - f(s, x) / f_prime(s)
        if abs(s_new - s) < tolerance:
            return s_new
        s = s_new
    return s


def ln_error(x, initial_guess=1.0):
   
    custom_result = newtons_method(x, initial_guess)

    
    builtin_result = math.log(x)

    relative_error = abs(custom_result - builtin_result) / builtin_result 

    print(f"Relative error for ln({x}): {relative_error:}")
    print(f"Newton's method result: {custom_result}")
    print(f"Built-in math.log result: {builtin_result}")


x = 2.5  
ln_error(x)