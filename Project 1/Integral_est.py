import numpy as np
from scipy.integrate import quad

## Riemann Sum
def riemann_sum(f, a, b, N):
    x = np.linspace(a, b, N)
    dx = (b - a) / N
    return np.sum(f(x)) * dx

## Trapezoidal Rule
def trapezoidal_rule(f, a, b, N):
    x = np.linspace(a, b, N)
    dx = (b - a) / (N - 1)
    return (f(a) + f(b) + 2 * np.sum(f(x[1:-1]))) * dx / 2

## Simpson's Rule
def simpsons_rule(f, a, b, N):
    if N % 2 == 1:  ## Simpson's rule requires an even number of intervals
        N += 1
    x = np.linspace(a, b, N)
    dx = (b - a) / N
    return (dx / 3) * (f(a) + 4 * np.sum(f(x[1:-1:2])) + 2 * np.sum(f(x[2:-1:2])) + f(b))

## Compare with SciPy Quad
def scipy_integral(f, a, b):
    result, error = quad(f, a, b)
    return result
