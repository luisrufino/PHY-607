import numpy as np
from integral_est import riemann_sum, trapezoidal_rule, simpsons_rule, scipy_integral

##  Gaussian function
def gaussian(x):
    return np.exp(-x**2)

##  Integration limits and number of points
a, b = -10, 10  ##  Finite limits for approximation
N = 1000  ##  Number of intervals

##  Perform integration using different methods
methods = ['Riemann Sum', 'Trapezoidal', 'Simpson', 'SciPy']
integrals = [
    riemann_sum(gaussian, a, b, N),
    trapezoidal_rule(gaussian, a, b, N),
    simpsons_rule(gaussian, a, b, N),
    scipy_integral(gaussian, a, b)
]

for method, result in zip(methods, integrals):
    print(f'{method}: {result}')
