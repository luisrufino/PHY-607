import numpy as np
import pylab as plt
import matplotlib.pyplot as plt
from Integral_est import riemann_sum, trapezoidal_rule, simpsons_rule, scipy_integral

## Gaussian function
def gaussian(x):
    '''
    Defines the Gaussian function, exp(-x^2)

    Parameters:
    x : float or ndarray
        The point(s) at which to evaluate the Gaussian function.

    Returns:
    float or ndarray
        The value of exp(-x^2) evaluated at x.
    '''
    return np.exp(-x**2)

## Integration limits and number of points
a, b = -10, 10  # Define the integration limits (approximating an infinite range for the Gaussian function)
N = 1000     # Number of bins

## Perform integration using different methods
methods = ['Riemann Sum', 'Trapezoidal', 'Simpson', 'SciPy']
"""
List of method names corresponding to the numerical integration techniques that
will be applied. This helps label the results when printing them out later.
"""

# Calculate the integrals using each method and store results in the 'integrals' list
integrals = [
    riemann_sum(gaussian, a, b, N),         # Riemann sum approximation
    trapezoidal_rule(gaussian, a, b, N),    # Trapezoidal rule approximation
    simpsons_rule(gaussian, a, b, N),       # Simpson's rule approximation
    scipy_integral(gaussian, a, b)          # SciPy's quad method (for comparison)
]

# Loop over the method names and the corresponding results, printing each one
for method, result in zip(methods, integrals):
    """
    The zip function pairs each method name with the corresponding integral result.
    This loop iterates through these pairs and prints the method name along with 
    the calculated integral result.
    """
    print(f'{method}: {result}')  # Print the method name and the computed integral

# NUM = 50
# slices = np.arange(1, NUM, 1)
# rs = []
# tr = []
# sr = []
#
# for N in slices:
#     ## Perform integration using different methods
#     methods = ['Riemann Sum', 'Trapezoidal', 'Simpson', 'SciPy']
#     """
#     List of method names corresponding to the numerical integration techniques that
#     will be applied. This helps label the results when printing them out later.
#     """
#
#     # Calculate the integrals using each method and store results in the 'integrals' list
#     integrals = [
#         riemann_sum(gaussian, a, b, N),  # Riemann sum approximation
#         trapezoidal_rule(gaussian, a, b, N),  # Trapezoidal rule approximation
#         simpsons_rule(gaussian, a, b, N),  # Simpson's rule approximation
#         #scipy_integral(gaussian, a, b)  # SciPy's quad method (for comparison)
#     ]
#
#     # Loop over the method names and the corresponding results, printing each one
#     for method, result in zip(methods, integrals):
#         """
#         The zip function pairs each method name with the corresponding integral result.
#         This loop iterates through these pairs and prints the method name along with
#         the calculated integral result.
#         """
#         if method == 'Riemann Sum':
#             rs.append(result)
#         if method == 'Trapezoidal':
#             tr.append(result)
#         if method == 'Simpson':
#             sr.append(result)
#         print(f'{method}: {result} for N: {N}')
#
# plt.plot(slices, tr, label = "Trapzodial rule")
# plt.plot(slices, rs, label = "Riemann Sum")
# plt.plot(slices, sr, label = "Simpson")
# plt.hlines(np.sqrt(np.pi),1, NUM, label = "Exact value", colors = "Black")
# plt.xlabel("Number of integration points")
# plt.title("Convergence of integration")
# plt.legend()
# plt.show()