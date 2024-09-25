import numpy as np
from scipy.integrate import quad


## Riemann Sum
def riemann_sum(f, a, b, N):
    """
    Approximates the integral of function f over the interval [a, b]
    using the Riemann sum method.

    Parameters:
    f : function
        The function to integrate.
    a : float
        The start of the interval.
    b : float
        The end of the interval.
    N : int
        The number of subintervals to use.

    Returns:
    float
        Approximate integral value using the Riemann sum.
    """
    # Create N evenly spaced points between a and b
    x = np.linspace(a, b, N)

    # Calculate the width of each subinterval
    dx = (b - a) / N

    # Evaluate the function at each point and sum the values, multiplied by dx
    return np.sum(f(x)) * dx


## Trapezoidal Rule
def trapezoidal_rule(f, a, b, N):
    """
    Approximates the integral of function f over the interval [a, b]
    using the trapezoidal rule.

    Parameters:
    f : function
        The function to integrate.
    a : float
        The start of the interval.
    b : float
        The end of the interval.
    N : int
        The number of subintervals to use (N points including boundaries).

    Returns:
    float
        Approximate integral value using the trapezoidal rule.
    """
    # Create N evenly spaced points between a and b
    x = np.linspace(a, b, N)

    # Calculate the width of each subinterval
    dx = (b - a) / (N - 1)

    # Apply the trapezoidal rule: the first and last points contribute half,
    # while the internal points contribute fully. Multiply by dx / 2 at the end.
    return (f(a) + f(b) + 2 * np.sum(f(x[1:-1]))) * dx / 2


## Simpson's Rule
def simpsons_rule(f, a, b, N):
    """
    Approximates the integral of function f over the interval [a, b]
    using Simpson's rule.

    Simpson's rule provides more accurate results than the trapezoidal rule,
    but it requires that N (number of subintervals) be even.

    Parameters:
    f : function
        The function to integrate.
    a : float
        The start of the interval.
    b : float
        The end of the interval.
    N : int
        The number of subintervals to use (must be even for Simpson's rule).

    Returns:
    float
        Approximate integral value using Simpson's rule.
    """
    # Ensure N is even for Simpson's rule; if odd, increment by 1
    if N % 2 == 1:
        N += 1

    # Create N evenly spaced points between a and b
    x = np.linspace(a, b, N)

    # Calculate the width of each subinterval
    dx = (b - a) / N

    # Apply Simpson's rule: alternating coefficients 4 and 2 are applied to
    # odd and even points, respectively, except the first and last points.
    return (dx / 3) * (f(a) + 4 * np.sum(f(x[1:-1:2])) + 2 * np.sum(f(x[2:-1:2])) + f(b))


## Compare with SciPy Quad
def scipy_integral(f, a, b):
    """
    Uses SciPy's quad function to compute the integral of function f
    over the interval [a, b] with high accuracy.

    This serves as a comparison with the manual numerical integration methods.

    Parameters:
    f : function
        The function to integrate.
    a : float
        The start of the interval.
    b : float
        The end of the interval.

    Returns:
    float
        The computed integral value using SciPy's quad method.
    """
    # The quad function returns both the integral result and an estimate of the error.
    result, error = quad(f, a, b)

    # Return only the integral result (ignoring the error estimate)
    return result
