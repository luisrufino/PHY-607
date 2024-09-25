import numpy as np
import matplotlib.pyplot as plt
from ODE_Solvers import solve_oscillator

## Parameters for the Damped Harmonic Oscillator
k, m, b = 1, 1, 0.1  # k: spring constant, m: mass, b: damping coefficient
q0, p0 = 1, 0        # Initial conditions: initial position (q0) and initial momentum (p0)
t = np.linspace(0, 10 * np.pi, 1000)  # Time array from 0 to 10π with 1000 points

## Solve using different methods
methods = ['euler', 'rk4', 'scipy']  # List of numerical methods to compare

for method in methods:
    """
    For each method ('euler', 'rk4', or 'scipy'), solve the damped harmonic oscillator
    and obtain the arrays of position (q) and momentum (p) over time.
    """
    q, p = solve_oscillator(method, q0, p0, t, k, m, b)  # Solve using the selected method

    # Plot the phase space trajectory (position vs. momentum)
    plt.plot(q, p, label=method)  # Plot with a label for the method

## Plotting phase space for the damped harmonic oscillator
plt.title("Phase Space of Damped Harmonic Oscillator")  # Title for the plot
plt.xlabel('Position (q)')  # Label for the x-axis (Position)
plt.ylabel('Momentum (p)')  # Label for the y-axis (Momentum)
plt.legend()  # Display a legend showing the method used for each curve
plt.show()  # Show the plot
