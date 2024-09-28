import numpy as np
import matplotlib.pyplot as plt
from ODE_Solvers import solve_oscillator

## Parameters for the Damped Harmonic Oscillator
k, m = 1, 1  # Spring constant (k) and mass (m)
q0, p0 = 2, 0  # Initial conditions: initial position (q0) and momentum (p0)
v0 = p0 / m  # Initial velocity
t = np.linspace(0, 10 * np.pi, 1000)  # Time array from 0 to 10π


# Exact solution for under-damped case
def ext_under(q0, v0, b, m, k, t):
    beta = b / (2 * m)
    w0 = np.sqrt(k / m)
    w1 = np.sqrt(w0 ** 2 - beta ** 2)

    # Calculate the constants based on initial conditions
    A = q0
    B = (v0 + beta * q0) / w1
    return A * np.exp(-beta * t) * np.cos(w1 * t) + B * np.exp(-beta * t) * np.sin(w1 * t)


# Exact solution for over-damped case
def ext_over(q0, v0, b, m, k, t):
    beta = b / (2 * m)
    w0 = np.sqrt(k / m)
    r1 = -beta + np.sqrt(beta ** 2 - w0 ** 2)
    r2 = -beta - np.sqrt(beta ** 2 - w0 ** 2)

    # Solve for the constants C1 and C2 based on initial conditions
    C1 = (v0 - r2 * q0) / (r1 - r2)
    C2 = q0 - C1
    return C1 * np.exp(r1 * t) + C2 * np.exp(r2 * t)


# Exact solution for critically-damped case
def ext_critical(q0, v0, b, m, k, t):
    beta = b / (2 * m)

    # Calculate the constants for the critically damped case
    A = q0
    B = v0 + beta * q0
    return (A + B * t) * np.exp(-beta * t)


## Plot for each damping scenario
b_u = 0.01
b_c = 2
b_o = 10
damping_cases = {
    f'Under-damped (b = {b_u})': {'b': b_u, 'exact': ext_under},
    f'Critically-damped (b = {b_c})': {'b': b_c, 'exact': ext_critical},
    f'Over-damped (b = {b_o})': {'b': b_o, 'exact': ext_over}
}

methods = ['euler', 'rk4', 'scipy']  # List of numerical methods to compare

for label, case in damping_cases.items():
    b = case['b']
    exact_solution = case['exact'](q0, v0, b, m, k, t)

    plt.figure(figsize=(10, 6))

    # Plot the exact solution
    plt.plot(t, exact_solution, label="Exact solution", color="black", linestyle='dashed')

    # Plot the solutions from different methods
    for method in methods:
        q, p = solve_oscillator(method, q0, p0, t, k, m, b)
        plt.plot(t, q, label=method)

    # Customize the plot
    plt.title(f"Damped Harmonic Oscillator - {label}")
    plt.xlabel('Time (s)')
    plt.ylabel('Position (q)')
    plt.legend()

    # Save the plot for LaTeX
    plt.savefig(f"{label.replace(' ', '_').replace('=', '_')}.png")  # Save each plot as a PNG file

## Plot for each damping scenario
b_values = [0.01, 2, 10]  # Underdamped, critically damped, overdamped
damping_labels = ['Underdamped', 'Critically Damped', 'Overdamped']
methods = ['euler', 'rk4', 'scipy']  # List of numerical methods to compare

# Prepare for plotting differences
plt.figure(figsize=(15, 15))

for idx, b in enumerate(b_values):
    t_exact = np.linspace(0, 30, 1000)  # Time array for exact solution

    # Select the appropriate exact solution function based on the damping case
    if idx == 0:
        exact_solution = ext_under(q0, v0, b, m, k, t_exact)
    elif idx == 1:
        exact_solution = ext_critical(q0, v0, b, m, k, t_exact)
    elif idx == 2:
        exact_solution = ext_over(q0, v0, b, m, k, t_exact)

    # Loop through each numerical method to calculate and plot differences
    for method in methods:
        t = np.arange(0, 30, 0.1)  # Define a time array for numerical methods
        q, _ = solve_oscillator(method, q0, p0, t, k, m, b)

        # Interpolate exact solution for the current time points
        exact_interp = np.interp(t, t_exact, exact_solution)  # Interpolating exact solution
        difference = np.abs(q - exact_interp)  # Calculate absolute difference

        # Plot the differences
        plt.subplot(3, 1, idx + 1)
        plt.plot(t, difference, label=f'{method} Difference')

    # Customize the plot for the current damping case
    plt.title(f'Difference Between Exact Solution and Numerical Methods - {damping_labels[idx]} Damping')
    plt.xlabel('Time (s)')
    plt.ylabel('Absolute Difference')
    plt.legend()
    plt.grid()

# Save the figure
plt.savefig('damped_harmonic_oscillator_differences.png', dpi=300)  # Save the plot as a PNG file
plt.close()
