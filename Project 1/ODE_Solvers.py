import numpy as np
from scipy.integrate import solve_ivp

## Euler's method for Damped Harmonic Oscillator
def euler_step(q, p, dt, k, m, b):
    """
    Performs one step of Euler's method to update the position (q) and momentum (p)
    for a damped harmonic oscillator.

    Parameters:
    q : float
        Current position.
    p : float
        Current momentum.
    dt : float
        Time step size.
    k : float
        Spring constant (controls the restoring force).
    m : float
        Mass of the oscillator.
    b : float
        Damping coefficient (controls the friction or resistance).

    Returns:
    q_new, p_new : float, float
        Updated position and momentum after one Euler step.
    """
    # Update position and momentum using Euler's method
    q_new = q + dt * (p / m)
    p_new = p - dt * (k * q + b * p)

    return q_new, p_new


## Runge-Kutta 4th order (RK4) for Damped Harmonic Oscillator
def rk4_step(q, p, dt, k, m, b):
    """
    Performs one step of the Runge-Kutta 4th order method to update position (q)
    and momentum (p) for a damped harmonic oscillator.

    Parameters:
    q : float
        Current position.
    p : float
        Current momentum.
    dt : float
        Time step size.
    k : float
        Spring constant (controls the restoring force).
    m : float
        Mass of the oscillator.
    b : float
        Damping coefficient (controls the friction or resistance).

    Returns:
    q_new, p_new : float, float
        Updated position and momentum after one RK4 step.
    """

    # Define a helper function for the oscillator's equations of motion
    def f(q, p):
        return np.array([p / m, -k * q - b * p])

    # Compute the four intermediate steps (k1, k2, k3, k4) for RK4 method
    k1 = dt * f(q, p)
    k2 = dt * f(q + k1[0] / 2, p + k1[1] / 2)
    k3 = dt * f(q + k2[0] / 2, p + k2[1] / 2)
    k4 = dt * f(q + k3[0], p + k3[1])

    # Update position and momentum using the weighted average of k1, k2, k3, k4
    q_new = q + (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
    p_new = p + (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6

    return q_new, p_new


## Solve using SciPy ODE solver
def scipy_ode_solver(t, y0, k, m, b):
    """
    Solves the damped harmonic oscillator equation using SciPy's `solve_ivp` ODE solver.

    Parameters:
    t : array-like
        Array of time points at which the solution is evaluated.
    y0 : list or array-like
        Initial conditions [q0, p0] for position and momentum.
    k : float
        Spring constant (controls the restoring force).
    m : float
        Mass of the oscillator.
    b : float
        Damping coefficient (controls the friction or resistance).

    Returns:
    q, p : ndarray, ndarray
        Arrays of positions and momenta at the time points in t.
    """

    # Define the system of differential equations for the damped harmonic oscillator
    def damped_oscillator(t, y):
        q, p = y
        dqdt = p / m
        dpdt = -k * q - b * p
        return [dqdt, dpdt]

    # Solve the ODE system over the time interval [t[0], t[-1]]
    sol = solve_ivp(damped_oscillator, [t[0], t[-1]], y0, t_eval=t)

    # Return the arrays of positions (q) and momenta (p) from the solution
    return sol.y[0], sol.y[1]


## General ODE Solver to run all methods
def solve_oscillator(method, q0, p0, t, k=1, m=1, b=0.1):
    """
    Solves the damped harmonic oscillator problem using the specified numerical method
    (Euler, RK4, or SciPy).

    Parameters:
    method : str
        The method to use ('euler', 'rk4', or 'scipy').
    q0 : float
        Initial position.
    p0 : float
        Initial momentum.
    t : array-like
        Array of time points at which the solution is evaluated.
    k : float, optional
        Spring constant (default is 1).
    m : float, optional
        Mass of the oscillator (default is 1).
    b : float, optional
        Damping coefficient (default is 0.1).

    Returns:
    q, p : ndarray, ndarray
        Arrays of positions and momenta at the time points in t.
    """
    # Create arrays to store the solution for position (q) and momentum (p)
    q, p = np.zeros(len(t)), np.zeros(len(t))  # Arrays to hold the results
    q[0], p[0] = q0, p0  # Set initial conditions for q and p

    # Loop through each time step and update q and p based on the selected method
    for i in range(1, len(t)):
        if method == 'euler':
            # Use Euler's method to update q and p
            q[i], p[i] = euler_step(q[i - 1], p[i - 1], t[i] - t[i - 1], k, m, b)
        if method == 'rk4':
            # Use RK4 method to update q and p
            q[i], p[i] = rk4_step(q[i - 1], p[i - 1], t[i] - t[i - 1], k, m, b)
        if method == 'scipy':
            # Use SciPy's ODE solver to compute the entire solution
            q, p = scipy_ode_solver(t, [q0, p0], k, m, b)
            #break  # SciPy solver computes everything in one step, no need for a loop

    return q, p  # Return the arrays of positions (q) and momenta (p)
