import numpy as np
from scipy.integrate import solve_ivp


##  Euler's method for Damped Harmonic Oscillator
def euler_step(q, p, dt, k, m, b):

    q_new = q + dt * (p / m)
    p_new = p - dt * (k * q + b * p)
    return q_new, p_new


##  Runge-Kutta 4th order for Damped Harmonic Oscillator
def rk4_step(q, p, dt, k, m, b):
    def f(q, p):
        return np.array([p / m, -k * q - b * p])

    k1 = dt * f(q, p)
    k2 = dt * f(q + k1[0] / 2, p + k1[1] / 2)
    k3 = dt * f(q + k2[0] / 2, p + k2[1] / 2)
    k4 = dt * f(q + k3[0], p + k3[1])

    q_new = q + (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
    p_new = p + (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6

    return q_new, p_new


##  Solve using Scipy ODE solver
def scipy_ode_solver(t, y0, k, m, b):
    def damped_oscillator(t, y):
        q, p = y
        dqdt = p / m
        dpdt = -k * q - b * p
        return [dqdt, dpdt]

    sol = solve_ivp(damped_oscillator, [t[0], t[-1]], y0, t_eval=t)
    return sol.y[0], sol.y[1]


##  General ODE Solver to run all methods
def solve_oscillator(method, q0, p0, t, k=1, m=1, b=0.1):
    q, p = np.zeros(len(t)), np.zeros(len(t)) ## ##  Creates arrays for both p and q of equal size
    q[0], p[0] = q0, p0 ## ##  Alocates the inital conditions

    for i in range(1, len(t)):
        if method == 'euler':
            q[i], p[i] = euler_step(q[i - 1], p[i - 1], t[i] - t[i - 1], k, m, b)
        if method == 'rk4':
            q[i], p[i] = rk4_step(q[i - 1], p[i - 1], t[i] - t[i - 1], k, m, b)
        if method == 'scipy':
            q, p = scipy_ode_solver(t, [q0, p0], k, m, b)
            break ## ##  Break here because the scipy solver will alocate all the elements of q and p when it's done.

    return q, p
