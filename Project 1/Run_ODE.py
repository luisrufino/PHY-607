import numpy as np
import matplotlib.pyplot as plt
from ODE_Solvers import solve_oscillator

##  Parameters
k, m, b = 1, 1, 0.1
q0, p0 = 1, 0  ##  Initial conditions
t = np.linspace(0, 10 * np.pi, 1000)

##  Solve using different methods
methods = ['euler', 'rk4', 'scipy']
for method in methods:
    q, p = solve_oscillator(method, q0, p0, t, k, m, b)

    plt.plot(q, p, label=method)

plt.title("Phase Space of Damped Harmonic Oscillator")
plt.xlabel('Position (q)')
plt.ylabel('Momentum (p)')
plt.legend()
plt.show()
