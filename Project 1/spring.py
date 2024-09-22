import numpy as np
import scipy as sp
import pylab as plt

##Spring ODE

def Force(x):
    return -k/m * x

def explicit(x, v, t):
    a = Force(x)
    vf = (v + a*t)*m
    xf = x + (v * t)
    return xf, vf

def derivs(q,p,t, type): ## Derivatives from hamiltons equations
    if type == 1: ## Update the position
        return q + t*(p/m)
    if type == 0: ## Update the momentum
        return p - t*(k*q)

if __name__ == "__main__":
    dt = np.pi/10
    delta = 0.1
    time = np.arange(0, 2*np.pi, dt)

    ## Init conditions
    k = 1
    m = 1
    ## Lest see if this conserves the area
    s0 = np.zeros((len(time),2)) ## Row: time, column: q, and P
    s0[0] = np.array([1,0])
    s1 = np.zeros((len(time), 2))
    s1[0] = s0[0] + np.array([delta, 0])
    s2 = np.zeros((len(time), 2))
    s2[0] = s0[0] + np.array([0, delta])
    s3 = np.zeros((len(time), 2))
    s3[0] = s0[0] + np.array([delta, delta])

    ps = [s0, s1, s2, s3]
    for i in range(len(ps)):
        p_f = ps[i][:,1]
        q_f = ps[i][:,0]
        for t in range(1,len(time)):
            p_f[t] = derivs(q_f[t - 1], p_f[t - 1], dt, 0) ## Updating the momentum
            q_f[t] = derivs(q_f[t - 1], p_f[t], dt, 1) ## Updating the position



    particles = np.array([s0, s1, s2, s3, s0])
    for i in range(len(time)):
        a = np.array([particles[p][i] for p in range(particles.shape[0])])

        plt.plot(a[:, 0], a[:, 1], 'bo')  # Draw points

        # Connect only the outer points in the correct order
        plt.plot([a[0, 0], a[1, 0]], [a[0, 1], a[1, 1]], 'r-')  # Connect point 0 to point 1
        plt.plot([a[1, 0], a[3, 0]], [a[1, 1], a[3, 1]], 'r-')  # Connect point 1 to point 3
        plt.plot([a[3, 0], a[2, 0]], [a[3, 1], a[2, 1]], 'r-')  # Connect point 3 to point 2
        plt.plot([a[2, 0], a[0, 0]], [a[2, 1], a[0, 1]], 'r-')  # Connect point 2 to point 0
    plt.title(f"Symplectic method")
    plt.xlabel('Position')
    plt.ylabel('Momentum')
    plt.show()


    ## Do the something simlar for the Classic implicit method
    s0 = np.zeros((len(time),2)) ## Row: time, column: x, and v
    s0[0] = np.array([1,0])
    s1 = np.zeros((len(time), 2))
    s1[0] = s0[0] + np.array([delta, 0])
    s2 = np.zeros((len(time), 2))
    s2[0] = s0[0] + np.array([0, delta])
    s3 = np.zeros((len(time), 2))
    s3[0] = s0[0] + np.array([delta, delta])
    ps = [s0, s1, s2, s3]
    ## Applying the forward method
    for i in range(len(ps)):
        x_f = ps[i][:,0]
        v_f = ps[i][:,1]
        for t in range(1,len(time)):
            x_f[t], v_f[t] = calc_forward(x_f[t-1], v_f[t-1], dt)

    particles = np.array([s0, s1, s2, s3, s0])
    for i in range(len(time)):
        a = np.array([particles[p][i] for p in range(particles.shape[0])])
        print(a)

        plt.plot(a[:, 0], a[:, 1], 'bo')  # Draw points

        # Connect only the outer points in the correct order
        plt.plot([a[0, 0], a[1, 0]], [a[0, 1], a[1, 1]], 'r-')  # Connect point 0 to point 1
        plt.plot([a[1, 0], a[3, 0]], [a[1, 1], a[3, 1]], 'r-')  # Connect point 1 to point 3
        plt.plot([a[3, 0], a[2, 0]], [a[3, 1], a[2, 1]], 'r-')  # Connect point 3 to point 2
        plt.plot([a[2, 0], a[0, 0]], [a[2, 1], a[0, 1]], 'r-')  # Connect point 2 to point 0
    plt.title(f"explicit method")
    plt.xlabel('Position')
    plt.ylabel('Momentum')
    plt.show()


