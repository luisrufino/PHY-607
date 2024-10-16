from scipy.integrate import solve_ivp


def second_order_ode(t, y):
    """
    Define the second-order ODE for the particle's motion.

    Parameters
    ----------
    t : float
        Time variable.
    y : np.ndarray
        Position and velocity combined as a single array.

    Returns
    -------
    np.ndarray
        Derivative of position and velocity.
    """
    position = y[0]
    velocity = y[1]
    d_position_dt = velocity
    d_velocity_dt = -position  # Simple harmonic motion
    return [d_position_dt, d_velocity_dt]


def solve_motion(y0, t_span):
    """
    Solve the second-order ODE for a particle's trajectory.

    Parameters
    ----------
    y0 : list
        Initial position and velocity.
    t_span : tuple
        Start and end times.

    Returns
    -------
    sol : OdeSolution
        Solution of the ODE.
    """
    sol = solve_ivp(second_order_ode, t_span, y0, method='RK45')
    return sol
