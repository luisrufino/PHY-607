import numpy as np
import scipy as sp
import pylab as plt


def calc_changes(pos, vel, force, time_step):
    """
    Calculate the changes in position and velocity.

    :param pos: Position vector (x, y)
    :param vel: Velocity vector (x, y)
    :param force: Force vector (x, y)
    :param time_step: Time step for the update
    :return: Updated position and velocity
    """
    pos_final = pos + vel * time_step + 0.5 * force * time_step ** 2
    vel_final = vel + force * time_step
    return pos_final, vel_final


def calc_total_force(vel):
    """
    Calculate total force (gravity + drag) on the object.

    :param vel: Velocity vector (x, y)
    :param g: Gravitational constant
    :param drag_const: Drag coefficient
    :return: Total force vector and its magnitude
    """
    force_grav = np.array([0, -g])
    force_drag = drag_const * vel ** 2
    total_force = force_grav + force_drag
    force_mag = np.linalg.norm(total_force)
    return total_force, force_mag


def calc_energy(pos, vel):
    """
    Calculate potential and kinetic energy of the system.

    :param pos: Current position (x, y)
    :param vel: Current velocity (x, y)
    ## The above parameters contain the changes in position and velocity
    :param mass: Mass of the object
    :param g: Gravitational constant
    :return: Potential energy, Kinetic energy
    """
    pe = g * mass * pos  # Using only vertical displacement for PE
    ke = 0.5 * mass * vel ** 2
    return pe, ke


if __name__ == "__main__":
    ## Define Global Variables, all units are in SI units
    g = 9.8
    mass = 1000
    drag_const = 0.01
    ## Create time array
    time = np.arange(0, 10, 0.00001)
    ## Create init parameters
    init_pos = np.array([0, 10])
    init_vel = np.array([0, 0])
    ## pos and vel containg the changes in the position and velocity
    pos = init_pos
    vel = init_vel

    pos_hist, vel_hist, pe_hist, ke_hist = [init_pos], [init_vel], [], []

    for t in time:
        force_t, _ = calc_total_force(vel)
        pos, vel = calc_changes(pos, vel, force_t, t)
        pe, ke = calc_energy(pos, vel,)

        pos_hist.append(pos)
        vel_hist.append(vel)
        pe_hist.append(pe)
        ke_hist.append(ke)

        print(f"At time step: {t}\n"
              f"Total force: {force_t}\n"
              f"Position: {pos}, Velocity: {vel}\n"
              f"Potential energy: {pe}, Kinetic energy: {ke}\n"
              f"---------------------------------\n")

        if pos[1] < 0:  # Cow hits the floor
            print(f"Cow passed the floor at time {t}. Simulation ends.")
            print(f"Final velocity: {vel_hist[-2]}")
            break

## Plot history of changes of position and velocity

