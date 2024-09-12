import numpy as np
import scipy as sp
import pylab as plt
import matplotlib.pyplot as plt2
import pandas as pd

def change_sign(tmp, tmp1):
    if tmp > 0:
        return -1 * abs(tmp1) ## Makes it negative
    if tmp < 0:
        return abs(tmp1) ## Makes it positive


def calc_changes(pos, vel, force, time_step):
    """
    Calculate the changes in position and velocity.

    :param pos: Position vector (x, y)
    :param vel: Velocity vector (x, y)
    :param force: Force vector (x, y)
    :param time_step: Time step for the update
    :return: Updated position and velocity
    """
    pos = np.array(pos)
    vel = np.array(vel)
    force = np.array(force)
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
    vel = np.array(vel)
    force_grav = np.array([0, -g])

    ## Drag should always be aposing the motion of the object
    force_drag = (drag_const/mass) * vel ** 2
    if vel[0] > 0:
        force_drag[0] = change_sign(vel[0], force_drag[0])
        # print(f"Positive velocity in x: {vel[0]}")
        # print(force_drag[1])
    else:
        force_drag[0] = change_sign(vel[0],force_drag[0])
        # print(f"Negative velocity in x: {vel[0]}")
        # print(force_drag[1])
    if vel[1] > 0:
        force_drag[1] = change_sign(vel[1], force_drag[1])
    else:
        force_drag[1] = change_sign(vel[1], force_drag[1])

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
    pos = np.array(pos)
    vel = np.array(vel)
    pe = g * mass * pos[1]  # Using only vertical displacement for PE
    ke = 0.5 * mass * vel ** 2
    return pe, ke





if __name__ == "__main__":
    ## Define Global Variables, all units are in SI units
    g = 9.8
    mass = 1000
    drag_const = 0.0
    ## Create time array
    time = np.arange(0, 1, 0.001)
    ## Create init parameters
    init_pos = np.array([0, 1000])
    init_vel = np.array([1, 100])
    ## pos and vel containg the changes in the position and velocity
    pos = init_pos
    vel = init_vel

    pos_hist, vel_hist, pe_hist, ke_hist = [init_pos], [init_vel], [], []
    time_hist = []
    vel_yarr = []
    vel_xarr = []
    count = 0
    for t in time:
        force_t, _ = calc_total_force(vel)
        pos, vel = calc_changes(pos, vel, force_t, t)
        pe, ke = calc_energy(pos, vel,)

        pos_hist.append(pos)
        vel_hist.append(vel)
        pe_hist.append(pe)
        ke_hist.append(ke)
        vel_xarr.append(vel[0])
        vel_yarr.append(vel[1])


        count += 1
        time_hist.append(t)
        print(f"At time step: {t}\n"
              f"Total force: {force_t}\n"
              f"Position: {pos}, Velocity: {vel}\n"
              f"Potential energy: {pe}, Kinetic energy: {ke}\n"
              f"---------------------------------\n")

        if pos[1] < 0:  # Cow hits the floor
            print(f"Cow passed the floor at time {t}. Simulation ends.")
            print(f"Final velocity: {vel_hist[count - 1]}")
            break

## Writing time, px, py to a file
pos_hist = np.array(pos_hist)
ke_hist = np.array(ke_hist)
time_hist = np.array(time_hist)


with open('data.txt', 'w') as f:
    for i in range(len(time_hist)):
        time = time_hist[i]
        x = pos_hist[:,0][i]
        y = pos_hist[:,1][i]
        stuff = f"{time:.4f}, {x}, {y}\n"
        f.write(stuff)

## Plot the potential energy and kinetic energy over time

plt.plot(time_hist[:count - 1], pe_hist[:count - 1], 'o', label = "Potential Energy")
plt.plot(time_hist[:count - 1], ke_hist[:,0][:count - 1], 'o', label = "Kinetic Energy - X")
plt.plot(time_hist[:count - 1], ke_hist[:,1][:count - 1], 'o', label = "Kinetic Energy - Y")
plt.legend()
plt.title("Energy vs Time")
plt.ylabel("Energy")
plt.xlabel("Time")
plt.show()


plt.plot(pos_hist[:,1][:count - 1], pe_hist[:count - 1], 'o', label = "Potential Energy")
#plt.plot(pos_hist[:,0][:count - 1], ke_hist[:,0][:count - 1], 'o', label = "Kinetic Energy - X")
plt.plot(pos_hist[:,1][:count - 1], ke_hist[:,1][:count - 1], 'o', label = "Kinetic Energy - Y")
plt.legend()
plt.title("Energy vs Distance")
plt.ylabel("Energy")
plt.xlabel("Distance")
plt.show()

plt.plot(time_hist,vel_yarr, 'x')
plt.title("Velocity in Y vs Time")
plt.ylabel("Velocity")
plt.xlabel("Time")
plt.show()


plt.plot(time_hist,np.multiply(mass,vel_xarr), 'x')
plt.title("Momentum in X vs Time")
plt.ylabel("Momentum")
plt.xlabel("Time")
plt.show()


plt.plot(time_hist,np.multiply(mass,vel_yarr), 'x')
plt.title("Momentum in Y vs Time")
plt.ylabel("Momentum")
plt.xlabel("Time")
plt.show()


## Init conditions
