from .box import Box
from .particle import Particle


def initialize_particles(n_particles, box_size):
    """
    Initialize particles randomly within the box.

    Parameters
    ----------
    n_particles : int
        The number of particles to initialize.
    box_size : float
        The size of the cubic box.

    Returns
    -------
    list
        A list of Particle objects initialized with random positions and velocities.
    """
    particles = []
    for _ in range(n_particles):
        position = np.random.rand(3) * box_size  # Random position within the box
        velocity = np.random.randn(3)  # Random velocity
        particle = Particle(position, velocity)
        particles.append(particle)
    return particles


def run_simulation(n_steps, dt, box_size, n_particles):
    """
    Run the simulation for a specified number of steps.

    Parameters
    ----------
    n_steps : int
        The number of time steps to run the simulation.
    dt : float
        Time step size.
    box_size : float
        Size of the simulation box.
    n_particles : int
        Number of particles in the simulation.
    """
    # Initialize the simulation box and particles
    box = Box(box_size)
    particles = initialize_particles(n_particles, box_size)

    for particle in particles:
        box.add_particle(particle)

    # Run the simulation
    for step in range(n_steps):
        box.step(dt)
        if step % 100 == 0:
            print(f"Step {step} completed")


import numpy as np


# Define the Maxwell-Boltzmann distribution for velocity components
def maxwell_boltzmann(v, m, T):
    """
    Maxwell-Boltzmann distribution for a given velocity component.

    Parameters
    ----------
    v : float
        Velocity component (v_x, v_y, or v_z).
    m : float
        Mass of the particle.
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Value of the Maxwell-Boltzmann distribution at velocity v.
    """
    k_B = 1.38e-23  # Boltzmann constant (J/K)
    return np.sqrt(m / (2 * np.pi * k_B * T)) * np.exp(-m * v ** 2 / (2 * k_B * T))


# Metropolis-Hastings algorithm for generating 3D velocities
def metropolis_hastings_3d(n_samples, m, T, proposal_std=1.0):
    """
    Use the Metropolis-Hastings algorithm to sample from the 3D Maxwell-Boltzmann distribution.

    Parameters
    ----------
    n_samples : int
        Number of samples to generate.
    m : float
        Mass of the particle (in kg).
    T : float
        Temperature (in Kelvin).
    proposal_std : float
        Standard deviation for the Gaussian proposal distribution.

    Returns
    -------
    np.ndarray
        Array of shape (n_samples, 3) containing the sampled velocities (v_x, v_y, v_z).
    """
    samples = np.zeros((n_samples, 3))  # To store sampled velocities
    current_v = np.zeros(3)  # Start at the origin (v_x = 0, v_y = 0, v_z = 0)

    for i in range(n_samples):
        # Propose a new velocity for each component
        proposed_v = current_v + np.random.normal(0, proposal_std, size=3)

        # Compute the acceptance ratio (target distribution at proposed / current)
        current_pdf = np.prod([maxwell_boltzmann(v, m, T) for v in current_v])
        proposed_pdf = np.prod([maxwell_boltzmann(v, m, T) for v in proposed_v])

        acceptance_ratio = proposed_pdf / current_pdf

        # Accept or reject the proposal
        if np.random.rand() < acceptance_ratio:
            current_v = proposed_v

        # Store the current velocity (either the proposed one if accepted, or the old one)
        samples[i] = current_v

    return samples


# Example usage
m = 4.65e-26  # Mass of an oxygen molecule in kg
T = 300  # Temperature in Kelvin
n_samples = 1000

velocities = metropolis_hastings_3d(n_samples, m, T)

# Plot the sampled velocities
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(velocities[:, 0], velocities[:, 1], velocities[:, 2], alpha=0.6)
ax.set_xlabel("v_x")
ax.set_ylabel("v_y")
ax.set_zlabel("v_z")
plt.title("Sampled Velocities Using Metropolis-Hastings")
plt.show()
