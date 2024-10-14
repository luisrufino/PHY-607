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
