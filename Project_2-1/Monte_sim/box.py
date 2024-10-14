import numpy as np
from .particle import Particle

class Box:
    """
    Class representing the 3D box in which particles are simulated.

    Attributes
    ----------
    size : float
        The size of the cubic box.
    particles : list
        A list of Particle objects in the box.
    """

    def __init__(self, size=10.0):
        self.size = size
        self.particles = []

    def add_particle(self, particle):
        """
        Add a particle to the simulation box.

        Parameters
        ----------
        particle : Particle
            The particle to add to the box.
        """
        self.particles.append(particle)

    def handle_wall_collisions(self, particle):
        """
        Handle collisions with the walls of the box.

        If a particle hits a wall, its velocity in the corresponding direction is reversed.

        Parameters
        ----------
        particle : Particle
            The particle to check for wall collisions.
        """
        for i in range(3):  # Check x, y, z directions
            if particle.position[i] - particle.radius <= 0 or particle.position[i] + particle.radius >= self.size:
                particle.velocity[i] *= -1  # Reverse velocity direction

    def handle_particle_collisions(self, particle1, particle2):
        """
        Handle collisions between two particles using elastic collision rules.

        Parameters
        ----------
        particle1 : Particle
            The first particle.
        particle2 : Particle
            The second particle.
        """
        distance = np.linalg.norm(particle1.position - particle2.position)
        if distance <= particle1.radius + particle2.radius:
            # Elastic collision between two particles
            v1, v2 = particle1.velocity, particle2.velocity
            m1, m2 = particle1.mass, particle2.mass
            pos_diff = particle1.position - particle2.position
            vel_diff = v1 - v2

            # Update velocities (elastic collision formula)
            particle1.velocity = v1 - (2 * m2 / (m1 + m2)) * (np.dot(vel_diff, pos_diff) / np.dot(pos_diff, pos_diff)) * pos_diff
            particle2.velocity = v2 - (2 * m1 / (m1 + m2)) * (np.dot(-vel_diff, -pos_diff) / np.dot(pos_diff, pos_diff)) * -pos_diff

    def step(self, dt):
        """
        Perform one step of the simulation by updating particle positions and handling collisions.

        Parameters
        ----------
        dt : float
            Time step for the simulation.
        """
        for particle in self.particles:
            particle.update_position(dt)
            self.handle_wall_collisions(particle)

        # Check for particle-particle collisions
        for i in range(len(self.particles)):
            for j in range(i + 1, len(self.particles)):
                self.handle_particle_collisions(self.particles[i], self.particles[j])
