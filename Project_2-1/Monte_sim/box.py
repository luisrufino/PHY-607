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
    temperature : float
        The temperature of the system (in Kelvin), which determines the strength of random forces for Brownian motion.
    """

    def __init__(self, size=10.0, temperature=300.0):
        self.size = size
        self.particles = []
        self.temperature = temperature
        self.volume = size ** 3  # Volume of the box
        self.k_B = 1.38e-23  # Boltzmann constant

        # Lists to store measurements over time
        self.pressure_over_time = []
        self.temperature_over_time = []
        self.energy_over_time = []
        self.entropy_over_time = []
    def add_particle(self, particle):
        """
        Add a particle to the simulation box.

        Parameters
        ----------
        particle : Particle
            The particle to add to the box.
        """
        self.particles.append(particle)

    def apply_brownian_motion(self, particle, dt):
        """
        Apply random forces to simulate Brownian motion.

        Parameters
        ----------
        particle : Particle
            The particle to which Brownian motion is applied.
        dt : float
            The time step of the simulation.
        """
        # Boltzmann constant (J/K)
        k_B = 1.38e-23

        # Standard deviation of random velocity perturbation
        std_dev = np.sqrt(2 * k_B * self.temperature * dt / particle.mass)

        # Add a random velocity perturbation (Gaussian-distributed)
        random_force = np.random.normal(0, std_dev, size=3)
        particle.velocity += random_force

    def handle_wall_collisions(self, particle):
        """
        Handle collisions with the walls of the box.

        If a particle hits a wall, its velocity in the corresponding direction is reversed.
        Adjust the particle's position if it exceeds the boundary (taking its radius into account).

        Returns True if a collision occurs, False otherwise.

        Parameters
        ----------
        particle : Particle
            The particle to check for wall collisions.
        """
        collision = False
        for i in range(3):  # Check x, y, z directions
            if particle.position[i] - particle.radius <= 0:
                particle.position[i] = particle.radius  # Keep particle inside the box
                particle.velocity[i] *= -1  # Reverse velocity direction
                collision = True
            elif particle.position[i] + particle.radius >= self.size:
                particle.position[i] = self.size - particle.radius  # Keep particle inside the box
                particle.velocity[i] *= -1  # Reverse velocity direction
                collision = True
        return collision  # Return whether collision occurred

    def handle_particle_collisions(self, particle1, particle2):
        """
        Handle collisions between two particles using elastic collision rules.

        Returns True if a collision occurs, False otherwise.

        Parameters
        ----------
        particle1 : Particle
            The first particle.
        particle2 : Particle
            The second particle.
        """
        distance = np.linalg.norm(particle1.position - particle2.position)

        # Check if the distance between particles is less than or equal to the sum of their radii
        if distance <= particle1.radius + particle2.radius:
            # Elastic collision between two particles
            v1, v2 = particle1.velocity, particle2.velocity
            m1, m2 = particle1.mass, particle2.mass
            pos_diff = particle1.position - particle2.position
            vel_diff = v1 - v2

            # Update velocities (elastic collision formula)
            particle1.velocity = v1 - (2 * m2 / (m1 + m2)) * (
                        np.dot(vel_diff, pos_diff) / np.dot(pos_diff, pos_diff)) * pos_diff
            particle2.velocity = v2 - (2 * m1 / (m1 + m2)) * (
                        np.dot(-vel_diff, -pos_diff) / np.dot(pos_diff, pos_diff)) * -pos_diff

            return True  # Collision occurred
        return False  # No collision

    def calculate_kinetic_energy(self):
        return sum(0.5 * p.mass * (p.velocity) ** 2 for p in self.particles)

    def calculate_temperature(self):
        avg_kinetic_energy = self.calculate_kinetic_energy() / len(self.particles)
        return (2 / 3) * avg_kinetic_energy / self.k_B

    def calculate_pressure(self):
        avg_kinetic_energy = self.calculate_kinetic_energy() / len(self.particles)
        return (2 / 3) * avg_kinetic_energy / self.volume

    def calculate_entropy(self):
        # Using velocity distribution to approximate entropy
        velocities = np.array([np.linalg.norm(p.velocity) for p in self.particles])
        hist, _ = np.histogram(velocities, bins=30, density=True)
        hist = hist[hist > 0]  # Remove zero probabilities
        return -self.k_B * np.sum(hist * np.log(hist))

    def step(self, dt):
        collisions = [False] * len(self.particles)

        for particle in self.particles:
            self.apply_brownian_motion(particle, dt)
            particle.update_position(dt)
            if self.handle_wall_collisions(particle):
                collisions[self.particles.index(particle)] = True

        for i in range(len(self.particles)):
            for j in range(i + 1, len(self.particles)):
                if self.handle_particle_collisions(self.particles[i], self.particles[j]):
                    collisions[i] = True
                    collisions[j] = True

        # Measure pressure, temperature, and energy
        pressure = self.calculate_pressure()
        temperature = self.calculate_temperature()
        energy = self.calculate_kinetic_energy()
        entropy = self.calculate_entropy()

        # Store measurements over time
        self.pressure_over_time.append(pressure)
        self.temperature_over_time.append(temperature)
        self.energy_over_time.append(energy)
        self.entropy_over_time.append(entropy)

        return collisions