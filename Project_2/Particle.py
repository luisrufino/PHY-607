class Particle:
    """
    Class representing a particle on the lattice.

    Parameters
    ----------
    position : tuple
        Initial position (x, y) of the particle.
    kinetic_energy : float
        Initial kinetic energy of the particle.
    """

    def __init__(self, position, kinetic_energy=0.0):
        self.position = position
        self.kinetic_energy = kinetic_energy

    def move(self, new_position):
        """Move the particle to a new position."""
        self.position = new_position

    def update_energy(self, new_energy):
        """Update the particle's kinetic energy."""
        self.kinetic_energy = new_energy
