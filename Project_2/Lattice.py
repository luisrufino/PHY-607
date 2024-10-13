import numpy as np

class Lattice2D:
    """
    Class representing a 2D lattice for simulating particle motion.

    Parameters
    ----------
    lattice_size : int
        Size of the NxN lattice grid.
    """

    def __init__(self, lattice_size):
        self.lattice_size = lattice_size
        self.lattice = np.zeros((lattice_size, lattice_size))  # Initialize with zeros
        self.initialize_particles()

    def initialize_particles(self):
        """Initialize particles at each node with random kinetic energy."""
        self.lattice = np.random.rand(self.lattice_size, self.lattice_size)

    def get_neighbors(self, i, j):
        """Get the neighboring nodes of a given node (i, j)."""
        directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]  # Up, Down, Right, Left
        neighbors = [(i + di, j + dj) for di, dj in directions if 0 <= i + di < self.lattice_size and 0 <= j + dj < self.lattice_size]
        return neighbors

    def update_particle_energy(self, i, j, energy):
        """Update the kinetic energy of a particle at position (i, j)."""
        self.lattice[i, j] = energy

    def get_total_energy(self):
        """Return the total kinetic energy of the lattice."""
        return np.sum(self.lattice)
