import numpy as np

class MonteCarloSimulation:
    """
    Class for running a Monte Carlo simulation on the lattice.
    """

    def __init__(self, lattice):
        self.lattice = lattice

    def metropolis_hastings(self, i, j, temperature=1.0):
        """
        Perform the Metropolis-Hastings algorithm on a particle at position (i, j).
        """
        neighbors = self.lattice.get_neighbors(i, j)
        selected_neighbor = neighbors[np.random.randint(len(neighbors))]

        current_energy = self.lattice.lattice[i, j]
        proposed_energy = np.random.rand()  # Randomly propose a new energy
        delta_energy = proposed_energy - current_energy

        if delta_energy < 0 or np.random.rand() < np.exp(-delta_energy / temperature):
            self.lattice.update_particle_energy(i, j, proposed_energy)
