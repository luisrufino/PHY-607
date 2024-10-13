import numpy as np

class MonteCarloSimulation:
    def __init__(self, lattice):
        self.lattice = lattice
        self.energy_history = []  # To track energy over time

    def metropolis_hastings(self, i, j, temperature=1.0):
        """Perform the Metropolis-Hastings algorithm."""
        neighbors = self.lattice.get_neighbors(i, j)
        selected_neighbor = neighbors[np.random.randint(len(neighbors))]

        current_energy = self.lattice.lattice[i, j]
        proposed_energy = np.random.rand()  # Randomly propose a new energy
        delta_energy = proposed_energy - current_energy

        if delta_energy < 0 or np.random.rand() < np.exp(-delta_energy / temperature):
            self.lattice.update_particle_energy(i, j, proposed_energy)

        # Track energy after each move
        total_energy = self.lattice.get_total_energy()
        self.energy_history.append(total_energy)

    def time_to_thermalization(self, tolerance=1e-3):
        """
        Measure the time to thermalize by tracking energy fluctuations.

        Parameters
        ----------
        tolerance : float
            The energy fluctuation threshold to consider the system thermalized.

        Returns
        -------
        int
            The number of steps to thermalize.
        """
        for step in range(1, len(self.energy_history)):
            energy_diff = abs(self.energy_history[step] - self.energy_history[step - 1])
            if energy_diff < tolerance:
                return step  # Return the step where thermalization occurs
        return -1  # Return -1 if not yet thermalized
