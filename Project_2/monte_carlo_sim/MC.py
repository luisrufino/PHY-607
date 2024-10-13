import numpy as np


class MonteCarloSimulation:
    def __init__(self, lattice):
        self.lattice = lattice
        self.energy_history = []  # To track energy over time

    def metropolis_hastings(self, temperature=1.0):
        """
        Perform the Metropolis-Hastings algorithm to update the energy of a random particle.
        This method involves energy exchange between a particle and a random neighbor.

        Parameters
        ----------
        temperature : float, optional
            Temperature of the system (default is 1.0).
        """
        # Randomly select a particle on the lattice
        i, j = np.random.randint(0, self.lattice.lattice_size, size=2)

        # Get neighbors of the current particle
        neighbors = self.lattice.get_neighbors(i, j)

        # Select a random neighboring node to exchange energy with
        ni, nj = neighbors[np.random.randint(len(neighbors))]

        # Current and neighbor energies
        current_energy = self.lattice.lattice[i, j]
        neighbor_energy = self.lattice.lattice[ni, nj]

        # Propose new energy values (e.g., a random fluctuation in energy)
        proposed_current_energy = current_energy + np.random.uniform(-0.1, 0.1)
        proposed_neighbor_energy = neighbor_energy + np.random.uniform(-0.1, 0.1)

        # Calculate energy differences for both particles
        delta_energy_current = proposed_current_energy - current_energy
        delta_energy_neighbor = proposed_neighbor_energy - neighbor_energy

        # Metropolis acceptance criterion (for both particles)
        if (delta_energy_current < 0 or np.random.rand() < np.exp(-delta_energy_current / temperature)) and \
                (delta_energy_neighbor < 0 or np.random.rand() < np.exp(-delta_energy_neighbor / temperature)):
            # Accept the new energies and update the lattice
            self.lattice.update_particle_energy(i, j, proposed_current_energy)
            self.lattice.update_particle_energy(ni, nj, proposed_neighbor_energy)

        # Track the total energy after the update
        total_energy = self.lattice.get_total_energy()
        self.energy_history.append(total_energy)

    def run_full_update(self, temperature=1.0):
        """
        Update every particle in the lattice during each Monte Carlo step.

        Parameters
        ----------
        temperature : float, optional
            Temperature of the system (default is 1.0).
        """
        for i in range(self.lattice.lattice_size):
            for j in range(self.lattice.lattice_size):
                self.metropolis_hastings(i, j, temperature=temperature)
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
