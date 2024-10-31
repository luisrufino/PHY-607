from monte_carlo_sim.Lattice import Lattice2D
from monte_carlo_sim.MC import MonteCarloSimulation

if __name__ == "__main__":
    # Create a 100x100 lattice (or larger lattice)
    lattice_size = 100
    lattice = Lattice2D(lattice_size)

    # Set up Monte Carlo simulation
    monte_carlo = MonteCarloSimulation(lattice)

    # Run Monte Carlo simulation for 1000 steps, updating all particles
    for _ in range(1000):
        monte_carlo.run_full_update(temperature=1.5)  # Update all particles at each step

    # Measure time to thermalization
    steps_to_thermalization = monte_carlo.time_to_thermalization(tolerance=1e-5)
    print(f"Time to thermalization: {steps_to_thermalization}")

    # Plot the total energy history over time
    monte_carlo.plot_energy_history()
