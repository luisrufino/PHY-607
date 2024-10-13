from monte_carlo_simulation.lattice import Lattice2D
from monte_carlo_simulation.monte_carlo import MonteCarloSimulation

if __name__ == "__main__":
    # Create a 10x10 lattice
    lattice_size = 10
    lattice = Lattice2D(lattice_size)

    # Set up Monte Carlo simulation
    monte_carlo = MonteCarloSimulation(lattice)

    # Run Monte Carlo simulation for 1000 steps, updating all particles
    for _ in range(1000):
        monte_carlo.run_full_update()  # Update all particles at each step

    # Measure time to thermalization
    steps_to_thermalization = monte_carlo.time_to_thermalization()
    print(f"Time to thermalization: {steps_to_thermalization}")
