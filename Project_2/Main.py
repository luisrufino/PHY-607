from monte_carlo_simulation.lattice import Lattice2D
from monte_carlo_simulation.monte_carlo import MonteCarloSimulation

if __name__ == "__main__":
    # Create a 10x10 lattice
    lattice_size = 10
    lattice = Lattice2D(lattice_size)

    # Set up Monte Carlo simulation
    monte_carlo = MonteCarloSimulation(lattice)

    # Example: Measure distance traveled for a particle starting at (0,0)
    start_position = (0, 0)
    current_position = (5, 5)  # Assume the particle moved to (5, 5)
    distance_traveled = lattice.measure_average_distance(current_position, start_position)
    print(f"Distance traveled from {start_position} to {current_position}: {distance_traveled}")

    # Run the Metropolis-Hastings algorithm
    monte_carlo.metropolis_hastings(5, 5)  # Run for a particle at position (5,5)

    # Measure time to thermalization
    steps_to_thermalization = monte_carlo.time_to_thermalization()
    print(f"Time to thermalization: {steps_to_thermalization}")
