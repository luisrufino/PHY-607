import sys
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Add project root to sys.path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

# Import from Monte_sim
from Monte_sim.box import Box  # Import from box.py
from Monte_sim.particle import Particle  # Import from particle.py

# Convert particle's radius to marker size for the plot
def radius_to_marker_size(radius):
    return radius * 1000  # Customize this scaling factor to match the plot size

def run_simulation(n_steps, dt, box_size, n_particles, plot=True):
    # Initialize the box and particles
    box = Box(box_size)
    particles = [Particle(np.random.rand(3) * box_size, np.random.randn(3), radius=0.02) for _ in range(n_particles)]

    for particle in particles:
        box.add_particle(particle)

    # If plotting is enabled, set up the 3D plot
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim([0, box_size])
        ax.set_ylim([0, box_size])
        ax.set_zlim([0, box_size])

    # Simulation loop
    for step in range(n_steps):
        collisions = box.step(dt)  # Get collision info from box step

        # Plotting mode
        if plot:
            # Extract updated particle positions
            x_data = [particle.position[0] for particle in box.particles]
            y_data = [particle.position[1] for particle in box.particles]
            z_data = [particle.position[2] for particle in box.particles]

            # Set marker size based on the particle radius in the simulation
            sizes = [radius_to_marker_size(particle.radius) for particle in box.particles]
            colors = ['blue'] * n_particles

            # Check for interactions and adjust color/size on collisions
            for i in range(n_particles):
                if collisions[i]:  # If a collision happened
                    sizes[i] = radius_to_marker_size(box.particles[i].radius) * 2  # Increase size on collision
                    colors[i] = 'red'  # Change color on collision

            # Clear the plot and re-draw
            ax.clear()
            ax.set_xlim([0, box_size])
            ax.set_ylim([0, box_size])
            ax.set_zlim([0, box_size])

            # Plot the particles with updated sizes and colors
            ax.scatter(x_data, y_data, z_data, s=sizes, c=colors)
            plt.draw()
            plt.pause(0.01)  # Pause for a brief moment to visualize movement

    # If not plotting, output physical measurements
    if not plot:
        plt.figure()
        plt.plot(box.energy_over_time)
        plt.title("Average Energy Over Time")
        plt.xlabel("Time Step")
        plt.ylabel("Energy (J)")
        plt.show()

        plt.figure()
        plt.plot(box.temperature_over_time)
        plt.title("Temperature Over Time")
        plt.xlabel("Time Step")
        plt.ylabel("Temperature (K)")
        plt.show()

        plt.figure()
        plt.plot(box.pressure_over_time)
        plt.title("Pressure Over Time")
        plt.xlabel("Time Step")
        plt.ylabel("Pressure (Pa)")
        plt.show()

        plt.figure()
        plt.plot(box.entropy_over_time)
        plt.title("Entropy Over Time")
        plt.xlabel("Time Step")
        plt.ylabel("Entropy (J/K)")
        plt.show()

if __name__ == "__main__":
    # Define parameters
    N_PARTICLES = 100  # Number of particles
    BOX_SIZE = 10.0  # Size of the simulation box
    N_STEPS = 2000  # Number of simulation steps
    DT = 0.05  # Time step size

    # Run the simulation with plotting enabled
    run_simulation(N_STEPS, DT, BOX_SIZE, N_PARTICLES, plot=True)

    # Alternatively, run the simulation without plotting and output measurements
    # run_simulation(N_STEPS, DT, BOX_SIZE, N_PARTICLES, plot=False)
