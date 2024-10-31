import numpy as np
import heapq  # For priority queue in Dijkstra's algorithm

class Lattice2D:
    def __init__(self, lattice_size):
        self.lattice_size = lattice_size  # The size of the NxN lattice
        self.lattice = np.empty((lattice_size, lattice_size), dtype=object)  # Initialize an empty lattice
        self.initialize_particles()  # Initialize each node with a label and random kinetic energy

    def initialize_particles(self):
        """Initialize particles at each node with a unique label and random kinetic energy."""
        for i in range(self.lattice_size):
            for j in range(self.lattice_size):
                node_label = f"Node({i},{j})"  # Create a unique node label
                kinetic_energy = np.random.rand()  # Assign random kinetic energy
                # Store the node as a dictionary with label and energy
                self.lattice[i, j] = {"label": node_label, "energy": kinetic_energy}
    def get_neighbors(self, i, j):
        """Get the neighboring nodes for a given node at (i, j)."""
        neighbors = []
        directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]  # Up, Down, Right, Left
        for di, dj in directions:
            ni, nj = i + di, j + dj
            # Ensure neighbors are within bounds of the lattice
            if 0 <= ni < self.lattice_size and 0 <= nj < self.lattice_size:
                neighbors.append(self.lattice[ni, nj])  # Append the neighboring node (label and energy)
        return neighbors

    def update_particle_energy(self, i, j, new_energy):
        """Update the kinetic energy of a given node (i, j)."""
        self.lattice[i, j]["energy"] = new_energy

    def get_total_energy(self):
        """Return the total kinetic energy of the lattice."""
        total_energy = 0
        for i in range(self.lattice_size):
            for j in range(self.lattice_size):
                total_energy += self.lattice[i, j]["energy"]
        return total_energy


    def dijkstra(self, start, target):
        """
        Use Dijkstra's algorithm to find the shortest path from start to target node.

        Parameters
        ----------
        start : tuple
            Starting coordinates (i, j).
        target : tuple
            Target coordinates (i, j).

        Returns
        -------
        float
            Shortest distance from start to target node.
        """
        # Priority queue to store (distance, node) and set for visited nodes
        pq = [(0, start)]
        distances = {start: 0}
        visited = set()

        while pq:
            current_distance, current_node = heapq.heappop(pq)
            if current_node in visited:
                continue

            visited.add(current_node)
            i, j = current_node

            # If we reached the target, return the distance
            if current_node == target:
                return current_distance

            # Check all neighbors
            for neighbor in self.get_neighbors(i, j):
                if neighbor in visited:
                    continue
                new_distance = current_distance + 1  # Distance between any two nodes is 1
                if neighbor not in distances or new_distance < distances[neighbor]:
                    distances[neighbor] = new_distance
                    heapq.heappush(pq, (new_distance, neighbor))

        return float('inf')  # Return infinity if no path found

    def measure_average_distance(self, particle_position, start_position):
        """
        Measure the average distance traveled by a particle using Dijkstra's algorithm.

        Parameters
        ----------
        particle_position : tuple
            Current position of the particle (i, j).
        start_position : tuple
            Starting position of the particle (i, j).

        Returns
        -------
        float
            Shortest distance from start to current position.
        """
        return self.dijkstra(start_position, particle_position)
