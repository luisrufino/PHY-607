import numpy as np
import heapq  # For priority queue in Dijkstra's algorithm

class Lattice2D:
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
