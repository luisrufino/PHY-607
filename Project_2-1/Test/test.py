import unittest
from Monte_sim.box import Box
from Monte_sim.particle import Particle


class TestSimulation(unittest.TestCase):

    def test_particle_wall_collision(self):
        box = Box(size=10)
        particle = Particle(position=[9.5, 5, 5], velocity=[1, 0, 0], radius=0.5)
        box.add_particle(particle)
        box.handle_wall_collisions(particle)
        self.assertEqual(particle.velocity[0], -1)


if __name__ == '__main__':
    unittest.main()
