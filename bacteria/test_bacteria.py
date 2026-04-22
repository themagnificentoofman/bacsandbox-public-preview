import math
import unittest
import sys
import os

if __package__ is None and not hasattr(sys, 'frozen'):
    path = os.path.realpath(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(os.path.dirname(path)))

from bacteria import Bacterium, BacteriaSpecies, FoodNode, Simulation, FPS

class MockSpecies:
    def __init__(self, **kwargs):
        self.food_capacity = kwargs.get('food_capacity', 10.0)
        self.metabolic_rate = kwargs.get('metabolic_rate', 0.1)
        self.consumption_speed = kwargs.get('consumption_speed', 1.0)
        self.speed = kwargs.get('speed', 40.0)
        self.sensor_radius = kwargs.get('sensor_radius', 200.0)
        self.greediness = kwargs.get('greediness', 0.5)
        self.aggression = kwargs.get('aggression', 0.5)
        self.defense = kwargs.get('defense', 0.5)
        self.rep_thresh = kwargs.get('rep_thresh', 0.8)
        self.share_thresh = kwargs.get('share_thresh', 0.6)
        self.predation_thresh = kwargs.get('predation_thresh', 0.9)
        self.base_size = kwargs.get('base_size', 1.0)

class TestBacteriumMath(unittest.TestCase):
    def test_bacterium_scaling_laws(self):
        # Base size 1.0
        sp_base = MockSpecies(food_capacity=10.0, metabolic_rate=0.1, consumption_speed=1.0, speed=40.0)
        b_base = Bacterium(sp_base, 0, 0, size=1.0)
        
        self.assertAlmostEqual(b_base.food_capacity, 10.0)
        self.assertAlmostEqual(b_base.metabolic_rate, 0.1)
        self.assertAlmostEqual(b_base.consumption_speed, 1.0)
        self.assertAlmostEqual(b_base.speed, 40.0)
        self.assertAlmostEqual(b_base.turn_rate, math.pi)

        # Base size 2.0
        b_large = Bacterium(sp_base, 0, 0, size=2.0)
        # Volumetric scaling: 10 * (2^3) = 80
        self.assertAlmostEqual(b_large.food_capacity, 80.0)
        # Kleiber's Law: 0.1 * (2^2.25)
        self.assertAlmostEqual(b_large.metabolic_rate, 0.1 * (2 ** 2.25))
        # Surface area: 1.0 * (2^2) = 4.0
        self.assertAlmostEqual(b_large.consumption_speed, 4.0)
        # Kinematic scaling: 40 / 2 = 20
        self.assertAlmostEqual(b_large.speed, 20.0)
        self.assertAlmostEqual(b_large.turn_rate, math.pi / 2.0)

    def test_metabolism_decay(self):
        sp = MockSpecies(food_capacity=10.0, metabolic_rate=1.44) # 1.44 / 144 FPS = 0.01 per step
        b = Bacterium(sp, 0, 0, size=1.0)
        b.age = 0 # reset age for exactness
        initial_food = b.food_level
        b.step([], -100, 100)
        
        # calculate expected reduction with age factor
        # age factor = 1.0 + (age / (180.0 * lifespan_variance))
        # for age=0, factor=1.0
        # Wait, age is incremented in step before food reduction
        expected_food = initial_food - (1.44 / FPS) * (1.0 + ((0.0 + 1.0/FPS) / (180.0 * b.lifespan_variance)))
        self.assertAlmostEqual(b.food_level, expected_food)

    def test_chemotaxis_score(self):
        # A food node at distance 100 with amount 50
        sp = MockSpecies(sensor_radius=200.0)
        b = Bacterium(sp, 0, 0, size=1.0)
        fnode = FoodNode(100, 0, 50)
        # The node should be visible
        b.step([fnode], -1000, 1000)
        # After step, if target is found, it will turn towards it.
        # It's hard to test the pure score without intercepting, but making sure no crash
        

if __name__ == '__main__':
    unittest.main()
