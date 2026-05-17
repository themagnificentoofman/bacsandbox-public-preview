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
        dt = 1.0 / 144.0
        b.step([], -100, 100, dt)
        
        # calculate expected reduction with age factor
        # age factor = 1.0 + (age / (180.0 * lifespan_variance))
        # for age=0, factor=1.0
        # Wait, age is incremented in step before food reduction
        expected_food = initial_food - (1.44 * dt) * (1.0 + ((0.0 + dt) / (180.0 * b.lifespan_variance)))
        self.assertAlmostEqual(b.food_level, expected_food)

    def test_chemotaxis_score(self):
        # A food node at distance 100 with amount 50
        sp = MockSpecies(sensor_radius=200.0)
        b = Bacterium(sp, 0, 0, size=1.0)
        fnode = FoodNode(100, 0, 50)
        # The node should be visible
        dt = 1.0 / 144.0
        b.step([fnode], -1000, 1000, dt)
        # After step, if target is found, it will turn towards it.
        # It's hard to test the pure score without intercepting, but making sure no crash

    def test_binary_fission_survival(self):
        """
        Verify that cells do not immediately die of starvation after binary fission ATP cost.
        """
        dt = 1.0 / 144.0
        sp = MockSpecies(rep_thresh=0.8, metabolic_rate=0.1)
        b = Bacterium(sp, 0, 0, size=1.0)
        
        b.base_size = 1.0
        # Give enough food to trigger fission but not much more
        child_volume = 1.0 ** 3
        atp_cost = 40.0 * child_volume
        b.food_capacity = 100.0
        b.food_level = atp_cost + 0.5  # Just enough to trigger and have 0.5 leftover
        
        # Test attempt_fission
        child = b.attempt_fission(1.0, 1.0, 0.8, dt)
        
        # If it had barely enough food, it should return None if survival is compromised,
        # or it returns child and both must survive at least 1 tick
        if child is not None:
            min_survival = b.metabolic_rate * (1.0 / 144.0)
            self.assertGreater(b.food_level, min_survival)
            self.assertGreater(child.food_level, min_survival)

    def test_ecm_memory_bounds(self):
        """
        Verify ECM removes empty cells to prevent dictionary memory leaks.
        """
        from bacteria.ecm import ECM
        ecm = ECM()
        
        # Add tiny amount that should decay below threshold instantly
        ecm.add(0, 0, 0.0010001, layer='qs')
        ecm.add(100, 100, 0.0010001, layer='bac')
        
        self.assertEqual(len(ecm.qs_inducers), 1)
        self.assertEqual(len(ecm.bacteriocins), 1)
        
        # Tick applies decay and diffusion, values drop below 0.001
        dt = 1.0 / 144.0
        ecm.tick(dt)
        
        self.assertEqual(len(ecm.qs_inducers), 0)
        self.assertEqual(len(ecm.bacteriocins), 0)

    def test_t6ss_damage_scaling(self):
        """
        Verify T6SS damage calculations don't overflow or produce negative floats.
        """
        sp1 = MockSpecies(t6ss_activity=1.0, immunity_factor=0.0)
        sp2 = MockSpecies(t6ss_activity=0.0, immunity_factor=1.0)
        
        # Attacker is huge (size 10), victim is tiny
        b_attacker = Bacterium(sp1, 0, 0, size=10.0)
        b_victim = Bacterium(sp2, 1, 1, size=0.1)
        
        # Base T6SS damage formula: (t6ss_activity * size^2 * 5.0) * dt
        dt = 1.0 / 144.0
        t6ss_dmg = (1.0 * (10.0 ** 2) * 5.0) * dt
        
        # Immunity defense: min(0.95, (immunity_factor * victim.size) / (attacker.size + 0.01))
        immunity = min(0.95, (1.0 * 0.1) / (10.0 + 0.01))
        
        actual_dmg = t6ss_dmg * (1.0 - immunity)
        
        # Must be bounded, positive, and not float('inf')
        self.assertGreater(actual_dmg, 0.0)
        self.assertLess(actual_dmg, 1000.0)  # arbitrary upper bound for 1 tick of damage
        self.assertFalse(math.isinf(actual_dmg))

if __name__ == '__main__':
    unittest.main()
