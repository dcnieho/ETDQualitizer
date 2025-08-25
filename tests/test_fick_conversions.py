import unittest
import math
from ETDQualitizer import Fick_to_vector, vector_to_Fick

class TestFickConversions(unittest.TestCase):
    def test_fick_to_vector_default_radius(self):
        x, y, z = Fick_to_vector(0, 0)
        self.assertAlmostEqual(x, 0)
        self.assertAlmostEqual(y, 0)
        self.assertAlmostEqual(z, 1)

    def test_fick_to_vector_custom_radius(self):
        x, y, z = Fick_to_vector(90, 0, 2)
        self.assertAlmostEqual(x, 2)
        self.assertAlmostEqual(y, 0)
        self.assertAlmostEqual(z, 0)

    def test_vector_to_fick(self):
        azi, ele = vector_to_Fick(0, 0, 1)
        self.assertAlmostEqual(azi, 0)
        self.assertAlmostEqual(ele, 0)

    def test_round_trip(self):
        azi, ele, r = 45, 30, 1.5
        x, y, z = Fick_to_vector(azi, ele, r)
        azi2, ele2 = vector_to_Fick(x, y, z)
        self.assertAlmostEqual(azi, azi2, places=5)
        self.assertAlmostEqual(ele, ele2, places=5)

    def test_vector_to_fick_negative_components(self):
        azi, ele = vector_to_Fick(-1, -1, -1)
        expected_azi = math.degrees(math.atan2(-1, -1))
        expected_ele = math.degrees(math.atan2(-1, math.hypot(-1, -1)))
        self.assertAlmostEqual(azi, expected_azi)
        self.assertAlmostEqual(ele, expected_ele)
