import unittest
import math
from ETDQualitizer import ScreenConfiguration

class TestScreenConfiguration(unittest.TestCase):
    def setUp(self):
        self.config = ScreenConfiguration(500, 300, 1920, 1080, 600)

    def test_constructor(self):
        sc = ScreenConfiguration(400, 250, 1600, 900, 500)
        self.assertEqual(sc.screen_size_x_mm, 400)
        self.assertEqual(sc.screen_size_y_mm, 250)
        self.assertEqual(sc.screen_res_x_pix, 1600)
        self.assertEqual(sc.screen_res_y_pix, 900)
        self.assertEqual(sc.viewing_distance_mm, 500)

    def test_pix_to_mm(self):
        x_mm, y_mm = self.config.pix_to_mm(960, 540)
        self.assertAlmostEqual(x_mm, 250)
        self.assertAlmostEqual(y_mm, 150)

    def test_mm_to_pix(self):
        x_pix, y_pix = self.config.mm_to_pix(250, 150)
        self.assertAlmostEqual(x_pix, 960)
        self.assertAlmostEqual(y_pix, 540)

    def test_mm_to_deg(self):
        azi, ele = self.config.mm_to_deg(250, 0)
        self.assertAlmostEqual(azi, math.degrees(math.atan2(250, 600)))
        self.assertAlmostEqual(ele, 0)

    def test_deg_to_mm(self):
        azi = math.degrees(math.atan2(250, 600))
        ele = 0
        x_mm, y_mm = self.config.deg_to_mm(azi, ele)
        self.assertAlmostEqual(x_mm, 250)
        self.assertAlmostEqual(y_mm, 0)

    def test_pix_to_deg(self):
        azi, ele = self.config.pix_to_deg(960, 540)
        self.assertIsInstance(azi, float)
        self.assertIsInstance(ele, float)

    def test_deg_to_pix(self):
        azi = math.degrees(math.atan2(250, 600))
        ele = math.degrees(math.atan2(150, math.hypot(600, 250)))
        x_pix, y_pix = self.config.deg_to_pix(azi, ele)
        self.assertAlmostEqual(x_pix, 960, places=1)
        self.assertAlmostEqual(y_pix, 540, places=1)

    def test_pix_to_mm_and_mm_to_pix_consistency(self):
        x, y = self.config.pix_to_mm(960, 540)
        xp, yp = self.config.mm_to_pix(x, y)
        self.assertAlmostEqual(xp, 960, places=10)
        self.assertAlmostEqual(yp, 540, places=10)

    def test_mm_to_deg_and_deg_to_mm_consistency(self):
        azi, ele = self.config.mm_to_deg(250, 0)
        x, y = self.config.deg_to_mm(azi, ele)
        self.assertAlmostEqual(x, 250, places=10)
        self.assertAlmostEqual(y, 0, places=10)

    def test_screen_extents(self):
        x_deg, y_deg = self.config.screen_extents()
        expected_x_deg = math.degrees(2 * math.atan2(250, 600))
        expected_y_deg = math.degrees(2 * math.atan2(150, 600))
        self.assertAlmostEqual(x_deg, expected_x_deg)
        self.assertAlmostEqual(y_deg, expected_y_deg)
