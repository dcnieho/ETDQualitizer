import unittest
import numpy as np
from ETDQualitizer import accuracy, std, bcea, rms_s2s, data_loss, data_loss_from_expected, effective_frequency

class TestDataQualityMetrics(unittest.TestCase):
    def test_accuracy(self):
        x = np.array([0, 1, -1])
        y = np.array([0, 1, -1])
        offset, offset_x, offset_y = accuracy(x, y, 0, 0)
        self.assertAlmostEqual(offset, 0)
        self.assertAlmostEqual(offset_x, 0)
        self.assertAlmostEqual(offset_y, 0)

    def test_accuracy_custom_central_tendency(self):
        x = np.array([0, 1, -1])
        y = np.array([0, 1, -1])
        offset, offset_x, offset_y = accuracy(x, y, 0, 0, np.nanmedian)
        self.assertAlmostEqual(offset, 0)
        self.assertAlmostEqual(offset_x, 0)
        self.assertAlmostEqual(offset_y, 0)

    def test_std(self):
        x = np.array([1, 2, 3])
        y = np.array([4, 5, 6])
        s, sx, sy = std(x, y)
        self.assertAlmostEqual(sx, np.std(x))
        self.assertAlmostEqual(sy, np.std(y))
        self.assertAlmostEqual(s, np.hypot(sx, sy))

    def test_std_with_nan(self):
        x = np.array([1, 2, np.nan, 3])
        y = np.array([4, 5, np.nan, 6])
        s, sx, sy = std(x, y)
        self.assertAlmostEqual(sx, np.std(x))
        self.assertAlmostEqual(sy, np.std(y))
        self.assertAlmostEqual(s, np.hypot(sx, sy))

    def test_bcea(self):
        x = np.random.randn(100000)
        y = np.random.randn(100000)
        area, orientation, ax1, ax2, aspect_ratio = bcea(x, y)
        self.assertGreater(area, 0)
        self.assertAlmostEqual(aspect_ratio, 1, places=1)
        self.assertAlmostEqual(area, 2*np.pi*ax1*ax2, places=3)

    def test_rms_s2s(self):
        x = np.array([1, 2, 3])
        y = np.array([4, 5, 6])
        rms, rms_x, rms_y = rms_s2s(x, y)
        self.assertGreaterEqual(rms, 0)
        self.assertAlmostEqual(rms_x, np.sqrt(np.mean(np.diff(x)**2)))
        self.assertAlmostEqual(rms_y, np.sqrt(np.mean(np.diff(y)**2)))

    def test_data_loss(self):
        x = np.array([1, np.nan, 3])
        y = np.array([4, 5, np.nan])
        loss = data_loss(x, y)
        self.assertAlmostEqual(loss, 2/3*100)

    def test_data_loss_from_expected(self):
        x = np.array([1, np.nan, 3])
        y = np.array([4, 5, np.nan])
        loss = data_loss_from_expected(x, y, 1, 3)
        self.assertAlmostEqual(loss, (1 - 1/3)*100)

    def test_effective_frequency(self):
        x = np.array([1, np.nan, 3])
        y = np.array([4, 5, np.nan])
        freq = effective_frequency(x, y, 1)
        self.assertEqual(freq, 1)
