import unittest
import numpy as np
from ETDQualitizer import accuracy, std, bcea, rms_s2s, data_loss_from_invalid, data_loss_from_expected, effective_frequency, Fick_to_vector

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

    def test_accuracy_nanmedian_uses_frechet_median(self):
        x = np.array([7.42, 73.96, 53.70, -84.53])
        y = np.array([28.59, -37.31, 18.37, -25.95])

        grid_azi = np.linspace(-85, 85, 121)
        grid_ele = np.linspace(-40, 40, 101)
        sample_vectors = np.column_stack(Fick_to_vector(x, y))

        def objective(target_azi, target_ele):
            target_vector = np.array(Fick_to_vector(target_azi, target_ele), dtype=float)
            dots = np.clip(sample_vectors @ target_vector, -1.0, 1.0)
            return float(np.sum(np.arccos(dots)))

        best_value = float('inf')
        best_target = (np.nan, np.nan)
        for target_azi in grid_azi:
            for target_ele in grid_ele:
                value = objective(target_azi, target_ele)
                if value < best_value:
                    best_value = value
                    best_target = (target_azi, target_ele)

        offset, offset_x, offset_y = accuracy(x, y, *best_target, np.nanmedian)
        self.assertLess(offset, 1.0)
        self.assertLess(abs(offset_x), 1.0)
        self.assertLess(abs(offset_y), 1.0)

        legacy_vector = np.array([
            np.nanmedian(sample_vectors[:, 0]),
            np.nanmedian(sample_vectors[:, 1]),
            np.nanmedian(sample_vectors[:, 2]),
        ], dtype=float)
        legacy_vector = legacy_vector / np.linalg.norm(legacy_vector)
        legacy_value = float(np.sum(np.arccos(np.clip(sample_vectors @ legacy_vector, -1.0, 1.0))))
        self.assertLess(best_value, legacy_value)

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
        self.assertAlmostEqual(sx, np.nanstd(x))
        self.assertAlmostEqual(sy, np.nanstd(y))
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

    def test_data_loss_from_invalid(self):
        x = np.array([1, np.nan, 3])
        y = np.array([4, 5, np.nan])
        loss = data_loss_from_invalid(x, y)
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
