import unittest
import numpy as np
from ETDQualitizer import DataQuality, ScreenConfiguration

class TestDataQualityClass(unittest.TestCase):
    def setUp(self):
        self.duration = 1000        # 1000 s
        self.freq = 100             # 100 Hz
        self.timestamps = np.arange(0, self.duration, 1./self.freq)
        self.n_samples = len(self.timestamps)
        self.azi = np.random.randn(self.n_samples)
        self.ele = np.random.randn(self.n_samples)
        self.screen = ScreenConfiguration(500, 300, 1920, 1080, 600)

    def test_constructor_degrees(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        np.testing.assert_array_equal(dq.azi, self.azi)
        np.testing.assert_array_equal(dq.ele, self.ele)
        np.testing.assert_array_equal(dq.timestamps, self.timestamps)

    def test_constructor_pixels(self):
        x_pix = 960 + self.azi * 10
        y_pix = 540 + self.ele * 10
        dq = DataQuality(x_pix, y_pix, self.timestamps, 'pixels', self.screen)
        azi_deg, ele_deg = self.screen.pix_to_deg(x_pix, y_pix)
        np.testing.assert_array_almost_equal(dq.azi, azi_deg)
        np.testing.assert_array_almost_equal(dq.ele, ele_deg)

    def test_accuracy(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        offset, offset_x, offset_y = dq.accuracy(0, 0)
        self.assertAlmostEqual(offset, 0, places=1)
        self.assertAlmostEqual(offset_x, 0, places=1)
        self.assertAlmostEqual(offset_y, 0, places=1)

    def test_precision_rms_lots_of_data(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        rms, rms_x, rms_y = dq.precision_RMS_S2S()
        self.assertAlmostEqual(rms, 2, places=1)
        self.assertAlmostEqual(rms_x, np.sqrt(2), places=1)
        self.assertAlmostEqual(rms_y, np.sqrt(2), places=1)

    def test_precision_rms_one_axis(self):
        azi = np.array([1., 2., 4.])
        ele = np.array([1., 1., 1.])
        dq = DataQuality(azi, ele, np.array([0., 1., 2.]), 'degrees')
        rms, rms_x, rms_y = dq.precision_RMS_S2S()
        expected_rms = np.sqrt(np.mean(np.array([1., 2.])**2))
        self.assertGreaterEqual(rms, expected_rms)
        self.assertGreaterEqual(rms_x, expected_rms)
        self.assertGreaterEqual(rms_y, 0)

    def test_precision_rms_two_axes(self):
        azi = np.array([1., 2., 4.])
        ele = np.array([1., 2., 4.])
        dq = DataQuality(azi, ele, np.array([0., 1., 2.]), 'degrees')
        rms, rms_x, rms_y = dq.precision_RMS_S2S()
        expected_rms_xy = np.sqrt(np.mean(np.array([1., 2.])**2))
        expected_rms = np.hypot(expected_rms_xy,expected_rms_xy)
        self.assertGreaterEqual(rms, expected_rms)
        self.assertGreaterEqual(rms_x, expected_rms_xy)
        self.assertGreaterEqual(rms_y, expected_rms_xy)

    def test_precision_std(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        s, sx, sy = dq.precision_STD()
        self.assertAlmostEqual(s, np.sqrt(2), places=1)
        self.assertAlmostEqual(sx, 1, places=1)
        self.assertAlmostEqual(sy, 1, places=1)

    def test_precision_bcea(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        area, _, ax1, ax2, aspect_ratio = dq.precision_BCEA()
        self.assertGreater(area, 0)
        self.assertAlmostEqual(aspect_ratio, 1, places=1)
        self.assertAlmostEqual(area, 2*np.pi*ax1*ax2, places=3)

    def test_precision_std_moving_window(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        s = dq.precision_using_moving_window(50, 'STD')
        self.assertLessEqual(s, np.sqrt(2))     # std of whole sequence is about sqrt(2), when taking it in smaller windows it'll come out a bit smaller even though generation process is stationary

    def test_data_loss_from_invalid(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        self.assertEqual(dq.data_loss_from_invalid(), 0)

    def test_data_loss_from_expected(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        loss = dq.data_loss_from_expected(self.freq)
        self.assertAlmostEqual(loss, 0)

    def test_effective_frequency(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        freq = dq.effective_frequency()
        self.assertEqual(freq, self.freq)

    def test_get_duration(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        duration = dq.get_duration()
        self.assertAlmostEqual(duration, self.duration)
