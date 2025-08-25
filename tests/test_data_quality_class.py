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
        self.assertAlmostEqual(offset, 0, places=2)
        self.assertAlmostEqual(offset_x, 0, places=2)
        self.assertAlmostEqual(offset_y, 0, places=2)

    def test_precision_rms(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        rms, rms_x, rms_y = dq.precision_RMS_S2S()
        self.assertGreaterEqual(rms, 0)
        self.assertGreaterEqual(rms_x, 0)
        self.assertGreaterEqual(rms_y, 0)

    def test_precision_std(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        s, sx, sy = dq.precision_STD()
        self.assertGreaterEqual(s, 0)
        self.assertAlmostEqual(s, np.hypot(sx, sy))

    def test_precision_bcea(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        area, _, _, _, aspect_ratio = dq.precision_BCEA()
        self.assertGreater(area, 0)
        self.assertGreaterEqual(aspect_ratio, 1)

    def test_precision_std_moving_window(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        s = dq.precision_using_moving_window(10, 'STD')
        self.assertGreater(s, 0)

    def test_data_loss(self):
        dq = DataQuality(self.azi, self.ele, self.timestamps, 'degrees')
        self.assertEqual(dq.data_loss(), 0)

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
