classdef TestDataQualityClass < matlab.unittest.TestCase
    properties
        duration
        freq
        timestamps
        azi
        ele
        screen
    end

    methods (TestMethodSetup)
        function setupData(testCase)
            testCase.duration = 1000;                   % 1000 s
            testCase.freq = 100;                        % 100 Hz
            testCase.timestamps = (0:1/testCase.freq:1000-1/testCase.freq)';
            n_samples = length(testCase.timestamps);
            testCase.azi = randn(n_samples,1);          % random azimuths
            testCase.ele = randn(n_samples,1);          % random elevations
            testCase.screen = ScreenConfiguration(500, 300, 1920, 1080, 600);
        end
    end

    methods (Test)
        function testConstructorWithDegrees(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            testCase.verifyEqual(dq.azi, testCase.azi);
            testCase.verifyEqual(dq.ele, testCase.ele);
            testCase.verifyEqual(dq.timestamps, testCase.timestamps);
        end

        function testConstructorWithPixels(testCase)
            x_pix = 960 + testCase.azi * 10;
            y_pix = 540 + testCase.ele * 10;
            dq = DataQuality(x_pix, y_pix, testCase.timestamps, 'pixels', testCase.screen);
            [azi_deg, ele_deg] = testCase.screen.pix_to_deg(x_pix, y_pix);
            testCase.verifyEqual(dq.azi, azi_deg, 'AbsTol', 1e-10);
            testCase.verifyEqual(dq.ele, ele_deg, 'AbsTol', 1e-10);
        end

        function testAccuracy(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            [offset, offset_x, offset_y] = dq.accuracy(0, 0);
            testCase.verifyEqual(offset, 0, 'AbsTol', 1e-2);
            testCase.verifyEqual(offset_x, 0, 'AbsTol', 1e-2);
            testCase.verifyEqual(offset_y, 0, 'AbsTol', 1e-2);
        end

        function testPrecisionRMS(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            [rms, rms_x, rms_y] = dq.precision_RMS_S2S();
            testCase.verifyGreaterThanOrEqual(rms, 0);
            testCase.verifyGreaterThanOrEqual(rms_x, 0);
            testCase.verifyGreaterThanOrEqual(rms_y, 0);
        end

        function testPrecisionSTD(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            [s, sx, sy] = dq.precision_STD();
            testCase.verifyGreaterThanOrEqual(s, 0);
            testCase.verifyEqual(s, hypot(sx, sy));
        end

        function testPrecisionBCEA(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            [area, orientation, ax1, ax2, aspect_ratio] = dq.precision_BCEA();
            testCase.verifyGreaterThan(area, 0);
            testCase.verifyGreaterThanOrEqual(aspect_ratio, 1);
        end

        function testPrecisionMovingWindow(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            s = dq.precision_using_moving_window(10, 'STD');
            testCase.verifyGreaterThanOrEqual(s, 0);
        end

        function testDataLoss(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            loss = dq.data_loss();
            testCase.verifyEqual(loss, 0);
        end

        function testDataLossFromExpected(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            loss = dq.data_loss_from_expected(testCase.freq);
            testCase.verifyEqual(loss, 0, 'AbsTol', 1e-10);
        end

        function testEffectiveFrequency(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            freq = dq.effective_frequency();
            testCase.verifyEqual(freq, testCase.freq, 'AbsTol', 1e-10);
        end

        function testGetDuration(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            duration = dq.get_duration();
            testCase.verifyEqual(duration, testCase.duration, 'AbsTol', 1e-10);
        end
    end
end
