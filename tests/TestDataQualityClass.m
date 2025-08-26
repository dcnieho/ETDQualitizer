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
            testCase.timestamps = (0:1/testCase.freq:testCase.duration-1/testCase.freq)';
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

        function testPrecisionRMSLotsOfData(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            [rms, rms_x, rms_y] = dq.precision_RMS_S2S();
            testCase.verifyEqual(rms, 2, 'AbsTol', 1e-2);
            testCase.verifyEqual(rms_x, sqrt(2), 'AbsTol', 1e-2);
            testCase.verifyEqual(rms_y, sqrt(2), 'AbsTol', 1e-2);
        end

        function testPrecisionRMSOneAxis(testCase)
            az = [1 2 4];
            el = [1 1 1];
            dq = DataQuality(az, el, [0 1 2], 'degrees');
            [rms, rms_x, rms_y] = dq.precision_RMS_S2S();
            expected_rms = sqrt(mean([1 2].^2));
            testCase.verifyEqual(rms, expected_rms, 'AbsTol', 1e-10);
            testCase.verifyEqual(rms_x, expected_rms, 'AbsTol', 1e-10);
            testCase.verifyEqual(rms_y, 0);
        end

        function testPrecisionRMSTwoAxes(testCase)
            az = [1 2 4];
            el = [1 2 4];
            dq = DataQuality(az, el, [0 1 2], 'degrees');
            [rms, rms_x, rms_y] = dq.precision_RMS_S2S();
            expected_rms_xy = sqrt(mean([1 2].^2));
            expected_rms = hypot(expected_rms_xy,expected_rms_xy);
            testCase.verifyEqual(rms, expected_rms, 'AbsTol', 1e-10);
            testCase.verifyEqual(rms_x, expected_rms_xy, 'AbsTol', 1e-10);
            testCase.verifyEqual(rms_y, expected_rms_xy, 'AbsTol', 1e-10);
        end

        function testPrecisionSTD(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            [s, sx, sy] = dq.precision_STD();
            testCase.verifyEqual(sx, 1, 'AbsTol', 1e-2);
            testCase.verifyEqual(sy, 1, 'AbsTol', 1e-2);
            testCase.verifyEqual(s, sqrt(2), 'AbsTol', 1e-2);
        end

        function testPrecisionBCEA(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            [area, ~, ax1, ax2, aspect_ratio] = dq.precision_BCEA();
            testCase.verifyGreaterThan(area, 0);
            testCase.verifyEqual(aspect_ratio, 1, 'AbsTol', 1e-2);
            testCase.verifyEqual(area, 2*pi*ax1*ax2, 'AbsTol', 1e-4);
        end

        function testPrecisionMovingWindow(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            s = dq.precision_using_moving_window(50, 'STD');
            testCase.verifyLessThanOrEqual(s, sqrt(2)); % std of whole sequence is about sqrt(2), when taking it in smaller windows it'll come out a bit smaller even though generation process is stationary
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
            fs = dq.effective_frequency();
            testCase.verifyEqual(fs, testCase.freq, 'AbsTol', 1e-10);
        end

        function testGetDuration(testCase)
            dq = DataQuality(testCase.azi, testCase.ele, testCase.timestamps, 'degrees');
            dur = dq.get_duration();
            testCase.verifyEqual(dur, testCase.duration, 'AbsTol', 1e-10);
        end
    end
end
