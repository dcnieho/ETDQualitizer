classdef TestDataQualityMetrics < matlab.unittest.TestCase
    methods (Test)
        function testAccuracyBasic(testCase)
            x = [0; 1; -1];
            y = [0; 1; -1];
            target_x = 0;
            target_y = 0;
            [offset, offset_x, offset_y] = accuracy(x, y, target_x, target_y);
            testCase.verifyEqual(offset, 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(offset_x, 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(offset_y, 0, 'AbsTol', 1e-10);
        end

        function testAccuracyWithCustomCentralTendency(testCase)
            x = [0; 1; -1];
            y = [0; 1; -1];
            [offset, offset_x, offset_y] = accuracy(x, y, 0, 0, @median);
            testCase.verifyGreaterThanOrEqual(offset, 0);
            testCase.verifyGreaterThanOrEqual(offset_x, 0);
            testCase.verifyGreaterThanOrEqual(offset_y, 0);
        end

        function testStdFunction(testCase)
            x = [1; 2; 3];
            y = [4; 5; 6];
            [s, sx, sy] = std_(x, y);
            testCase.verifyEqual(sx, std(x, 1, 'omitnan'));
            testCase.verifyEqual(sy, std(y, 1, 'omitnan'));
            testCase.verifyEqual(s, hypot(sx, sy));
        end

        function testStdFunctionWithNan(testCase)
            x = [1; 2; nan; 3];
            y = [4; 5; nan; 6];
            [s, sx, sy] = std_(x, y);
            testCase.verifyEqual(sx, std(x, 1, 'omitnan'));
            testCase.verifyEqual(sy, std(y, 1, 'omitnan'));
            testCase.verifyEqual(s, hypot(sx, sy));
        end

        function testBCEA(testCase)
            x = randn(100000,1);
            y = randn(100000,1);
            [area, ~, ax1, ax2, aspect_ratio] = bcea(x, y, 0.68);
            testCase.verifyGreaterThan(area, 0);
            testCase.verifyEqual(aspect_ratio, 1, 'AbsTol', 1e-2);
            testCase.verifyEqual(area, 2*pi*ax1*ax2, 'AbsTol', 1e-4);
        end

        function testRmsS2S(testCase)
            x = [1; 2; 3; 4];
            y = [4; 5; 6; 7];
            [rms, rms_x, rms_y] = rms_s2s(x, y);
            testCase.verifyGreaterThanOrEqual(rms, 0);
            testCase.verifyEqual(rms_x, sqrt(mean(diff(x).^2)));
            testCase.verifyEqual(rms_y, sqrt(mean(diff(y).^2)));
        end

        function testDataLossFromInvalid(testCase)
            x = [1; NaN; 3];
            y = [4; 5; NaN];
            loss = data_loss_from_invalid(x, y);
            testCase.verifyEqual(loss, 2/3*100);
        end

        function testDataLossFromExpected(testCase)
            x = [1; NaN; 3];
            y = [4; 5; NaN];
            duration = 1;
            frequency = 3;
            loss = data_loss_from_expected(x, y, duration, frequency);
            testCase.verifyEqual(loss, (1 - 1/3)*100);
        end

        function testEffectiveFrequency(testCase)
            x = [1; NaN; 3];
            y = [4; 5; NaN];
            duration = 1;
            freq = effective_frequency(x, y, duration);
            testCase.verifyEqual(freq, 1);  % One valid sample
        end
    end
end
