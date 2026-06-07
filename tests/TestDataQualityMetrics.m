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

        function testAccuracyMedianUsesFrechetMedian(testCase)
            x = [7.42; 73.96; 53.70; -84.53];
            y = [28.59; -37.31; 18.37; -25.95];
            grid_azi = linspace(-85, 85, 121);
            grid_ele = linspace(-40, 40, 101);
            [sx, sy, sz] = Fick_to_vector(x, y);
            sample_vectors = [sx, sy, sz];

            best_value = inf;
            best_target = [NaN, NaN];
            for target_azi = grid_azi
                for target_ele = grid_ele
                    [tx, ty, tz] = Fick_to_vector(target_azi, target_ele);
                    value = sum(acos(max(-1, min(1, sample_vectors * [tx; ty; tz]))));
                    if value < best_value
                        best_value = value;
                        best_target = [target_azi, target_ele];
                    end
                end
            end

            [offset, offset_x, offset_y] = accuracy(x, y, best_target(1), best_target(2), @median);
            testCase.verifyLessThan(offset, 1);
            testCase.verifyLessThan(abs(offset_x), 1);
            testCase.verifyLessThan(abs(offset_y), 1);

            legacy_vector = [median(sample_vectors(:,1)), median(sample_vectors(:,2)), median(sample_vectors(:,3))];
            legacy_vector = legacy_vector / norm(legacy_vector);
            legacy_value = sum(acos(max(-1, min(1, sample_vectors * legacy_vector.'))));
            testCase.verifyLessThan(best_value, legacy_value);
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
