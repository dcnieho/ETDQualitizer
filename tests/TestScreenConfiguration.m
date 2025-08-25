classdef TestScreenConfiguration < matlab.unittest.TestCase
    properties
        config
    end

    methods (TestMethodSetup)
        function createConfig(testCase)
            % Create a standard screen configuration for testing
            testCase.config = ScreenConfiguration(500, 300, 1920, 1080, 600);
        end
    end

    methods (Test)
        function testConstructor(testCase)
            % Verify properties are set correctly
            sc = ScreenConfiguration(400, 250, 1600, 900, 500);
            testCase.verifyEqual(sc.screen_size_x_mm, 400);
            testCase.verifyEqual(sc.screen_size_y_mm, 250);
            testCase.verifyEqual(sc.screen_res_x_pix, 1600);
            testCase.verifyEqual(sc.screen_res_y_pix, 900);
            testCase.verifyEqual(sc.viewing_distance_mm, 500);
        end

        function testPixToMm(testCase)
            [x_mm, y_mm] = testCase.config.pix_to_mm(960, 540);
            testCase.verifyEqual(x_mm, 250, 'AbsTol', 1e-10);
            testCase.verifyEqual(y_mm, 150, 'AbsTol', 1e-10);
        end

        function testMmToPix(testCase)
            [x_pix, y_pix] = testCase.config.mm_to_pix(250, 150);
            testCase.verifyEqual(x_pix, 960, 'AbsTol', 1e-10);
            testCase.verifyEqual(y_pix, 540, 'AbsTol', 1e-10);
        end

        function testMmToDeg(testCase)
            [azi, ele] = testCase.config.mm_to_deg(250, 0);
            testCase.verifyEqual(azi, atan2d(250, 600), 'AbsTol', 1e-10);
            testCase.verifyEqual(ele, 0, 'AbsTol', 1e-10);
        end

        function testDegToMm(testCase)
            azi = atan2d(250, 600);
            ele = 0;
            [x_mm, y_mm] = testCase.config.deg_to_mm(azi, ele);
            testCase.verifyEqual(x_mm, 250, 'AbsTol', 1e-10);
            testCase.verifyEqual(y_mm, 0, 'AbsTol', 1e-10);
        end

        function testPixToDeg(testCase)
            [azi, ele] = testCase.config.pix_to_deg(960, 540);
            expected_azi = atan2d(250, 600);
            expected_ele = atan2d(150, hypot(600, 250));
            testCase.verifyEqual(azi, expected_azi, 'AbsTol', 1e-10);
            testCase.verifyEqual(ele, expected_ele, 'AbsTol', 1e-10);
        end

        function testDegToPix(testCase)
            azi = atan2d(250, 600);
            ele = atan2d(150, hypot(600, 250));
            [x_pix, y_pix] = testCase.config.deg_to_pix(azi, ele);
            testCase.verifyEqual(x_pix, 960, 'AbsTol', 1e-10);
            testCase.verifyEqual(y_pix, 540, 'AbsTol', 1e-10);
        end

        function testScreenExtents(testCase)
            [x_deg, y_deg] = testCase.config.screen_extents();
            expected_x_deg = 2 * atan2d(250, 600);
            expected_y_deg = 2 * atan2d(150, 600);
            testCase.verifyEqual(x_deg, expected_x_deg, 'AbsTol', 1e-10);
            testCase.verifyEqual(y_deg, expected_y_deg, 'AbsTol', 1e-10);
        end
    end
end
