classdef TestFickConversions < matlab.unittest.TestCase
    methods (Test)
        function testFickToVectorDefaultRadius(testCase)
            azi = 0;
            ele = 0;
            [x, y, z] = Fick_to_vector(azi, ele);
            testCase.verifyEqual(x, 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(y, 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(z, 1, 'AbsTol', 1e-10);
        end

        function testFickToVectorCustomRadius(testCase)
            azi = 90;
            ele = 0;
            r = 2;
            [x, y, z] = Fick_to_vector(azi, ele, r);
            testCase.verifyEqual(x, 2, 'AbsTol', 1e-10);
            testCase.verifyEqual(y, 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(z, 0, 'AbsTol', 1e-10);
        end

        function testVectorToFickZeroVector(testCase)
            x = 0;
            y = 0;
            z = 1;
            [azi, ele] = vector_to_Fick(x, y, z);
            testCase.verifyEqual(azi, 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(ele, 0, 'AbsTol', 1e-10);
        end

        function testRoundTripConversion(testCase)
            azi = 45;
            ele = 30;
            r = 1.5;
            [x, y, z] = Fick_to_vector(azi, ele, r);
            [azi2, ele2] = vector_to_Fick(x, y, z);
            testCase.verifyEqual(azi2, azi, 'AbsTol', 1e-10);
            testCase.verifyEqual(ele2, ele, 'AbsTol', 1e-10);
        end

        function testVectorToFickNegativeComponents(testCase)
            x = -1;
            y = -1;
            z = -1;
            [azi, ele] = vector_to_Fick(x, y, z);
            expected_azi = atan2d(-1, -1);
            expected_ele = atan2d(-1, hypot(-1, -1));
            testCase.verifyEqual(azi, expected_azi, 'AbsTol', 1e-10);
            testCase.verifyEqual(ele, expected_ele, 'AbsTol', 1e-10);
        end
    end
end
