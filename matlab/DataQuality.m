classdef DataQuality
    %DATAQUALITY Class for calculating Data Quality from a gaze data segment
    %
    %   Provides methods for assessing the quality of gaze data, including accuracy,
    %   precision, data loss, and effective sampling frequency.
    %
    %   Notes:
    %   - Missing data should be coded as NaN, not as special values like (0,0) or (-xres,-yres).
    %   - Missing samples should not be removed, or RMS calculations will be incorrect.
    %   - Timestamps should be in seconds.
    %   - All angular positions are expected to be expressed in Fick angles.
    %
    %   Example:
    %       sc = ScreenConfiguration(500, 300, 1920, 1080, 600);
    %       dq = DataQuality([0;1;-1], [0;1;-1], [0;1;2], 'pixels', sc);
    %       dq.accuracy(0, 0);
    %       dq.precision_RMS_S2S();
    %       dq.data_loss();
    properties (SetAccess=private)
        timestamps      % Vector of timestamps in seconds
        azi             % Azimuth angles in degrees (Fick angles)
        ele             % Elevation angles in degrees (Fick angles)
    end

    methods
        function obj = DataQuality(x, y, timestamps, unit, screen)
            % DATAQUALITY Construct a DataQuality object from gaze data and timestamps.
            %
            % Inputs:
            %   x, y        - Horizontal and vertical gaze positions (pixels or degrees)
            %   timestamps  - Vector of timestamps in seconds
            %   unit        - Unit of gaze data: 'pixels' or 'degrees'
            %   screen      - Optional ScreenConfiguration object (required if unit is 'pixels')
            arguments
                x           (:,1) {mustBeNumeric}
                y           (:,1) {mustBeNumeric}
                timestamps  (:,1) {mustBeNumeric}
                unit        (1,:) char {mustBeMember(unit,{'pixels','degrees'})}
                screen      {mustBeScalarOrEmpty,mustBeA(screen,{'ScreenConfiguration','double'})} = []
            end
            obj.timestamps = timestamps;
            if strcmp(unit, 'pixels')
                if isempty(screen)
                    error('If unit is "pixels", a screen configuration must be supplied')
                end
                [x,y] = screen.pix_to_deg(x, y);
            end
            obj.azi = x;
            obj.ele = y;
        end


        function [offset, offset_azi, offset_ele] = accuracy(obj, target_azi, target_ele, central_tendency_fun)
            % ACCURACY Calculates the accuracy of gaze data relative to a known target location.
            %
            % Inputs:
            %   target_azi - Target azimuth in degrees
            %   target_ele - Target elevation in degrees
            %   central_tendency_fun - Function to compute central tendency (default: mean)
            %
            % Outputs:
            %   offset     - Total angular offset in degrees
            %   offset_azi - Horizontal offset in degrees
            %   offset_ele - Vertical offset in degrees
            arguments
                obj
                target_azi          (1,1) {mustBeNumeric}
                target_ele          (1,1) {mustBeNumeric}
                central_tendency_fun(1,1) {mustBeA(central_tendency_fun,'function_handle')} = @(x) mean(x,'omitnan')
            end
            [offset, offset_azi, offset_ele] = accuracy(obj.azi, obj.ele, target_azi, target_ele, central_tendency_fun);
        end

        function [rms, rms_azi, rms_ele] = precision_RMS_S2S(obj, central_tendency_fun)
            % PRECISION_RMS_S2S Calculates precision as root mean square of sample-to-sample distances.
            %
            % Input:
            %   central_tendency_fun - Function to compute central tendency (default: mean)
            %
            % Outputs:
            %   rms     - Total RMS precision in degrees
            %   rms_azi - RMS of azimuthal component
            %   rms_ele - RMS of elevation component
            arguments
                obj
                central_tendency_fun(1,1) {mustBeA(central_tendency_fun,'function_handle')} = @(x) mean(x,'omitnan')
            end
            [rms, rms_azi, rms_ele] = rms_s2s(obj.azi, obj.ele, central_tendency_fun);
        end

        function [std__, std_azi, std_ele] = precision_STD(obj)
            % PRECISION_STD Calculates precision as standard deviation of gaze positions.
            %
            % Outputs:
            %   std__   - Total standard deviation in degrees
            %   std_azi - Standard deviation of azimuth
            %   std_ele - Standard deviation of elevation
            [std__, std_azi, std_ele] = std_(obj.azi, obj.ele);
        end

        function [area, orientation, ax1, ax2, aspect_ratio] = precision_BCEA(obj, P)
            % PRECISION_BCEA Calculates the Bivariate Contour Ellipse Area (BCEA) and ellipse parameters.
            %
            % Input:
            %   P - Proportion of data to include in the ellipse (default: 0.68)
            %
            % Outputs:
            %   area         - BCEA value
            %   orientation  - Orientation of the ellipse in degrees
            %   ax1          - Length of major axis
            %   ax2          - Length of minor axis
            %   aspect_ratio - Ratio of major to minor axis
            arguments
                obj
                P   (1,1) {mustBeNumeric} = 0.68    % for BCEA: cumulative probability of area under the multivariate normal
            end
            [area, orientation, ax1, ax2, aspect_ratio] = bcea(obj.azi, obj.ele, P);
        end

        function precision = precision_using_moving_window(obj, window_length, metric, input_args, aggregation_fun)
            % PRECISION_USING_MOVING_WINDOW Calculates precision using a moving window approach.
            %
            % Inputs:
            %   window_length   - Length of the moving window in number of samples
            %   metric          - Precision metric: 'RMS-S2S', 'STD', or 'BCEA'
            %   input_args      - Additional arguments passed to the metric function
            %   aggregation_fun - Function to aggregate windowed precision values (default: median)
            %
            % Output:
            %   precision - Aggregated precision value
            arguments
                obj
                window_length   (1,1) {mustBeInteger}
                metric          (1,:) char {mustBeMember(metric,{'RMS-S2S','STD','BCEA'})}
                input_args      (1,:) cell = {}
                aggregation_fun (1,1) {mustBeA(aggregation_fun,'function_handle')} = @(x) median(x,'omitnan')
            end

            precision = precision_using_moving_window(obj.azi, obj.ele, window_length, metric, input_args, aggregation_fun);
        end

        function loss_percentage = data_loss_from_invalid(obj)
            % DATA_LOSS_FROM_INVALID Calculates the proportion of missing data (coded as NaN).
            %
            % Output:
            %   loss_percentage - Percentage of missing samples
            loss_percentage = data_loss_from_invalid(obj.azi, obj.ele);
        end

        function loss_percentage = data_loss_from_expected(obj, frequency)
            % DATAS_LOSS_FROM_EXPECTED Estimates data loss based on expected number of samples.
            %
            % Input:
            %   frequency - Expected sampling frequency in Hz
            %
            % Output:
            %   loss_percentage - Percentage of missing samples
            loss_percentage = data_loss_from_expected(obj.azi, obj.ele, obj.get_duration(), frequency);
        end

        function freq = effective_frequency(obj)
            % EFFECTIVE_FREQUENCY Calculates the effective sampling frequency based on timestamps.
            %
            % Output:
            %   freq - Effective frequency in Hz (rate of valid samples)
            freq = effective_frequency(obj.azi, obj.ele, obj.get_duration());
        end

        function duration = get_duration(obj)
            % GET_DURATION Computes the total duration of the gaze recording, including the last sample.
            %
            % Output:
            %   duration - Duration in seconds

            % to get duration right, we need to include duration of last
            % sample
            isi     = median(diff(obj.timestamps));
            duration= obj.timestamps(end)-obj.timestamps(1)+isi;
        end
    end
end
