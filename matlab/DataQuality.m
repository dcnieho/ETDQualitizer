classdef DataQuality
    % N.B: for this module it is assumed that any missing data are not coded with some special value
    % such as (0,0) or (-xres,-yres) but as nan. Missing data should also not be removed, or the RMS
    % calculation would be incorrect.
    %
    % timestamps should be in seconds.
    %
    % all angular positions are expected to be expressed in Fick angles.
    properties (SetAccess=private)
        timestamps
        x
        y
    end

    methods
        function obj = DataQuality(x, y, timestamps, unit, screen)
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
            obj.x = x;
            obj.y = y;
        end


        function [offset, offset_x, offset_y] = accuracy(obj, target_x_deg, target_y_deg, central_tendency_fun)
            arguments
                obj
                target_x_deg        (1,1) {mustBeNumeric}
                target_y_deg        (1,1) {mustBeNumeric}
                central_tendency_fun(1,1) {mustBeA(central_tendency_fun,'function_handle')} = @(x) mean(x,'omitnan')
            end
            [offset, offset_x, offset_y] = accuracy(obj.x, obj.y, target_x_deg, target_y_deg, central_tendency_fun);
        end

        function [rms, rms_x, rms_y] = precision_RMS_S2S(obj, central_tendency_fun)
            arguments
                obj
                central_tendency_fun(1,1) {mustBeA(central_tendency_fun,'function_handle')} = @(x) mean(x,'omitnan')
            end
            [rms, rms_x, rms_y] = rms_s2s(obj.x, obj.y, central_tendency_fun);
        end

        function [std__, std_x, std_y] = precision_STD(obj)
            [std__, std_x, std_y] = std_(obj.x, obj.y);
        end

        function [area, orientation, ax1, ax2, aspect_ratio] = precision_BCEA(obj, P)
            arguments
                obj
                P   (1,1) {mustBeNumeric} = 0.68    % for BCEA: cumulative probability of area under the multivariate normal
            end
            [area, orientation, ax1, ax2, aspect_ratio] = bcea(obj.x, obj.y, P);
        end

        function precision = precision_using_moving_window(obj, window_length, metric, input_args, aggregation_fun)
            arguments
                obj
                window_length   (1,1) {mustBeInteger}
                metric          (1,:) char {mustBeMember(metric,{'RMS_S2S','STD','BCEA'})}
                input_args      (1,:) cell = {}
                aggregation_fun (1,1) {mustBeA(aggregation_fun,'function_handle')} = @(x) median(x,'omitnan')
            end

            switch metric
                case 'RMS_S2S'
                    fun =  @rms_s2s;
                case 'STD'
                    fun =  @std_;
                case 'BCEA'
                    fun =  @bcea;
            end
            
            % get number of samples in data
            ns  = length(obj.x);

            if window_length < ns % if number of samples in data exceeds window size
                values = nan(1,ns-window_length); % pre-allocate
                for p=1:ns-window_length+1
                    values(p) = fun(obj.x(p:p+window_length-1), obj.y(p:p+window_length-1), input_args{:});
                end
                precision = aggregation_fun(values);
            else % if too few samples in data
                precision = NaN;
            end
        end

        function loss_percentage = data_loss(obj)
            loss_percentage = data_loss(obj.x, obj.y);
        end

        function loss_percentage = data_loss_nominal(obj, frequency)
            loss_percentage = data_loss_nominal(obj.x, obj.y, obj.get_duration(), frequency);
        end

        function freq = effective_frequency(obj)
            % rate of valid samples
            freq = effective_frequency(obj.x, obj.y, obj.get_duration());
        end

        function duration = get_duration(obj)
            % to get duration right, we need to include duration of last
            % sample
            isi     = median(diff(obj.timestamps));
            duration= obj.timestamps(end)-obj.timestamps(1)+isi;
        end
    end
end
