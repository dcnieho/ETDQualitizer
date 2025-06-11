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

            % get unit vectors for gaze and target
            [g_x,g_y,g_z] = Fick_to_cartesian(   obj.x    ,    obj.y);
            [t_x,t_y,t_z] = Fick_to_cartesian(target_x_deg, target_y_deg);
            % calculate angular offset for each sample using dot product
            offsets     = acos(dot([g_x,g_y,g_z], repmat([t_x,t_y,t_z],length(g_x),1),2));
            % calculate on-screen orientation so we can decompose offset into x and y
            direction   = atan2(g_y./g_z-t_y./t_z, g_x./g_z-t_x./t_z);   % compute direction on tangent screen (divide by z to project to screen at 1m)
            offsets_2D  = offsets.*[cos(direction), sin(direction)]*180/pi;
            % calculate mean horizontal and vertical offset
            offset_x    = central_tendency_fun(offsets_2D(:,1));
            offset_y    = central_tendency_fun(offsets_2D(:,2));
            % calculate offset of centroid
            offset      = hypot(offset_x, offset_y);
        end

        function [rms, rms_x, rms_y] = precision_RMS_S2S(obj, central_tendency_fun)
            arguments
                obj
                central_tendency_fun(1,1) {mustBeA(central_tendency_fun,'function_handle')} = @(x) mean(x,'omitnan')
            end
            [rms, rms_x, rms_y] = RMS_S2S_impl(obj.x, obj.y, central_tendency_fun);
        end

        function [std_, std_x, std_y] = precision_STD(obj)
            [std_, std_x, std_y] = STD_impl(obj.x, obj.y);
        end

        function [area, orientation, ax1, ax2, aspect_ratio] = precision_BCEA(obj, P)
            arguments
                obj
                P   (1,1) {mustBeNumeric} = 0.68    % for BCEA: cumulative probability of area under the multivariate normal
            end
            [area, orientation, ax1, ax2, aspect_ratio] = BCEA_impl(obj.x, obj.y, P);
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
                    fun =  @RMS_S2S_impl;
                case 'STD'
                    fun =  @STD_impl;
                case 'BCEA'
                    fun =  @BCEA_impl;
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

        function loss_percentage = data_loss_percentage(obj)
            missing         = isnan(obj.x) | isnan(obj.y);
            loss_percentage = sum(missing)/length(missing)*100;
        end

        function freq = effective_frequency(obj)
            % rate of valid samples
            valid   = ~(isnan(obj.x) | isnan(obj.y));
            % to get duration right, we need to include duration of last
            % sample
            isi     = median(diff(obj.timestamps));
            freq    = sum(valid)/(obj.timestamps(end)-obj.timestamps(1)+isi);
        end
    end
end


% helpers
function [x,y,z] = Fick_to_cartesian(azi, ele, r)
arguments
    azi (:,1) {mustBeNumeric}
    ele (:,1) {mustBeNumeric}
    r   (:,1) {mustBeNumeric} = 1.
end
    r_cos_ele = r .* cosd(ele);

    x = r_cos_ele .* sind(azi);
    y =         r .* sind(ele);
    z = r_cos_ele .* cosd(azi);
end

function [rms, rms_x, rms_y] = RMS_S2S_impl(x, y, central_tendency_fun)
arguments
    x                   (:,1) {mustBeNumeric}
    y                   (:,1) {mustBeNumeric}
    central_tendency_fun(1,1) {mustBeA(central_tendency_fun,'function_handle')} = @(x) mean(x,'omitnan')
end
    x_diff  = diff(x).^2;
    y_diff  = diff(y).^2;
    rms_x   = sqrt(central_tendency_fun(x_diff));
    rms_y   = sqrt(central_tendency_fun(y_diff));
    % N.B.: cannot simplify to hypot(rms_x, rms_y)
    % as that is only equivalent when mean() is used as central tendency estimator
    rms     = sqrt(central_tendency_fun(x_diff + y_diff));
end

function [std_, std_x, std_y] = STD_impl(x, y)
arguments
    x   (:,1) {mustBeNumeric}
    y   (:,1) {mustBeNumeric}
end
    std_x   = std(x,1,'omitnan');
    std_y   = std(y,1,'omitnan');
    std_    = hypot(std_x, std_y);
end

function [area, orientation, ax1, ax2, aspect_ratio] = BCEA_impl(x, y, P)
arguments
    x   (:,1) {mustBeNumeric}
    y   (:,1) {mustBeNumeric}
    P   (1,1) {mustBeNumeric} = 0.68    % cumulative probability of area under the multivariate normal
end
    k    = log(1/(1-P));

    valid   = ~(isnan(x) | isnan(y));

    std_x   = std(x,0,'omitnan');
    std_y   = std(y,0,'omitnan');
    xx      = corrcoef(x(valid), y(valid));
    rho     = xx(1,2);
    area    = 2*k*pi*std_x*std_y*sqrt(1-rho.^2);

    % calculate orientation of the bivariate normal ellipse
    % (see
    % https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Geometric_interpretation)
    % and aspect ratio of axes. Note that an axis is half the
    % diameter of ellipse along that direction.
    % note that v and d outputs below are reordered versions of a
    % and c in [a,b,c]=pca([obj.x(qValid) obj.y(qValid)]);
    [v,d] = eig(cov(x(valid), y(valid)));
    [~,i] = max(diag(d));
    orientation = atan2(v(2,i), v(1,i))/pi*180;
    ax1 = sqrt(k*d(i,i));
    ax2 = sqrt(k*d(3-i,3-i));
    aspect_ratio = max([ax1 ax2])/min([ax1 ax2]);
    % sanity check: this (formula for area of ellipse) should
    % closely match directly computed area from above
    % 2*pi*ax1*ax2
end
