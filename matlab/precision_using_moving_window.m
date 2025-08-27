function precision = precision_using_moving_window(x, y, window_length, metric, input_args, aggregation_fun)
%PRECISION_USING_MOVING_WINDOW Precision Using Moving Window
%
%   Computes gaze precision using a moving window and selected metric.
%
%   Syntax:
%       precision = precision_using_moving_window(x, y, window_length, metric)
%       precision = precision_using_moving_window(x, y, window_length, metric, input_args)
%       precision = precision_using_moving_window(x, y, window_length, metric, input_args, aggregation_fun)
%
%   Inputs:
%       x               - Azimuth values (column vector)
%       y               - Elevation values (column vector)
%       window_length   - Window size in samples (scalar integer)
%       metric          - Precision metric: 'RMS-S2S', 'STD', or 'BCEA'
%       input_args      - Cell array of additional arguments passed to metric function (optional)
%       aggregation_fun - Function handle to aggregate precision values (default: median ignoring NaNs)
%
%   Output:
%       precision - Aggregated precision value
%
%   Example:
%       precision = precision_using_moving_window(randn(100,1), randn(100,1), 10, 'STD')

arguments
    x               (:,1) {mustBeNumeric}
    y               (:,1) {mustBeNumeric}
    window_length   (1,1) {mustBeInteger}
    metric          (1,:) char {mustBeMember(metric,{'RMS-S2S','STD','BCEA'})}
    input_args      (1,:) cell = {}
    aggregation_fun (1,1) {mustBeA(aggregation_fun,'function_handle')} = @(x) median(x,'omitnan')
end

switch metric
    case 'RMS-S2S'
        fun =  @rms_s2s;
    case 'STD'
        fun =  @std_;
    case 'BCEA'
        fun =  @bcea;
end

% get number of samples in data
ns  = length(x);

if window_length < ns % if number of samples in data exceeds window size
    values = nan(1,ns-window_length); % pre-allocate
    for p=1:ns-window_length+1
        values(p) = fun(x(p:p+window_length-1), y(p:p+window_length-1), input_args{:});
    end
    precision = aggregation_fun(values);
else % if too few samples in data
    precision = NaN;
end
