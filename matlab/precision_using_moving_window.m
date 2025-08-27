function precision = precision_using_moving_window(x, y, window_length, metric, input_args, aggregation_fun)
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
end
