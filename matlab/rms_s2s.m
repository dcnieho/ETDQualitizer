function [rms, rms_x, rms_y] = rms_s2s(x, y, central_tendency_fun)
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
