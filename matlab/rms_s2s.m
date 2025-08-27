function [rms, rms_x, rms_y] = rms_s2s(x, y, central_tendency_fun)
%RMS_S2S RMS of Sample-to-Sample Differences
%
%   [rms, rms_a, rms_e] = RMS_S2S(azi, ele) computes the root mean square
%   (RMS) of differences between successive gaze samples in azimuth and elevation.
%
%   Inputs:
%       azi - Azimuth values in degrees (vector)
%       ele - Elevation values in degrees (vector)
%       central_tendency_fun - Function handle to compute central tendency
%                               (default: @mean)
%
%   Outputs:
%       rms    - Total RMS of sample-to-sample distances (degrees)
%       rms_a  - RMS of azimuthal component (degrees)
%       rms_e  - RMS of elevation component (degrees)
%
%   Example:
%       [rms, rms_a, rms_e] = rms_s2s([1, 2, 3], [1, 2, 3])

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
