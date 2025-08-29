function [rms, rms_azi, rms_ele] = rms_s2s(azi, ele, central_tendency_fun)
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
%       rms     - Total RMS of sample-to-sample distances (degrees)
%       rms_azi - RMS of azimuthal component (degrees)
%       rms_ele - RMS of elevation component (degrees)
%
%   Example:
%       [rms, rms_a, rms_e] = rms_s2s([1, 2, 3], [1, 2, 3])

arguments
    azi                 (:,1) {mustBeNumeric}
    ele                 (:,1) {mustBeNumeric}
    central_tendency_fun(1,1) {mustBeA(central_tendency_fun,'function_handle')} = @(x) mean(x,'omitnan')
end

a_diff  = diff(azi).^2;
e_diff  = diff(ele).^2;
rms_azi = sqrt(central_tendency_fun(a_diff));
rms_ele = sqrt(central_tendency_fun(e_diff));
% N.B.: cannot simplify to hypot(rms_azi, rms_ele)
% as that is only equivalent when mean() is used as central tendency estimator
rms     = sqrt(central_tendency_fun(a_diff + e_diff));
