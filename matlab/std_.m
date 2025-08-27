function [std_, std_x, std_y] = std_(x, y)
%STD_ Standard Deviation of Gaze Samples
%
%   result = STD_(azi, ele) computes the standard deviation of azimuth
%   and elevation gaze samples.
%
%   Inputs:
%       azi - Azimuth values in degrees (vector)
%       ele - Elevation values in degrees (vector)
%
%   Outputs:
%       result - Struct with fields:
%           std    - Total standard deviation (degrees)
%           std_a  - Standard deviation of azimuth (degrees)
%           std_e  - Standard deviation of elevation (degrees)
%
%   Example:
%       result = std_gaze([1, 2, 3], [1, 2, 3])

arguments
    x   (:,1) {mustBeNumeric}
    y   (:,1) {mustBeNumeric}
end

std_x   = std(x,1,'omitnan');
std_y   = std(y,1,'omitnan');
std_    = hypot(std_x, std_y);
