function [offset, offset_azi, offset_ele] = accuracy(azi, ele, target_azi, target_ele, central_tendency_fun)
%ACCURACY Compute Gaze Accuracy
%
%   [offset, offset_azi, offset_ele] = ACCURACY(azi, ele, target_azi, target_ele, central_tendency_fun)
%   calculates the angular offset between gaze and target directions.
%
%   Inputs:
%       azi, ele - Gaze azimuth and elevation in degrees
%       target_azi, target_ele - Target azimuth and elevation in degrees
%       central_tendency_fun - Function handle for central tendency (default: @mean)
%
%   Outputs:
%       offset - Total angular offset in degrees
%       offset_azi - Horizontal offset in degrees
%       offset_ele - Vertical offset in degrees
%
%   Example:
%       [offset, offset_azi, offset_ele] = accuracy([1,2], [1,2], 0, 0)

arguments
    azi                 (:,1) {mustBeNumeric}
    ele                 (:,1) {mustBeNumeric}
    target_azi          (1,1) {mustBeNumeric}
    target_ele          (1,1) {mustBeNumeric}
    central_tendency_fun(1,1) {mustBeA(central_tendency_fun,'function_handle')} = @(x) mean(x,'omitnan')
end

% convert gaze directions to unit vectors
[gx, gy, gz] = Fick_to_vector(azi, ele);

% compute central gaze direction in 3D (e.g. mean direction)
g = [central_tendency_fun(gx), ...
     central_tendency_fun(gy), ...
     central_tendency_fun(gz)];

% normalize to obtain a unit direction vector
g = g / norm(g);

% precompute trigonometric terms for target orientation
ca = cosd(target_azi); sa = sind(target_azi);
ce = cosd(target_ele); se = sind(target_ele);

% express central gaze direction in a target-centered frame
% (i.e., rotate such that the target lies straight ahead)
x_rel =  ca*g(1) - sa*g(3);
y_rel =  ce*g(2) - se*(sa*g(1) + ca*g(3));
z_rel =  se*g(2) + ce*(sa*g(1) + ca*g(3));

% decompose relative direction into Fick components
[offset_azi, offset_ele] = vector_to_Fick(x_rel, y_rel, z_rel);

% compute total angular offset (angle between gaze and target)
offset = atan2d(hypot(x_rel, y_rel), z_rel);
