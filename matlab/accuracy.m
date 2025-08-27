function [offset, offset_x, offset_y] = accuracy(x, y, target_x_deg, target_y_deg, central_tendency_fun)
%ACCURACY Compute Gaze Accuracy
%
%   [offset, offset_x, offset_y] = ACCURACY(x, y, target_x_deg, target_y_deg, central_tendency_fun)
%   calculates the angular offset between gaze and target directions.
%
%   Inputs:
%       x, y - Gaze azimuth and elevation in degrees
%       target_x_deg, target_y_deg - Target azimuth and elevation in degrees
%       central_tendency_fun - Function handle for central tendency (default: @mean)
%
%   Outputs:
%       offset - Total angular offset in degrees
%       offset_x - Horizontal offset in degrees
%       offset_y - Vertical offset in degrees
%
%   Example:
%       [offset, offset_x, offset_y] = accuracy([1,2], [1,2], 0, 0)

arguments
    x                   (:,1) {mustBeNumeric}
    y                   (:,1) {mustBeNumeric}
    target_x_deg        (1,1) {mustBeNumeric}
    target_y_deg        (1,1) {mustBeNumeric}
    central_tendency_fun(1,1) {mustBeA(central_tendency_fun,'function_handle')} = @(x) mean(x,'omitnan')
end

% get unit vectors for gaze and target
[g_x,g_y,g_z] = Fick_to_vector(       x    ,        y);
[t_x,t_y,t_z] = Fick_to_vector(target_x_deg, target_y_deg);
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

