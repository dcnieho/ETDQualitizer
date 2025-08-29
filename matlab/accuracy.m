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

% get unit vectors for gaze and target
[g_x,g_y,g_z] = Fick_to_vector(       azi    ,        ele);
[t_x,t_y,t_z] = Fick_to_vector(target_azi, target_ele);
% calculate angular offset for each sample using dot product
offsets     = acos(dot([g_x,g_y,g_z], repmat([t_x,t_y,t_z],length(g_x),1),2));
% calculate on-screen orientation so we can decompose offset into x and y
direction   = atan2(g_y./g_z-t_y./t_z, g_x./g_z-t_x./t_z);   % compute direction on tangent screen (divide by z to project to screen at 1m)
offsets_2D  = offsets.*[cos(direction), sin(direction)]*180/pi;
% calculate mean horizontal and vertical offset
offset_azi  = central_tendency_fun(offsets_2D(:,1));
offset_ele  = central_tendency_fun(offsets_2D(:,2));
% calculate offset of centroid
offset      = hypot(offset_azi, offset_ele);

