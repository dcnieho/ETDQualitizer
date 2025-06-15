function [offset, offset_x, offset_y] = accuracy(x, y, target_x_deg, target_y_deg, central_tendency_fun)
arguments
    x                   (:,1) {mustBeNumeric}
    y                   (:,1) {mustBeNumeric}
    target_x_deg        (1,1) {mustBeNumeric}
    target_y_deg        (1,1) {mustBeNumeric}
    central_tendency_fun(1,1) {mustBeA(central_tendency_fun,'function_handle')} = @(x) mean(x,'omitnan')
end

% get unit vectors for gaze and target
[g_x,g_y,g_z] = Fick_to_cartesian(       x    ,        y);
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
