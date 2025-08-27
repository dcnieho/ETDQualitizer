function [azi, ele] = vector_to_Fick(x, y, z)
%VECTOR_TO_FICK Convert 3D Vector to Fick Angles
%
%   [azi, ele] = VECTOR_TO_FICK(x, y, z) converts a 3D vector to azimuth and
%   elevation angles (in degrees).
%
%   Inputs:
%       x, y, z - Components of the 3D vector
%
%   Outputs:
%       azi - Azimuth angle in degrees
%       ele - Elevation angle in degrees
%
%   Example:
%       [azi, ele] = vector_to_Fick(0.5, 0.2, 0.8)

arguments
    x (:,1) {mustBeNumeric}
    y (:,1) {mustBeNumeric}
    z (:,1) {mustBeNumeric}
end

azi = atan2(x,         z)*180/pi;
ele = atan2(y,hypot(x,z))*180/pi;

