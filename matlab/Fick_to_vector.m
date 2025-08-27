function [x,y,z] = Fick_to_vector(azi, ele, r)
%FICK_TO_VECTOR Convert Fick Angles to 3D Vector
%
%   [x, y, z] = FICK_TO_VECTOR(azi, ele, rho) converts azimuth and
%   elevation angles (in degrees) to a 3D unit vector. The radius rho is
%   optional and defaults to 1.0.
%
%   Inputs:
%       azi - Azimuth angle in degrees (scalar or vector)
%       ele - Elevation angle in degrees (scalar or vector)
%       rho - Radius (default: 1.0)
%
%   Outputs:
%       x, y, z - Components of the 3D vector
%
%   Example:
%       [x, y, z] = Fick_to_vector(30, 10

arguments
    azi (:,1) {mustBeNumeric}
    ele (:,1) {mustBeNumeric}
    r   (:,1) {mustBeNumeric} = 1.
end

r_cos_ele = r .* cosd(ele);

x = r_cos_ele .* sind(azi);
y =         r .* sind(ele);
z = r_cos_ele .* cosd(azi);
