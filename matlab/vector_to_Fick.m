function [azi, ele] = vector_to_Fick(x, y, z)
arguments
    x (:,1) {mustBeNumeric}
    y (:,1) {mustBeNumeric}
    z (:,1) {mustBeNumeric}
end
    azi = atan2(x,         z)*180/pi;
    ele = atan2(y,hypot(x,z))*180/pi;
end
