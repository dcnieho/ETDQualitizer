function [x,y,z] = Fick_to_vector(azi, ele, r)
arguments
    azi (:,1) {mustBeNumeric}
    ele (:,1) {mustBeNumeric}
    r   (:,1) {mustBeNumeric} = 1.
end
    r_cos_ele = r .* cosd(ele);

    x = r_cos_ele .* sind(azi);
    y =         r .* sind(ele);
    z = r_cos_ele .* cosd(azi);
end
