function [std_, std_x, std_y] = std_(x, y)
arguments
    x   (:,1) {mustBeNumeric}
    y   (:,1) {mustBeNumeric}
end
    std_x   = std(x,1,'omitnan');
    std_y   = std(y,1,'omitnan');
    std_    = hypot(std_x, std_y);
end
