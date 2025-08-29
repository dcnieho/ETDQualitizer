function [area, orientation, ax1, ax2, aspect_ratio] = bcea(azi, ele, P)
%BCEA Bivariate Contour Ellipse Area
%
%   result = BCEA(azi, ele, P) computes the BCEA and ellipse parameters
%   for gaze precision.
%
%   Inputs:
%       azi - Azimuth values in degrees (vector)
%       ele - Elevation values in degrees (vector)
%       P   - Cumulative probability (default: 0.68)
%
%   Outputs:
%       result - Struct with fields:
%           area         - BCEA value
%           orientation  - Orientation of the ellipse (degrees)
%           ax1          - Length of major axis
%           ax2          - Length of minor axis
%           aspect_ratio - Ratio of major to minor axis
%
%   Example:
%       result = bcea(randn(100,1), randn(100,1))

arguments
    azi (:,1) {mustBeNumeric}
    ele (:,1) {mustBeNumeric}
    P   (1,1) {mustBeNumeric} = 0.68    % cumulative probability of area under the multivariate normal
end

k    = log(1/(1-P));

valid   = ~(isnan(azi) | isnan(ele));

std_x   = std(azi,0,'omitnan');
std_y   = std(ele,0,'omitnan');
xx      = corrcoef(azi(valid), ele(valid));
rho     = xx(1,2);
area    = 2*k*pi*std_x*std_y*sqrt(1-rho.^2);

% calculate orientation of the bivariate normal ellipse
% (see
% https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Geometric_interpretation)
% and aspect ratio of axes. Note that an axis is half the
% diameter of ellipse along that direction.
% note that v and d outputs below are reordered versions of a
% and c in [a,b,c]=pca([obj.x(qValid) obj.y(qValid)]);
[v,d] = eig(cov(azi(valid), ele(valid)));
[~,i] = max(diag(d));
orientation = atan2(v(2,i), v(1,i))/pi*180;
ax1 = sqrt(k*d(i,i));
ax2 = sqrt(k*d(3-i,3-i));
aspect_ratio = max([ax1 ax2])/min([ax1 ax2]);
% sanity check: this (formula for area of ellipse) should
% closely match directly computed area from above
% 2*pi*ax1*ax2

