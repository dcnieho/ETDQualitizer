function [area, orientation, ax1, ax2, aspect_ratio] = bcea(x, y, P)
arguments
    x   (:,1) {mustBeNumeric}
    y   (:,1) {mustBeNumeric}
    P   (1,1) {mustBeNumeric} = 0.68    % cumulative probability of area under the multivariate normal
end
    k    = log(1/(1-P));

    valid   = ~(isnan(x) | isnan(y));

    std_x   = std(x,0,'omitnan');
    std_y   = std(y,0,'omitnan');
    xx      = corrcoef(x(valid), y(valid));
    rho     = xx(1,2);
    area    = 2*k*pi*std_x*std_y*sqrt(1-rho.^2);

    % calculate orientation of the bivariate normal ellipse
    % (see
    % https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Geometric_interpretation)
    % and aspect ratio of axes. Note that an axis is half the
    % diameter of ellipse along that direction.
    % note that v and d outputs below are reordered versions of a
    % and c in [a,b,c]=pca([obj.x(qValid) obj.y(qValid)]);
    [v,d] = eig(cov(x(valid), y(valid)));
    [~,i] = max(diag(d));
    orientation = atan2(v(2,i), v(1,i))/pi*180;
    ax1 = sqrt(k*d(i,i));
    ax2 = sqrt(k*d(3-i,3-i));
    aspect_ratio = max([ax1 ax2])/min([ax1 ax2]);
    % sanity check: this (formula for area of ellipse) should
    % closely match directly computed area from above
    % 2*pi*ax1*ax2
end
