function [offset, offset_azi, offset_ele] = accuracy(azi, ele, target_azi, target_ele, central_tendency_fun)
%ACCURACY Compute Gaze Accuracy
%
%   [offset, offset_azi, offset_ele] = ACCURACY(azi, ele, target_azi, target_ele, central_tendency_fun)
%   calculates the angular offset between gaze and target directions.
%
%   Inputs:
%       azi, ele - Gaze azimuth and elevation in degrees
%       target_azi, target_ele - Target azimuth and elevation in degrees
%       central_tendency_fun - Function handle for central tendency (default: @mean).
%           Median-like functions are interpreted as a spherical Fréchet
%           median instead of a component-wise Cartesian median.
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

% convert gaze directions to unit vectors
[gx, gy, gz] = Fick_to_vector(azi, ele);

% compute central gaze direction in 3D (e.g. mean direction)
if uses_frechet_median(central_tendency_fun)
    g = frechet_median_on_sphere([gx, gy, gz]);
else
    g = [central_tendency_fun(gx), ...
         central_tendency_fun(gy), ...
         central_tendency_fun(gz)];
end

% normalize to obtain a unit direction vector
g = g / norm(g);

% precompute trigonometric terms for target orientation
ca = cosd(target_azi); sa = sind(target_azi);
ce = cosd(target_ele); se = sind(target_ele);

% express central gaze direction in a target-centered frame
% (i.e., rotate such that the target lies straight ahead)
x_rel =  ca*g(1) - sa*g(3);
y_rel =  ce*g(2) - se*(sa*g(1) + ca*g(3));
z_rel =  se*g(2) + ce*(sa*g(1) + ca*g(3));

% decompose relative direction into Fick components
[offset_azi, offset_ele] = vector_to_Fick(x_rel, y_rel, z_rel);

% compute total angular offset (angle between gaze and target)
offset = atan2d(hypot(x_rel, y_rel), z_rel);

end

function tf = uses_frechet_median(central_tendency_fun)
name = lower(regexprep(func2str(central_tendency_fun), '\s+', ''));
tf = strcmp(name, 'median') || strcmp(name, '@median') || ...
     strcmp(name, 'nanmedian') || strcmp(name, '@nanmedian') || ...
     strcmp(name, 'frechet_median') || strcmp(name, '@frechet_median') || ...
     strcmp(name, 'frechetmedian') || strcmp(name, '@frechetmedian') || ...
     ~isempty(regexp(name, '^@\(.*\)median\(', 'once')) || ...
     ~isempty(regexp(name, '^@\(.*\)nanmedian\(', 'once'));
end

function g = frechet_median_on_sphere(vectors, tol, max_iter)
if nargin < 2
    tol = 1e-12;
end
if nargin < 3
    max_iter = 128;
end

vectors = vectors(all(isfinite(vectors), 2), :);
if isempty(vectors)
    g = [NaN, NaN, NaN];
    return
end

g = mean(vectors, 1, 'omitnan');
g_norm = norm(g);
if g_norm < tol
    g = vectors(1, :);
else
    g = g / g_norm;
end

current_value = sphere_objective(g, vectors);

for iter = 1:max_iter
    dots = max(-1, min(1, vectors * g.'));
    angles = acos(dots);
    is_differentiable = angles > tol & angles < pi - tol;
    if ~any(is_differentiable)
        break
    end

    tangent = vectors(is_differentiable, :) - dots(is_differentiable) .* g;
    tangent = tangent ./ sin(angles(is_differentiable));
    direction = sum(tangent, 1);
    direction_norm = norm(direction);
    if direction_norm < tol
        break
    end
    direction = direction / direction_norm;

    step = pi / 4;
    improved = false;
    while step > tol
        candidate = cos(step) * g + sin(step) * direction;
        candidate = candidate / norm(candidate);
        candidate_value = sphere_objective(candidate, vectors);
        if candidate_value + tol < current_value
            g = candidate;
            if current_value - candidate_value < tol
                return
            end
            current_value = candidate_value;
            improved = true;
            break
        end
        step = step / 2;
    end

    if ~improved
        break
    end
end
end

function value = sphere_objective(candidate, vectors)
dots = max(-1, min(1, vectors * candidate.'));
value = sum(acos(dots));
end
