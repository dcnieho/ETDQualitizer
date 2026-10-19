function g = frechet_median_on_sphere(vectors, tol, max_iter)
%FRECHET_MEDIAN_ON_SPHERE Compute the spherical Frechet median of 3D vectors.
%
%   g = FRECHET_MEDIAN_ON_SPHERE(vectors) computes the unit direction that
%   minimizes the sum of angular distances to the rows of vectors. Non-finite
%   rows are ignored. Optional arguments set the convergence tolerance and
%   maximum number of iterations.
%
%   Inputs:
%       vectors - N-by-3 matrix of 3D direction vectors
%       tol - Convergence tolerance (default: 1e-12)
%       max_iter - Maximum number of iterations (default: 128)
%
%   Output:
%       g - 1-by-3 spherical Frechet median; [NaN, NaN, NaN] if no finite
%           vectors are provided
%
%   Example:
%       g = frechet_median_on_sphere([0, 0, 1; 0, 0, 1])

arguments
    vectors  (:,3) {mustBeNumeric}
    tol      (1,1) {mustBeNumeric, mustBePositive} = 1e-12
    max_iter (1,1) {mustBeInteger, mustBePositive} = 128
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