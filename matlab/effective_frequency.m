function freq = effective_frequency(x, y, duration)
%EFFECTIVE_FREQUENCY Compute Effective Sampling Frequency
%
%   Calculates effective frequency based on valid samples.
%
%   Syntax:
%       freq = effective_frequency(x, y, duration)
%
%   Inputs:
%       x        - Azimuth values (vector)
%       y        - Elevation values (vector)
%       duration - Duration in seconds (scalar)
%
%   Output:
%       freq - Effective frequency in Hz
%
%   Example:
%       freq = effective_frequency([1 NaN 3], [1 2 NaN], 1)

N_valid = sum(~(isnan(x) | isnan(y)));
freq    = N_valid/duration;
