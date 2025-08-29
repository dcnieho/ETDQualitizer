function freq = effective_frequency(azi, ele, duration)
%EFFECTIVE_FREQUENCY Compute Effective Sampling Frequency
%
%   Calculates effective frequency based on valid samples.
%
%   Syntax:
%       freq = effective_frequency(azi, ele, duration)
%
%   Inputs:
%       azi      - Azimuth values (vector)
%       ele      - Elevation values (vector)
%       duration - Duration in seconds (scalar)
%
%   Output:
%       freq - Effective frequency in Hz
%
%   Example:
%       freq = effective_frequency([1 NaN 3], [1 2 NaN], 1)

N_valid = sum(~(isnan(azi) | isnan(ele)));
freq    = N_valid/duration;
