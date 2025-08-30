function freq = effective_frequency(a, b, duration)
%EFFECTIVE_FREQUENCY Compute Effective Sampling Frequency
%
%   Calculates effective frequency based on valid samples.
%
%   Syntax:
%       freq = effective_frequency(a, b, duration)
%
%   Inputs:
%       a - Horizontal gaze values (vector, e.g. azimuth or horizontal coordinate in pixels or mm).
%       b - Vertical gaze values (vector, e.g. azimuth or horizontal coordinate in pixels or mm).
%       duration - Duration in seconds (scalar)
%
%   Output:
%       freq - Effective frequency in Hz
%
%   Example:
%       freq = effective_frequency([1 NaN 3], [1 2 NaN], 1)

N_valid = sum(~(isnan(a) | isnan(b)));
freq    = N_valid/duration;
