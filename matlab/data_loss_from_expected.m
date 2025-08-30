function loss_percentage = data_loss_from_expected(a, b, duration, frequency)
%DATA_LOSS_FROM_EXPECTED Compute Data Loss from Expected Sample Count
%
%   Calculates data loss based on expected number of samples.
%
%   Syntax:
%       loss_percentage = data_loss_from_expected(a, b, duration, frequency)
%
%   Inputs:
%       a - Horizontal gaze values (vector, e.g. azimuth or horizontal coordinate in pixels or mm).
%       b - Vertical gaze values (vector, e.g. azimuth or horizontal coordinate in pixels or mm).
%       duration  - Duration in seconds (scalar)
%       frequency - Sampling frequency in Hz (scalar)
%
%   Output:
%       loss_percentage - Percentage of data loss
%
%   Example:
%       loss_percentage = data_loss_from_expected([1 NaN 3], [1 2 NaN], 1, 3)

N_valid         = sum(~(isnan(a) | isnan(b)));
loss_percentage = (1-N_valid/(duration*frequency))*100;
