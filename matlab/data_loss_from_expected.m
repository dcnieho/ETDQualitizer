function loss_percentage = data_loss_from_expected(x, y, duration, frequency)
%DATA_LOSS_FROM_EXPECTED Compute Data Loss from Expected Sample Count
%
%   Calculates data loss based on expected number of samples.
%
%   Syntax:
%       loss_percentage = data_loss_from_expected(x, y, duration, frequency)
%
%   Inputs:
%       x         - Azimuth values (vector)
%       y         - Elevation values (vector)
%       duration  - Duration in seconds (scalar)
%       frequency - Sampling frequency in Hz (scalar)
%
%   Output:
%       loss_percentage - Percentage of data loss
%
%   Example:
%       loss_percentage = data_loss_from_expected([1 NaN 3], [1 2 NaN], 1, 3)

N_valid         = sum(~(isnan(x) | isnan(y)));
loss_percentage = (1-N_valid/(duration*frequency))*100;
