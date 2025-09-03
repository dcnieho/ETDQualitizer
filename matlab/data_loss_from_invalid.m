function loss_percentage = data_loss_from_invalid(a, b)
%DATA_LOSS_FROM_INVALID Compute Data Loss from number of invalid samples
%
%   Calculates the percentage of missing gaze samples.
%
%   Syntax:
%       loss_percentage = data_loss_from_invalid(a, b)
%
%   Inputs:
%       a - Horizontal gaze values (vector, e.g. azimuth or horizontal coordinate in pixels or mm).
%       b - Vertical gaze values (vector, e.g. azimuth or horizontal coordinate in pixels or mm).
%
%   Output:
%       loss_percentage - Percentage of missing samples
%
%   Example:
%       loss_percentage = data_loss_from_invalid([1 NaN 3], [1 2 NaN])

missing         = isnan(a) | isnan(b);
loss_percentage = sum(missing)/length(missing)*100;
