function loss_percentage = data_loss(a, b)
%DATA_LOSS Compute Data Loss Percentage
%
%   Calculates the percentage of missing gaze samples.
%
%   Syntax:
%       loss_percentage = data_loss(a, b)
%
%   Inputs:
%       a - Horizontal gaze values (vector, e.g. azimuth or horizontal coordinate in pixels or mm).
%       b - Vertical gaze values (vector, e.g. azimuth or horizontal coordinate in pixels or mm).
%
%   Output:
%       loss_percentage - Percentage of missing samples
%
%   Example:
%       loss_percentage = data_loss([1 NaN 3], [1 2 NaN])

missing         = isnan(a) | isnan(b);
loss_percentage = sum(missing)/length(missing)*100;
