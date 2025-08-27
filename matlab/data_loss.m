function loss_percentage = data_loss(x, y)
%DATA_LOSS Compute Data Loss Percentage
%
%   Calculates the percentage of missing gaze samples.
%
%   Syntax:
%       loss_percentage = data_loss(x, y)
%
%   Inputs:
%       x - Azimuth values (vector)
%       y - Elevation values (vector)
%
%   Output:
%       loss_percentage - Percentage of missing samples
%
%   Example:
%       loss_percentage = data_loss([1 NaN 3], [1 2 NaN])

missing         = isnan(x) | isnan(y);
loss_percentage = sum(missing)/length(missing)*100;
