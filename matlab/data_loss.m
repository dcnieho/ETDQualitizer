function loss_percentage = data_loss(azi, ele)
%DATA_LOSS Compute Data Loss Percentage
%
%   Calculates the percentage of missing gaze samples.
%
%   Syntax:
%       loss_percentage = data_loss(azi, ele)
%
%   Inputs:
%       azi - Azimuth values (vector)
%       ele - Elevation values (vector)
%
%   Output:
%       loss_percentage - Percentage of missing samples
%
%   Example:
%       loss_percentage = data_loss([1 NaN 3], [1 2 NaN])

missing         = isnan(azi) | isnan(ele);
loss_percentage = sum(missing)/length(missing)*100;
