
function [y, percent_range] = rangeLine(fig_handle,offset_percent,subtract)
% [y, percent_range] = rangeLine(fig_handle,offset_percent,subtract)
% get the y value for a point just below the 
% top of the graph -- used for plotting the
% range of a stimulus, like the laser onset or odor onset
% subtract = true (default) if the line wants to 
% be lower than the current ylims
% subtract = false if you want to go above current ylimits
% y =  the new y value for the given offset percentage
% percent_range = what is the range of the given percent offset (e.g. 5% of
% the current range)
%
% ES Dickinson, University of Washington, 2020


figure(fig_handle); %activate the desired figure
if nargin==1
    offset_percent = 5;
end
yMax = max(ylim);

yRange = range(ylim);
offset = (yRange/100)*offset_percent;

if nargin<3
    subtract = true;
end
    
if subtract
    y = yMax-(abs(offset));
else
    y = yMax+(abs(offset));
end

percent_range = offset;


end