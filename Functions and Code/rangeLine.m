
function [y, percent_range] = rangeLine(fig_handle, offset_percent, subtract)
% [y, percent_range] = rangeLine(fig_handle,offset_percent,subtract)
%
% Find an appropriate y-value for a point just below the top of the graph 
% -- used for plotting the range of a stimulus, like the laser onset or odor onset
%
% INPUTS
% 'fig_handle' : handle to the figure in question (must have the
%       appropriate subplot already activated if there are multiple)
% 'offset_percent' : percent of total y-axis height for the desired y-value offset
%       (default value = 5%)
% 'subtract' : find value lower or higher than the current Y-Limits
%       subtract = true (default) if the line wants to be within the current ylims
%       subtract = false if you want to go above current ylimits
%
% OUTPUTS
% 'y' :  the new y value for the given offset percentage
% 'percent_range' :  range of the given percent offset (e.g. 5% of
%       the current range)
%
% ES Dickinson, University of Washington, 2020

%%
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


