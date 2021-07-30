
function y = rangeLine(fig_handle)
% y = rangeLine(fig_handle)
% get the y value for a point just below the 
% top of the graph -- used for plotting the
% range of a stimulus, like the laser
% 
% ES Dickinson, University of Washington, 2020


figure(fig_handle); %activate the desired figure

offset_percent = 5;
yMax = max(ylim);

yRange = range(ylim);
offset = (yRange/100)*offset_percent;

y = yMax-(abs(offset));


end