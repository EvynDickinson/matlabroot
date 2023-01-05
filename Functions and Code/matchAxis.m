
function fig = matchAxis(fig,independent)
% fig = matchAxis(fig,independent)
% makes the x and y axes equal to the larger of the axes
% independent = true if you want to match the x and y axes independently; false if
% you want x&y to match across all subplots

[y_lim, x_lim] = deal([]);
figAxes = findall(fig,'type','axes');
% get axis limits:
for ii = 1:size(figAxes,1)
    x_lim = autoCat(x_lim,figAxes(ii).XLim);
    y_lim = autoCat(y_lim,figAxes(ii).YLim);
end
    
    
if nargin==2 && independent==true
    for ii = 1:size(figAxes,1)
        figAxes(ii).XLim = [min(min(x_lim)),max(max(x_lim))];
        figAxes(ii).YLim = [min(min(y_lim)),max(max(y_lim))];
    end
else
    % find min and max:
    low = min(min([x_lim; y_lim]));
    high = max(max([x_lim; y_lim]));

    for ii = 1:size(figAxes,1)
        figAxes(ii).XLim = [low,high];
        figAxes(ii).YLim = [low,high];
    end
end












