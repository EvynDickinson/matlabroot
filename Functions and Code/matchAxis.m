
function fig = matchAxis(fig,scale,independent,);
% fig = matchAxis(fig,~)
% makes the x and y axes equal to the larger of the axes
% if there are multiple subplots, it will align those as well 
% r_c = rows and columns in the figure
% sp = subplot indexes

% 

if nargin==1

    [y_lim, x_lim] = deal([]);
    figAxes = findall(fig,'type','axes');
    % get axis limits:
    for ii = 1:size(figAxes,1)
        x_lim = autoCat(x_lim,figAxes(ii).XLim);
        y_lim = autoCat(y_lim,figAxes(ii).YLim);
    end
    
    % find min and max:
    low = min(min([x_lim; y_lim]));
    high = max(max([x_lim; y_lim]));

    for ii = 1:size(figAxes,1)
        figAxes(ii).XLim = [low,high];
        figAxes(ii).YLim = [low,high];
    end

elseif nargin == 2  %normalize axes for each subplot independently
    % TODO
end

% % could also add a function here that would set all the subplot axes to the
% % same even if the x and y are different
% 
% if nargin==2
%     [y_lim, x_lim] = deal([]);
%     figAxes = findall(fig,'type','axes');
%     % get axis limits:
%     for ii = 1:size(figAxes,1)
%         x_lim = autoCat(x_lim,figAxes(ii).XLim);
%         y_lim = autoCat(y_lim,figAxes(ii).YLim);
%     end
%     
%     % find min and max:
%     xRange = 
%     low = min(min([x_lim; y_lim]));
%     high = max(max([x_lim; y_lim]));
% 
% 
% 











