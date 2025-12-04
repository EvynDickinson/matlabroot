function figHandle = formatFig(figHandle, BandWoption, sbplots, subplotInd)
% figHandle = formatFig(figHandle, BandWoption, sbplots, subplotInd);
% e.g. fig = formatFig(fig,true);
% turns a figure into a black background with with axes
% no BandWoption just separates the axes
% 
% inputs:
% figHandle --> handle for the figure
% BandWoption --> true||false for black background
% sbplots --> number of subplots *[row,columns]
% subplotInd --> subplot indexs if not 1:1
% e.g. subplotInd(1).idx = [1,3,4];
% 
% ES Dickinson,
% University of Washington, 2020

% Adjust the font sizes dynamically by the computer and version of matlab 
if ismac 
    labelSize = 15; 
else
    if strcmp(version('-release'),'2025b')
        labelSize = 15;
        axWidth = 2;
        tickDir = 'out';
    else
        labelSize = 10; 
        axWidth = 1;
        tickDir = 'out';
    end
end


% If black color is requested set background first:
% adjust colors to blk&white
backColor = 'w';
labelColor = 'k';
if nargin>=2
  if BandWoption==true
    set(figHandle, 'color', 'k') 
    labelColor = 'w';
    backColor = 'k';
  else
    set(figHandle, 'color', 'w') 
  end
else
    set(figHandle, 'color', 'w') 
end

% determine if there are subplots or just a single plot:
if nargin >= 3 
    %there are multiple plots
    ncol = sbplots(2);  %number of columns
    nrow = sbplots(1);  %number of rows
    nplots = length(findobj(figHandle,'type','axes'));
        % set properties for each subplot:
    for n = 1:nplots
        if nargin==4
            ii = subplotInd(n).idx;
        else
            ii = n;
        end
        subplot(nrow, ncol, ii)
        % get plot handles
        ax = gca(figHandle);
        ax.LineWidth = axWidth; %change axes lines to width of 1
        box off
        % set the subplot axes to the selected color scheme
        set(ax, 'color', backColor, 'YColor', labelColor, 'XColor', labelColor,'tickdir', tickDir)
        % set the font: 
        set(ax, 'FontName', 'Arial');

        % set the label sizes and color
        labelHandles = findall(ax, 'type', 'text', 'handlevisibility', 'off');
        set(labelHandles,'FontSize', labelSize, 'color', labelColor)
        set(gca, 'FontSize', labelSize)
    end    
else % just a single plot
   % get plot handles
    ax = gca(figHandle);
    ax.LineWidth = axWidth; %change axes lines to width of 1
    box off
    % set the subplot axes to the selected color scheme
    set(ax, 'color', backColor, 'YColor', labelColor, 'XColor', labelColor,'tickdir', tickDir)
    % set the font: 
    set(ax, 'FontName', 'Arial');
    
    % set the label sizes and color
    labelHandles = findall(ax, 'type', 'text', 'handlevisibility', 'off');
    set(labelHandles,'FontSize', labelSize, 'color', labelColor)
    set(gca, 'FontSize', labelSize)

end

clear ax


% 
% for n = 1:nplots
%     ax = subplot(nrow, ncol,n);
%     % separate figure axes:
%     xRange = ax.XLim;
%     xMin(n) = xRange(1);
%     xOffset = range(xRange)*0.03; % 3-percent range offset
%     ax.XLim(1) = xMin(n)-xOffset;
%     
%     yRange = ax.YLim;
%     yMin(n) = yRange(1);
%     yOffset = range(yRange)*0.05; % 3-percent range offset
%     ax.YLim(1) = yMin(n)-yOffset;
%     drawnow()
%      
% %     % set the X axis vertex start to its the orignial point
% %     origin = get(ax.XAxis.Axle,'VertexData');
% %     origin(1,1) = xMin(n);
% %     set(ax.XAxis.Axle,'VertexData',origin);
% % 
% %     % set the Y axis vertex start to its the orignial point
% %     origin = get(ax.YAxis.Axle,'VertexData');
% %     origin(2,1) = yMin(n);
% %     set(ax.YAxis.Axle,'VertexData',origin);
% 
% end


% % redraw the plot line
% for n = 1:nplots
%    subplot(nrow,ncol,n);
%    separateAxes
% end
% 


% 
% % redraw the plot line
% for n = 1:nplots
%     ax = subplot(nrow,ncol,n);
% %     separateAxes;  
% 
%     % set the X axis vertex start to its the orignial point
%     origin = get(ax.XAxis.Axle,'VertexData');
%     origin(1,1) = ax.XTick(1);
%     set(ax.XAxis.Axle,'VertexData',origin);
%     
%     % set the Y axis vertex start to its the orignial point
%     origin = get(ax.YAxis.Axle,'VertexData');
%     origin(2,1) = ax.YTick(1);
%     set(ax.YAxis.Axle,'VertexData',origin);
%     
% end
% 






% % adjust colors to blk&white
% if nargin>=2
%   if BandWoption==true
%     set(figHandle, 'color', 'k') 
%     set(ax, 'color', 'k', 'YColor', 'w', 'XColor', 'w')
%     labelColor = 'w';
%   end
% end
% 
% % set the tick labels, axes labels and title font and size
% textHandles = findall(figHandle,'-property','FontSize');
% set(textHandles,'FontSize', 10, 'fontname', 'Arial')
% 
% % set the label sizes and color
% labelHandles = findall(ax, 'type', 'text', 'handlevisibility', 'off');
% set(labelHandles,'FontSize', 12, 'color', labelColor)
% 
% 
% % separate figure axes:
% xMin = ax.XLim(1);
% yMin = ax.YLim(1);
% ax.XLim(1) = xMin - range(ax.XTick(1:2))/4; % offset for X axis
% ax.YLim(1) = yMin - range(ax.YTick(1:2))/3; % offset for Y axis
% drawnow()
% % set the X axis vertex start to its the orignial point
% origin = get(ax.XAxis.Axle,'VertexData');
% origin(1,1) = xMin;
% set(ax.XAxis.Axle,'VertexData',origin);
% % set the Y axis vertex start to its the orignial point
% origin = get(ax.YAxis.Axle,'VertexData');
% origin(2,1) = yMin;
% set(ax.YAxis.Axle,'VertexData',origin);


end

