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

% default parameters:
txtSize = 18; % text size
axisLW = 2; % axis line width
txtFont = 'Arial'; % all fonts
tickMarkSize = 15;    
tickDir = 'in'; 

% SET PLOT BACKGROUND COLOR
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
else % just adjust the line formats but no color
    set(figHandle, 'color', 'w') 
end

% SET AXIS SIZE & COLORS AND LABEL FONTS & SIZE
if nargin >= 3 %subplots
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
        ax.LineWidth = axisLW; %change axes lines to width of 1
        box off
        % set the subplot axes to the selected color scheme
        set(ax, 'color', backColor, 'YColor', labelColor, 'XColor', labelColor)
        % set the font: 
        set(ax, 'FontName', txtFont);

        % set the label sizes and color
        set(gca,'fontsize',tickMarkSize,'FontWeight','normal')
        labelHandles = findall(ax, 'type', 'text', 'handlevisibility', 'off');
        set(labelHandles,'FontSize', txtSize, 'color', labelColor)
    end    
else    % single plot
   % get plot handles
    ax = gca(figHandle);
    try
    ax.LineWidth = axisLW; %change axes lines to preset width
    box off
    catch
    end
    
    % set the subplot axes to the selected color scheme
    set(ax, 'color', backColor, 'YColor', labelColor, 'XColor', labelColor)
    % set the font: 
    set(ax, 'FontName', txtFont);
    
    % set the label sizes and color
    labelHandles = findall(ax, 'type', 'text', 'handlevisibility', 'off');
    set(labelHandles,'FontSize', txtSize, 'color', labelColor)
    set(gca,'fontsize',tickMarkSize,'FontWeight','normal')
    
    
end



