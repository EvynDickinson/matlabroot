
function figHandle = formatFig(figHandle, BandWoption, sbplots, subplotInd)
% figHandle = formatFig(figHandle, BandWoption, sbplots, subplotInd);
% 
% Beautifies a figure by setting basic parameters to herein defined defaults
% Also determines colors of the axes and background based on whether the
% desired figure is black or white background based
% 
% INPUTS
%    'figHandle' :  name of handle for the figure
%    'BandWoption' : logical for black background
%           true = black background
%           false = white background
%    'sbplots' : number of subplots *[row,columns]
%    'subplotInd' : subplot indexes if not 1 space per subplot
%           e.g. subplotInd(1).idx = [1,3,4]; 
%
% OUTPUT
%    'figHandle' : handle for the figure
% 
% EXAMPLE
%    fig = formatFig(fig,true);
% 
% ES Dickinson, 2020

%% 
% default parameters:
txtSize = 16; % text size
axisLW = 2; % axis line width
txtFont = 'Arial'; % font for all text 
tickMarkSize = 14;  % length of tick marks
tickDir = 'out'; % direction of tick marks

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
else % just adjust the line formats but not color
    set(figHandle, 'color', 'w') 
end

% SET AXIS SIZE & COLORS AND LABEL FONTS & SIZE
if nargin >= 3 % e.g. there are subplots
    %there are multiple plots
    ncol = sbplots(2);  %number of columns
    nrow = sbplots(1);  %number of rows
    nplots = length(findobj(figHandle,'type','axes'));
        % set properties for each subplot:
    for n = 1:nplots
        % find name of the field with the indexes for the subplots
        if nargin==4 % does not assume even subplot sizing
            sb_field = fieldnames(subplotInd);
            sb_field = sb_field{1};
            ii = subplotInd(n).(sb_field);
        else % assumes each subplot is a only 1 unit size
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
        set(ax,'tickdir',tickDir)
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
    set(ax,'tickdir',tickDir)
    
end



