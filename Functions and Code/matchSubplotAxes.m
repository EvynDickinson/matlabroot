function fig = matchSubplotAxes(fig, options)
% fig = matchSubplotAxes(fig, options)
%
% PURPOSE
% Finds all subplots in a figure and matches the axis scale across all
% subplots to either the largest or smallest range found across all axes
%
% INPUTS (fig required, all others optional)
%   'fig'    : figure handle containing subplots to match
%  INPUT OPTIONS
%   'X'         : match x axes across all subplots
%       true  = match x axis scale
%       false = leave x axis scale unchanged
%       (default : true)
%   'Y'         : match y axes across all subplots
%       true  = match y axis scale
%       false = leave y axis scale unchanged
%       (default : true)
%   'Method'    : which range to match across subplots
%       'largest'  = use the largest range found across all subplots
%       'smallest' = use the smallest range found across all subplots
%       (default : 'largest')
%
% OUTPUTS
%   'fig'    : figure handle with matched axis scales across all subplots
%
% EXAMPLE
%   fig = matchSubplotAxes(fig)                              % match both axes to largest range
%   fig = matchSubplotAxes(fig, X=false)                     % y axis only
%   fig = matchSubplotAxes(fig, Y=false)                     % x axis only
%   fig = matchSubplotAxes(fig, Method='smallest')           % match to smallest range
%   fig = matchSubplotAxes(fig, X=true, Y=false, Method='largest') % x only, largest
%
% ES DICKINSON, 2026

%% CODE:
arguments
    fig
    options.X  = true
    options.Y = true
    options.Method = 'largest'  % 'largest' or 'smallest'
end
% unpack options
matchX  = options.X;
matchY  = options.Y;
method  = options.Method;

% activate figure
figure(fig);

% find axes collection - excludes legends, colorbars, and other non-plot axes
allAx = findobj(fig, 'Type', 'axes');
isPlotAx = arrayfun(@(a) strcmp(get(a,'Tag'), ''), allAx);
allAx = allAx(isPlotAx);
if isempty(allAx)
    warning('matchSubplotAxes: no axes found in figure');
    return
end

% collect x and y limits across all subplots
xLims = zeros(numel(allAx), 2);
yLims = zeros(numel(allAx), 2);
for i = 1:numel(allAx)
    xLims(i,:) = xlim(allAx(i));
    yLims(i,:) = ylim(allAx(i));
end

% compute x limits to apply
if matchX
    xRanges = diff(xLims, 1, 2);
    if strcmpi(method, 'largest')
        [~, idx] = max(xRanges);
    else
        [~, idx] = min(xRanges);
    end
    xLimNew = xLims(idx,:);
end

% compute y limits to apply
if matchY
    yRanges = diff(yLims, 1, 2);
    if strcmpi(method, 'largest')
        [~, idx] = max(yRanges);
    else
        [~, idx] = min(yRanges);
    end
    yLimNew = yLims(idx,:);
end

% apply matched limits to all subplots
for i = 1:numel(allAx)
    if matchX
        xlim(allAx(i), xLimNew);
    end
    if matchY
        ylim(allAx(i), yLimNew); 
    end
end