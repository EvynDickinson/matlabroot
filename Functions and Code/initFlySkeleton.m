

function h = initFlySkeleton(ax, kolor, nodes, nodeSize)
% h = initFlySkeleton(ax, kolor, nodes, nodeSize)
%
% PURPOSE:
% creates a graphic handle that can be updated to create
% new 'plotted' fly skeletons without plotting a million new 
% objects -- this increases the speed and efficiency dramatically
% Paired with : updateFlySkeleton.m 
%
% INPUTS
%       'ax' : axes handle
%       'kolor' : color to plot (RBG format)
%       'nodes' : logical if plotting node circles
%       'nodeSize' : scatter point size for the fly skeleton nodes
%
% OUTPUTS
%       'h' : handle to the graphics
%

%% 

if nargin < 3
    nodes = true; % default to ON
end
if nargin < 4
    nodeSize = 50; % default node size
end

% which body points should be connected
skeleton = [1 2; 2 3; 2 4; 2 5];

axes(ax); hold(ax,'on') % activate the given axes

% --- Create edge lines (all 4 in one object)
% --- Create 4 edge line objects
for i = 1:size(skeleton,1)
    h.edges(i) = plot(ax, nan, nan, ...
        'Color', kolor, ...
        'HandleVisibility','off');
end
% this doesn't actually plot anything (since NaNs)
% it just creates basically a graphics template that will
% later be updated

% --- Create node markers
if nodes
    h.nodes = plot(ax, nan, nan, 'o', ...
        'MarkerSize', sqrt(nodeSize), ...
        'MarkerFaceColor', kolor, ...
        'MarkerEdgeColor', kolor, ...
        'HandleVisibility','off');
else
    h.nodes = [];
end

h.skeleton = skeleton;
h.color = kolor;


