
function fig = plotFlySkeleton(fig, x, y, kolor, nodes,nodeSize)
% fig = plotFlySkeleton(fig, x,y,kolor,nodes,nodeSize)
% 
% INPUTS
%       'fig' : figure handle
%       'x' : matrix with the 5 y values of the fly body position
%           rows are frames, columns are body positions
%       'y' : matrix with the 5 y values of the fly body position
%           rows are frames, columns are body positions
%       'kolor' : color to plot (RBG format)
%       'nodes' : logical if plotting node circles
%       'nodeSize' : scatter point size for the fly skeleton nodes
%
% OUTPUTS
%       'fig' : handle to the figure
%
% ES DICKINSON 2025

%%

if nargin < 3
    kolor = 'k';
end

if nargin < 4
    nodes = true;
end

if nargin < 5
    nodeSize = 50;
end

skeleton = [1,2; 2,3; 2,4; 2,5];

% edges
xs = x(skeleton)';
ys = y(skeleton)';
plot(xs, ys, 'Color', kolor, 'HandleVisibility','off')

% nodes
if nodes 
    plot(x, y, 'o', ...
     'MarkerSize', sqrt(nodeSize), ...
     'MarkerFaceColor', kolor, ...
     'MarkerEdgeColor', kolor, ...
     'HandleVisibility','off')
end


%%  Original slower code: 
% if ~exist('kolor', 'var')
%     kolor = 'k';
% end
% skeleton = [1,2; 2,3; 2,4; 2,5];
% 
% figure(fig)
% hold on
% % edges
% for i = 1:size(skeleton,1)
%     plot(x(skeleton(i,:)),y(skeleton(i,:)), 'color', kolor,'HandleVisibility','off')
% end
% % nodes
% if nodes 
%     if ~exist('nodeSize', 'var')
%         nodeSize = 50;
%     end
%     scatter(x, y, nodeSize, kolor, 'filled')
% end
