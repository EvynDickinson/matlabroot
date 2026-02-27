
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

if ~exist('nodeSize', 'var')
    nodeSize = 50;
end
 
skeleton = [1,2; 2,3; 2,4; 2,5]; % specific links for lines between body parts
nFlies = size(x, 1); % number of fly skeletons to plot

% Vectorized edge plotting
% Create all edge coordinates at once
x_edges = nan(3 * size(skeleton, 1), nFlies);
y_edges = nan(3 * size(skeleton, 1), nFlies);

for i = 1:size(skeleton, 1)
    row_start = (i-1)*3 + 1;
    % extract X edges
    x_edges(row_start, :) = x(:, skeleton(i, 1))';
    x_edges(row_start+1, :) = x(:, skeleton(i, 2))';
    % extract Y edges
    y_edges(row_start, :) = y(:, skeleton(i, 1))';
    y_edges(row_start+1, :) = y(:, skeleton(i, 2))';
end

% Plot all edges in one call
plot(x_edges(:), y_edges(:), 'Color', kolor, 'HandleVisibility', 'off');

% Plot all nodes
if nodes
    scatter(x(:), y(:), nodeSize, kolor, 'filled', 'HandleVisibility', 'off');
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
