

function updateFlySkeleton(h, x, y, kolor)
% updateFlySkeleton(h, x, y, kolor)
%
% PURPOSE: 
% update the fly body skeleton graphics object 
% that was initiated with : initFlySkeleton.m
% 
%
% INPUTS
%       'h' : graphics object handle
%       'x' : matrix with the 5 y values of the fly body position
%           rows are frames, columns are body positions
%       'y' : matrix with the 5 y values of the fly body position
%           rows are frames, columns are body positions
%       'kolor' : color to plot (RBG format)
%
% ES DICKINSON 2025

%%

% Update color if needed
if ~isequal(kolor, h.color)
    set(h.edges,'Color',kolor)
    if ~isempty(h.nodes)
        set(h.nodes,'MarkerFaceColor',kolor,...
                    'MarkerEdgeColor',kolor)
    end
    h.color = kolor;
end

% Update edges
for i = 1:size(h.skeleton,1)
    idx = h.skeleton(i,:); % each line updated alone
    set(h.edges(i), ...
        'XData', x(idx), ...
        'YData', y(idx));
end

% Update nodes if they are set to be plotted
if ~isempty(h.nodes)
    set(h.nodes,'XData',x,'YData',y)
end

