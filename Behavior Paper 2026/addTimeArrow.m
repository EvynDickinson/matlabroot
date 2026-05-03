

function addTimeArrow(ax, ArrowColor, y_offset)
% addTimeArrow(ax, foreColor)
%
% PURPOSE
% Adds a 'time' label with a rightward arrow annotation to the bottom
% left of a given axes, indicating the direction of time flow.
%
% INPUTS
%   'ax'  : target axes object
%       handle to the axes on which the annotation will be drawn
%   'foreColor' : color of the text and arrow
%       1-by-3 RGB vector, e.g. [1 1 1] for white
%   'y_offset' : offset from the main axes 
%           default -0.06
%
% OUTPUTS
%   none — draws directly onto ax
%
% EXAMPLE
%   addTimeArrow(ax_heat, [0 0 0])
%
% ES DICKINSON, YALE, 2026

%%

if nargin<3
    y_offset = -0.06;
end

xl = xlim(ax);
yl = ylim(ax);

text_y = yl(1) + y_offset * diff(yl);   % below x-axis

isReversedX = strcmp(ax.XDir, 'reverse');
if isReversedX
    text_x = xl(2) + 0.03 * diff(xl);
    arrow_str = 'time\rightarrow';
else
    text_x = xl(1) + 0.03 * diff(xl);
    arrow_str = 'time\rightarrow';% \Rightarrow is bigger
end

text(ax, text_x, text_y, arrow_str, ...   
    'Color', ArrowColor, 'FontSize', 11, ...
    'FontAngle', 'italic',...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'Clipping', 'off')
end