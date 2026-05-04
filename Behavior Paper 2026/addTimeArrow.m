

function addTimeArrow(ax, ArrowColor, y_offset, x_offset_norm, x_offset_reversed)
% addTimeArrow(ax, ArrowColor, y_offset, x_offset_norm, x_offset_reversed)
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
%   'x_offset_norm' : offset from the main x axes for normally directed  axes
%           default 0.01
%   'x_offset_norm' : offset from the main x axes for normally directed  axes
%           default  0.03
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

if ~exist('x_offset_norm', 'var')
    x_offset_norm = 0.01;
end

if ~exist('x_offset_reversed', 'var')
    x_offset_reversed = 0.03;
end


xl = xlim(ax);
yl = ylim(ax);

text_y = yl(1) + y_offset * diff(yl);   % below x-axis

isReversedX = strcmp(ax.XDir, 'reverse');
if isReversedX
    text_x = xl(2) + x_offset_reversed * diff(xl);
    arrow_str = 'time\rightarrow';
else
    text_x = xl(1) + x_offset_norm * diff(xl);
    arrow_str = 'time\rightarrow';% \Rightarrow is bigger
end

text(ax, text_x, text_y, arrow_str, ...   
    'Color', ArrowColor, 'FontSize', 11, ...
    'FontAngle', 'italic',...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'Clipping', 'off')

