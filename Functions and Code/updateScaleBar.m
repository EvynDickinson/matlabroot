function updateScaleBar(ax, foreColor, scaleBar_duration, scale_label_str)
% updateScaleBar(ax, foreColor, scaleBar_duration, scale_label_str)
%
% PURPOSE
% Draws a dynamic scale bar on a given axes object to replace the x-axis.
% The bar repositions automatically when called (e.g. on zoom/pan) by
% deleting and redrawing tagged graphics objects. Compatible with both
% standard line plot axes and imagesc axes (reversed y-direction).
%
% INPUTS
%   'ax'                : target axes object
%       handle to the axes on which the scale bar will be drawn
%   'foreColor'         : color of the scale bar and text
%       1-by-3 RGB vector, e.g. [1 1 1] for white
%   'scaleBar_duration' : length of the scale bar in x-axis units
%       scalar, e.g. 100 for 100 mins
%   'scale_label_str'   : label to display below the scale bar
%       string, e.g. '100 min' or '1000 frames'
%
% OUTPUTS
%   none — draws directly onto ax
%
% NOTES
%   Scale bar and label are tagged 'ScaleBar' so they can be found and
%   deleted on redraw without affecting other plot elements. Position is
%   computed as a fraction of the current xlim/ylim, so the bar stays
%   proportionally placed after zooming or panning. Automatically detects
%   reversed y-axis (imagesc) and adjusts bar/label placement accordingly.
%   Attach to an axes listener to trigger on zoom/pan:
%       addlistener(ax, 'XLim', 'PostSet', ...
%           @(~,~) updateScaleBar(ax, foreColor, scaleBar_duration, scale_label_str))
%
% EXAMPLE
%   updateScaleBar(ax4, [1 1 1], 100, '100 min')     % line plot
%   updateScaleBar(ax,  [1 1 1], 1000, '1000 frames') % imagesc
%
% ES DICKINSON, YALE, 2026
%%
delete(findobj(ax, 'Tag', 'ScaleBar'))

offset_percent = 0.005;

xl = xlim(ax);
yl = ylim(ax);

scaleBar_x_start = xl(1) + 0.03 * diff(xl);   % left-aligned
scaleBar_x_end   = scaleBar_x_start + scaleBar_duration;

% handle normal vs reversed y-axis (imagesc sets YDir to 'reverse')
isReversedY = strcmp(ax.YDir, 'reverse');
if isReversedY
    % expand the axes limits downward to make room for the scale bar
    room = offset_percent * abs(diff(yl));
    ylim(ax, [yl(1), yl(2) + room])   % push bottom limit down
    scaleBar_y = yl(2) + 0.04 * room;
    text_y     = scaleBar_y + 0.03 * room;
    textVA     = 'top';
else
    room = offset_percent * abs(diff(yl));
    ylim(ax, [yl(1) - room, yl(2)])
    scaleBar_y = yl(1) - 0.04 * room;
    text_y     = scaleBar_y - 0.03 * room;
    textVA     = 'top';
end

% % debugging code: 
% fprintf('YDir: %s\n', ax.YDir)
% fprintf('yl: [%.2f, %.2f]\n', yl(1), yl(2))
% fprintf('diff(yl): %.2f\n', diff(yl))
% fprintf('abs(diff(yl)): %.2f\n', abs(diff(yl)))
% fprintf('scaleBar_y: %.2f\n', scaleBar_y)
% fprintf('text_y: %.2f\n', text_y)
% fprintf('foreColor: [%.2f, %.2f, %.2f]\n', foreColor(1), foreColor(2), foreColor(3))

line(ax, [scaleBar_x_start, scaleBar_x_end], [scaleBar_y, scaleBar_y], ...
    'Color', foreColor, 'LineWidth', 2, 'Clipping', 'off', 'Tag', 'ScaleBar')

text(ax, scaleBar_x_start, text_y, ...
    scale_label_str, ...
    'Color', foreColor, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', textVA, 'FontSize', 18, ...
    'Clipping', 'off', 'Tag', 'ScaleBar')

end