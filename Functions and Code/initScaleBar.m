
function initScaleBar(ax, foreColor, scaleBar_duration, scale_label_str)
% initScaleBar(ax, foreColor, scaleBar_duration, scale_label_str)
%
% PURPOSE
% Initialises a scale bar on a given axes and attaches a listener so it
% redraws automatically on zoom or pan.
%
% INPUTS
%   'ax'                : target axes handle
%   'foreColor'         : 1-by-3 RGB color vector
%   'scaleBar_duration' : length of scale bar in x-axis units
%   'scale_label_str'   : label string, e.g. '100 min'
%
% EXAMPLE
%   initScaleBar(ax4, foreColor, 100, '100 min')
%
% ES DICKINSON, YALE, 2026

    axes(ax)
    set(ax, 'XColor', 'none')
    updateScaleBar(ax, foreColor, scaleBar_duration, scale_label_str)
    addlistener(ax, 'XLim', 'PostSet', ...
        @(~,~) updateScaleBar(ax, foreColor, scaleBar_duration, scale_label_str))
end