


function fig = plotVerticalPDF(fig, y, options)
% fig = plotVerticalPDF(fig, y, options)
%
% PURPOSE
% Plots a vertical kernel density estimate (KDE) as a smoothed PDF on the
% current axes of the provided figure, optionally with a filled region
% under the curve
%
% INPUTS (fig and y required, all others optional)
%   'fig'    : figure handle to plot into
%   'y'      : data vector to compute KDE from
%  INPUT OPTIONS
%   'Filled' : toggle filled region under the KDE curve
%       true  = fill the area between the curve and the zero axis
%       false = line only
%       (default : true)
%   'Color'  : RGB color vector for the fill and line
%       (default : Color('dodgerblue'))
%   'FaceAlpha'     : face alpha transparency of the filled region
%       (default : 0.6)
%   'LineWidth'     : line width of the KDE curve
%       (default : 1.5)
%
% OUTPUTS
%   'fig'    : figure handle with the vertical PDF plotted on current axes
%
% EXAMPLE
%   fig = figure;
%   y = randn(100,1);
%   fig = plotVerticalPDF(fig, y)                        % all defaults
%   fig = plotVerticalPDF(fig, y, LW=2.0)                % thicker line
%   fig = plotVerticalPDF(fig, y, filled=false, LW=1.5)  % line only
%   fig = plotVerticalPDF(fig, y, color=[1 0 0], FA=0.3) % red, transparent
%
% ES DICKINSON, 2026

%% CODE: 
arguments
    fig
    y
    options.Filled  = true
    options.Color = Color('dodgerblue')
    options.FaceAlpha  = 0.6
    options.LineWidth  = 1.5
end

% unpack options (for convience only)
filled = options.Filled;
kolor = options.Color;
FA     = options.FaceAlpha;
LW    = options.LineWidth;

% activate figure
figure(fig);

% density function distribution
[f, xi] = ksdensity(y);     % smoothed density of Y

% plot kernal density estimate
if filled
    fill([0; f(:); 0], [xi(1); xi(:); xi(end)], ...
        kolor, 'FaceAlpha', FA, 'EdgeColor', 'none');
end
plot(f, xi, 'Color', kolor, 'LineWidth', LW);
    


    