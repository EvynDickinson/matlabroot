function fig = plotVerticalHist(fig, y, options)
% fig = plotVerticalHist(fig, y, options)
%
% PURPOSE
% Plots a vertical smoothed histogram on the current axes of the provided
% figure, optionally with a filled region under the curve
%
% INPUTS (fig and y required, all others optional)
%   'fig'    : figure handle to plot into
%   'y'      : data vector to compute histogram from
%  INPUT OPTIONS
%   'Filled'       : toggle filled region under the histogram curve
%       true  = fill the area between the curve and the zero axis
%       false = line only
%       (default : true)
%   'Color'        : RGB color vector for the fill and line
%       (default : Color('dodgerblue'))
%   'FaceAlpha'    : face alpha transparency of the filled region
%       (default : 0.3)
%   'LineWidth'    : line width of the histogram curve
%       (default : 1.5)
%   'NumBins'      : number of histogram bins
%       (default : 30)
%   'SmoothWindow' : smoothing window size in number of bins for gaussian smoothing
%       (default : 10)
%
% OUTPUTS
%   'fig'    : figure handle with the vertical histogram plotted on current axes
%
% EXAMPLE
%   fig = figure;
%   y = randn(100,1);
%   fig = plotVerticalHist(fig, y)                             % all defaults
%   fig = plotVerticalHist(fig, y, NumBins=50)                 % more bins
%   fig = plotVerticalHist(fig, y, Filled=false, LineWidth=2)  % line only
%   fig = plotVerticalHist(fig, y, Color=[1 0 0], FaceAlpha=0.5, SmoothWindow=5)
%
% ES DICKINSON, 2026

%% CODE:
arguments
    fig
    y
    options.Filled = true
    options.Color = Color('dodgerblue')
    options.FaceAlpha = 0.3
    options.LineWidth = 1.5
    options.NumBins = 30
    options.SmoothWindow = 10
end
% unpack options (for convenience only)
filled = options.Filled;
kolor  = options.Color;
FA     = options.FaceAlpha;
LW     = options.LineWidth;
nBins  = options.NumBins;
smWin  = options.SmoothWindow;

% PLOTTING
% activate figure
figure(fig);

% histogram distribution
[counts, edges] = histcounts(y, nBins);
bin_centers = edges(1:end-1) + diff(edges)/2;

% smooth the counts
counts_smooth = smoothdata(counts, 'gaussian', smWin);

% plot vertically
if filled
    fill([0; counts_smooth(:); 0], ...
         [bin_centers(1); bin_centers(:); bin_centers(end)], ...
         kolor, 'FaceAlpha', FA, 'EdgeColor', 'none');
end
plot(counts_smooth, bin_centers, 'Color', kolor, 'LineWidth', LW);

hold off;