
function h = plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, FA)
% h = plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, FaceAlpha);
%
% PURPOSE
% Plot a shaded y-error region 
%  ** automatically removes nans for plotting ** 
%
% INPUTS
%    'plot_err' : logical for plot this error region or not [true | false]
%    'x' : vector of x values (like the mean)
%    'y' : vector of y values (like the mean)
%    'y_err' : vector of y-error values (this will be the shaded region around the 'y' vector
%    'kolor' : color to plot the shaded region
%    'fig_type' : '-pdf' or '-png' --> this changes whether the lines are
%           plotted for the shaded region or just the shaded region (since PDFs do
%           not do the semi-transparent regions well, so they get the outer lines for
%           the error region instead
%    'FA' : 0.2 (default shading alpha)
%      
% OUTPUT
%    'h' : handle for the shaded region
%
% ES Dickinson, Yale University 2023

%%
% Defaults
if nargin <= 6
    FA = 0.2; % FaceAlpha
end

% skip plotting anything if plot error is turned off
if ~plot_err
    return
end

% check for nans that would cancel the error...
% set the orientation to be the same :
z = size(x);
if ~(z(2)==1 && z(1)>1)
    x = x';
end
z = size(y);
if ~(z(2)==1 && z(1)>1)
    y = y';
end
z = size(y_err);
if ~(z(2)==1 && z(1)>1)
    y_err = y_err';
end

% remove nan locations across all values if present in any location
% required so that the error fill will plot 
loc = isnan(x) | isnan(y) | isnan(y_err);
x(loc) = [];
y(loc) = [];
y_err(loc) = [];

% plot shaded region
if plot_err && ~strcmpi(fig_type,'-pdf')
    fill_data = error_fill(x, y, y_err);
    h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
    set(h, 'facealpha', FA);
elseif strcmpi(getenv('COMPUTERNAME'),'EvynPC')
    fill_data = error_fill(x, y, y_err);
    h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none','HandleVisibility','off');
    set(h, 'facealpha', FA);
elseif plot_err && strcmpi(fig_type,'-pdf')
    plot(x,y-y_err,'color',kolor, 'linewidth', 0.5,'HandleVisibility','off');
    plot(x,y+y_err,'color',kolor, 'linewidth', 0.5,'HandleVisibility','off');
    h = [];
end


