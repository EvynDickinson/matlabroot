
function h = plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, FA)
% h = plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, FaceAlpha);
% h = plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type);
%
% ==== INPUTS ====
% plot_err = true | false
% x = vector
% y = vector
% y_err = vector
% kolor = color of line
% fig_type = '-pdf' or '-png'
% FA = 0.2 (default shading alpha)
% automatically removes nans for plotting
%
% ES Dickinson, Yale University 2023

% Defaults
if nargin <= 6
    FA = 0.2; % FaceAlpha
end

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

% remove nans so that the error fill will plot
loc = isnan(x) | isnan(y) | isnan(y_err);
x(loc) = [];
y(loc) = [];
y_err(loc) = [];

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


