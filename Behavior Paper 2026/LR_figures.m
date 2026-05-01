
blkbgd = false;
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

initial_vars = add_var(initial_vars, 'saveDir');
saveDir = createFolder([figDir, 'Paper Figures/']);

% 
% field_names = {'escape jump', 'escape ring', 'food quadrant', 'fly on food', 'courtship', 'sleep'};
% kolors = {'WongOrange', 'WongRed','WongBlue','WongLightBlue', 'WongPink', 'WongGreen'}; % colors for the diff behaviors

%% FIGURE: time course of all the different behaviors ordered by precidence
clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd); %get background colors

plot_err = true;
LW = 1.5;
sSpan = 360; % 2 minute smoothing 

% 
r = 7;
c = 1;
s(1).idx = 1;
s(2).idx = 2:3;
s(3).idx = 4:5;
s(4).idx = 6:7;

g = 1; % group index

x = grouped(g).time;   
kolor = foreColor;

fig = getfig('', 1, [ 648 680]);

%temp
ax1 = subplot(r,c,s(1).idx); hold on
    y = grouped(g).temp;
    plot(x,y,'LineWidth',LW,'Color',kolor)
    ylabel({'temp'; '(\circC)'})

% escape ring
ax2 = subplot(r,c,s(2).idx); hold on
    y = smooth(grouped(g).ring.avg, sSpan, 'moving');
    plot(x,y,'LineWidth',LW,'Color',kolor)
    ylabel({'escape ring'; '(% flies)'})

% food
ax3 = subplot(r,c,s(3).idx); hold on
    y = smooth(grouped(g).innerquad.food.avg, sSpan, 'moving');
    plot(x,y,'LineWidth',LW,'Color',kolor)
    ylabel({'food region'; '(% flies)'})

% sleep
ax4 = subplot(r,c,s(4).idx); hold on
    y = smooth(grouped(g).sleep.avg, sSpan, 'moving');
    plot(x,y,'LineWidth',LW,'Color',kolor)  
    ylabel({'sleeping'; '(% flies)'})

linkaxes([ax1, ax2, ax3, ax4], 'x') % link axes
xlim([100, 650])

% ===========  formatting ===========  
formatFig(fig, blkbgd, [r,c], s);

% set y lims & labels
subplot(r,c,s(1).idx) % temperature
    ylim([14.7, 35.3])
    set(ax1, 'YTick', [15,35])
 
subplot(r,c,s(2).idx) % escape ring
    ylim([0, 50])
    set(ax2, 'YTick', [0, 25, 50])    

subplot(r,c,s(3).idx) % at food
    ylim([0, 80])
    set(ax3, 'YTick', 0:40:80)  

subplot(r,c,s(4).idx) % at food
    ylim([-3, 43])
    set(ax4, 'YTick', 0:20:40)  
    
% Hide x-axis on all subplots
for ii = 1:4
    subplot(r,c,s(ii).idx)
    set(gca, 'XColor', 'none')  % hide x-axis line, ticks, labels
end

% Add scale bar to the bottom subplot (ax4)
axes(ax4)
xl = xlim;
yl = ylim;

% --- Scale bar parameters ---
scaleBar_duration = 100;        % length in x-axis units (e.g., min)
scaleBar_x_start = xl(1) + 0.03 * diff(xl);   % left-aligned with small margin
scaleBar_x_end   = scaleBar_x_start + scaleBar_duration;
scaleBar_y = yl(1) - 0.08 * diff(yl);       % just below the data

% Draw the scale bar line
annotation_ax = ax4;
line(annotation_ax, ...
    [scaleBar_x_start, scaleBar_x_end], ...
    [scaleBar_y, scaleBar_y], ...
    'Color', foreColor, 'LineWidth', 2, 'Clipping', 'off')

% Add label centered above the bar
scale_label_str = sprintf('%d min', scaleBar_duration);
text(ax4, scaleBar_x_start, scaleBar_y - 0.04*diff(yl), ...
    scale_label_str, ...
    'Color', foreColor, 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'FontSize', 18, ...
    'Clipping', 'off', 'Tag', 'ScaleBar')
addlistener(ax4, 'XLim', 'PostSet', @(~,~) updateScaleBar(ax4, foreColor, scaleBar_duration, scale_label_str));


% save figure: 
save_figure(fig, [saveDir, 'Fig 1 timecourse of behavior ', grouped(g).name]);

%% FIGURE: temperature tuning panel with all behaviors
clearvars('-except',initial_vars{:})

exp = 2; % which exp to plot

autoLim = false;
manual_xlims = [13, 37];
% manual_xlims = [15, 27];
% manual_xlims = [13, 25];

plot_err = true; % plot SEM
plot_high_null = true; % plot the low or high null occupancy for empty trials
foreColor = formattingColors(blkbgd,true); % get foreground color for current background configuration

narrow_fig = false;
narrow_fig_size = [676 680];

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end,quad_regions] = PlotParamSelection(false,false,true);
nTypes = length(pName);
% get plotting colors for the paper scheme:
colors = getPaperColors(pName);
plot_err = true;
  
% set region selection for trials that have subregions to select between
null_quad = 'high'; % highest occupancy quadrant
food_quad = 'food'; % quadrant with food

% Plotting Parameters
LW = 0.5; % line width
hLW = 1.5; % highlight line width for empty trials (to indicate that it is plotting highest occupancy region, not the food quadrant) 
buff = 0.5; % linewidth buffer to match the non null to null width...
xlimits = []; % initialize empty x for finding uniform x axis scaling later

% FIGURE:
r = 1;
c = 2; % heating and cooling
fig = getfig('', true);
ax_cool = subplot(r, c, 1); hold on; title('cooling', 'Color', foreColor)
ax_heat = subplot(r, c, 2); hold on; title('heating', 'Color', foreColor)

for roi = 1:nTypes
    kolor = colors(roi,:);
    % determine if this is a food or null trial
    if data(exp).emptytrial
        subfield = null_quad;
        hColor = foreColor;
        highlight = false; % change to true if desired outlines
    else
        subfield = food_quad;
        highlight = false;
    end

    % pull the subfield data structure
    if quad_regions(roi)
        yy = grouped(exp).(pName{roi}).(subfield);
    else
        yy = grouped(exp).(pName{roi});
    end

    x = yy.temps;
    YH = yy.increasing.raw  * scaler(roi);
    YC = yy.decreasing.raw  * scaler(roi);
    y_H = mean(YH, 2, 'omitnan');
    y_C = mean(YC, 2, 'omitnan');
    y_err_H = std(YH, 0, 2, 'omitnan') ./ sqrt(num.trial(exp));
    y_err_C = std(YC, 0, 2, 'omitnan') ./ sqrt(num.trial(exp));

    switch dType
        case 1 % single trial lines
            for trial = 1:num.trial(exp)
                y_H_trial = YH(:, trial);
                y_C_trial = YC(:, trial);
                if highlight
                    plot(ax_heat, x, y_H_trial, 'LineWidth', LW+hLW, 'Color', hColor)
                    plot(ax_cool, x, y_C_trial, 'LineWidth', LW+hLW, 'Color', hColor)
                end
                plot(ax_heat, x, y_H_trial, 'Color', kolor, 'LineWidth', LW+buff)
                plot(ax_cool, x, y_C_trial, 'Color', kolor, 'LineWidth', LW+buff)
            end

        case {2, 3} % avg lines — overlay all metrics on same axes
            if highlight
                plot(ax_heat, x, y_H, 'LineWidth', LW+hLW, 'Color', hColor)
                plot(ax_cool, x, y_C, 'LineWidth', LW+hLW, 'Color', hColor)
            end
            plot(ax_heat, x, y_H, 'Color', kolor, 'LineWidth', LW+buff)
            plot(ax_cool, x, y_C, 'Color', kolor, 'LineWidth', LW+buff)
            axes(ax_heat)
            plot_error_fills(plot_err, x, y_H, y_err_H, kolor, fig_type, 0.35)
            axes(ax_cool)
            plot_error_fills(plot_err, x, y_C, y_err_C, kolor, fig_type, 0.35)
    end
    xlimits = [xlimits, xlim(ax_heat), xlim(ax_cool)];
end

% FORMATTING AND LABELS
xlimits = [min(xlimits), max(xlimits)];
formatFig(fig, blkbgd, [r, c]);
matchSubplotAxes(fig, "Y", true);
set(ax_cool, 'XDir', 'reverse')
set(ax_heat, 'XDir', 'normal')
for ax = [ax_heat, ax_cool]
    axes(ax)
    xlabel(ax, 'temp (\circC)')
    if autoLim
        xlim(xlimits)
    else
        xlim(manual_xlims)
    end
    ylabel('% flies')
end
l = legend(y_lab, 'Box','off', 'Location', 'northwest');

% plot a line for the 'safe' vs 'threat' zones: 
y = rangeLine(fig, 0, false);
xrange = xlim;
x_less = [xrange(1), 25];
x_more = [25, xrange(2)];
axes(ax_heat)
    plot(x_more, [y,y], 'Color', Color('darkred'), 'LineWidth',3, 'HandleVisibility','off') % threat
    plot(x_less, [y,y], 'Color', Color('grey'), 'LineWidth',3, 'HandleVisibility','off') % safe
axes(ax_cool)
    plot(x_less, [y,y], 'Color', Color('darkred'), 'LineWidth',3, 'HandleVisibility','off') % threat
    plot(x_more, [y,y], 'Color', Color('grey'), 'LineWidth',3, 'HandleVisibility','off') % safe


% set the y limits manually
switch grouped(exp).name
    case 'Berlin LTS 15-35 no food'
        ylims = [-1, 50];
        y_ticks = 0:10:50;
        l.Position = [0.34461 0.8549 0.18891 0.095588];
    case 'Berlin LTS 15-35 caviar plate 1'
        ylims = [-1, 80];
        y_ticks = 0:20:80;
end
axes(ax_heat)
ylim(ylims)
set(ax_heat, 'YTick', y_ticks)
addTimeArrow(ax_heat, foreColor);
axes(ax_cool)
ylim(ylims)
set(ax_cool, 'YTick', y_ticks)
addTimeArrow(ax_cool, foreColor);


% save figure: 
save_figure(fig, [saveDir, 'Fig 1 temp tuning curve', grouped(exp).name]);


%% FIGURE: normalized temperature tuning panel with all behaviors
clearvars('-except',initial_vars{:})

exp = 1; % which exp to plot

autoLim = false;
manual_xlims = [13, 37];
% manual_xlims = [15, 27];
% manual_xlims = [13, 25];

plot_err = true; % plot SEM
plot_high_null = true; % plot the low or high null occupancy for empty trials
foreColor = formattingColors(blkbgd,true); % get foreground color for current background configuration

narrow_fig = false;
narrow_fig_size = [676 680];

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end,quad_regions] = PlotParamSelection(false,false,true);
nTypes = length(pName);
% get plotting colors for the paper scheme:
colors = getPaperColors(pName);

switch questdlg('Plot error?','','True','False', 'Cancel','True')
    case 'True'
        plot_err = true;
    case 'False'
        plot_err = false;
    case 'Cancel'
        return
    case ''
        return
end

% set region selection for trials that have subregions to select between
null_quad = 'high'; % highest occupancy quadrant
food_quad = 'food'; % quadrant with food

% Plotting Parameters
LW = 0.5; % line width
hLW = 1.5; % highlight line width for empty trials (to indicate that it is plotting highest occupancy region, not the food quadrant) 
buff = 0.5; % linewidth buffer to match the non null to null width...
xlimits = []; % initialize empty x for finding uniform x axis scaling later

% FIGURE:
r = 1;
c = 2; % heating and cooling
fig = getfig('', true);
ax_cool = subplot(r, c, 1); hold on; title('cooling', 'Color', foreColor)
ax_heat = subplot(r, c, 2); hold on; title('heating', 'Color', foreColor)

for roi = 1:nTypes
    kolor = colors(roi,:);
    % determine if this is a food or null trial
    if data(exp).emptytrial
        subfield = null_quad;
        hColor = foreColor;
        highlight = true;
    else
        subfield = food_quad;
        highlight = false;
    end

    % pull the subfield data structure
    if quad_regions(roi)
        yy = grouped(exp).(pName{roi}).(subfield);
    else
        yy = grouped(exp).(pName{roi});
    end

    x = yy.temps;
    YH = yy.increasing.raw  * scaler(roi);
    YC = yy.decreasing.raw  * scaler(roi);
    y_H = mean(YH, 2, 'omitnan');
    y_C = mean(YC, 2, 'omitnan');
    y_err_H = std(YH, 0, 2, 'omitnan') ./ sqrt(num.trial(exp));
    y_err_C = std(YC, 0, 2, 'omitnan') ./ sqrt(num.trial(exp));

    switch dType
        case 1 % single trial lines
            for trial = 1:num.trial(exp)
                y_H_trial = YH(:, trial);
                y_C_trial = YC(:, trial);
                if highlight
                    plot(ax_heat, x, y_H_trial, 'LineWidth', LW+hLW, 'Color', hColor)
                    plot(ax_cool, x, y_C_trial, 'LineWidth', LW+hLW, 'Color', hColor)
                end
                plot(ax_heat, x, y_H_trial, 'Color', kolor, 'LineWidth', LW+buff)
                plot(ax_cool, x, y_C_trial, 'Color', kolor, 'LineWidth', LW+buff)
            end

        case {2, 3} % avg lines — overlay all metrics on same axes
            if highlight
                plot(ax_heat, x, y_H, 'LineWidth', LW+hLW, 'Color', hColor)
                plot(ax_cool, x, y_C, 'LineWidth', LW+hLW, 'Color', hColor)
            end
            plot(ax_heat, x, y_H, 'Color', kolor, 'LineWidth', LW+buff)
            plot(ax_cool, x, y_C, 'Color', kolor, 'LineWidth', LW+buff)
            axes(ax_heat)
            plot_error_fills(plot_err, x, y_H, y_err_H, kolor, fig_type, 0.35)
            axes(ax_cool)
            plot_error_fills(plot_err, x, y_C, y_err_C, kolor, fig_type, 0.35)
    end
    xlimits = [xlimits, xlim(ax_heat), xlim(ax_cool)];
end


% FORMATTING AND LABELS
xlimits = [min(xlimits), max(xlimits)];
formatFig(fig, blkbgd, [r, c]);
matchSubplotAxes(fig, "Y", true);
set(ax_cool, 'XDir', 'reverse')
set(ax_heat, 'XDir', 'normal')
for ax = [ax_heat, ax_cool]
    axes(ax)
    xlabel(ax, 'temp (\circC)')
    if autoLim
        xlim(xlimits)
    else
        xlim(manual_xlims)
    end
    ylabel('% flies')
    addTimeArrow(ax, foreColor);
end
legend(y_lab, 'Box','off', 'Location', 'northwest')
axes(ax_cool)


% plot a line for the 'safe' vs 'threat' zones: 
y = rangeLine(fig, 3, true);
xrange = xlim;
x_less = [xrange(1), 25];
x_more = [25, xrange(2)];
axes(ax_heat)
    plot(x_more, [y,y], 'Color', Color('darkred'), 'LineWidth',3, 'HandleVisibility','off') % threat
    plot(x_less, [y,y], 'Color', Color('grey'), 'LineWidth',3, 'HandleVisibility','off') % safe
axes(ax_cool)
    plot(x_less, [y,y], 'Color', Color('darkred'), 'LineWidth',3, 'HandleVisibility','off') % threat
    plot(x_more, [y,y], 'Color', Color('grey'), 'LineWidth',3, 'HandleVisibility','off') % safe


% save figure: 
save_figure(fig, [saveDir, 'Fig 1 temp tuning curve', grouped(exp).name]);






