
blkbgd = false;

initial_vars = add_var(initial_vars, 'paper_figs');
paper_figs = createFolder([saveDir, 'Paper Figures/']);

% 
% field_names = {'escape jump', 'escape ring', 'food quadrant', 'fly on food', 'courtship', 'sleep'};
% kolors = {'WongOrange', 'WongRed','WongBlue','WongLightBlue', 'WongPink', 'WongGreen'}; % colors for the diff behaviors

%% FIGURE: time course of all the different behaviors ordered by precidence [single exp type]
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
save_figure(fig, [paper_figs, 'Fig 1 timecourse of behavior ', grouped(g).name figcolor(blkbgd)]);

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
save_figure(fig, [paper_figs, 'Fig 1 temp tuning curve', grouped(exp).name]);


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
save_figure(fig, [paper_figs, 'Fig 1 temp tuning curve', grouped(exp).name]);



%% FIGURE: Plot multiple tuning curves -- select the parameters
clearvars('-except',initial_vars{:})

% blkbgd = false;
% fig_type = '-pdf';

autoLim = false;
switch expGroup
    case  'Berlin LTS 15-35 caviar vs empty'
        manual_xlims = [13, 37];
end
% manual_xlims = [13, 37];
% manual_xlims = [15, 27];
% manual_xlims = [13, 25];

plot_err = true; % plot SEM
plot_high_null = true; % plot the low or high null occupancy for empty trials
foreColor = formattingColors(blkbgd,true); % get foreground color for current background configuration

narrow_fig = false;
narrow_fig_size = [676 680];

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end,quad_regions] = PlotParamSelection(true,false,true);
nPlots = length(nullD); % how many parameters will be plotted?
if nPlots == 0
    return
end

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

% set figure folder
fig_dir = createFolder([paper_figs, dir_end]);

% set region selection for trials that have subregions to select between
null_quad = 'high'; % highest occupancy quadrant
food_quad = 'food'; % quadrant with food

% Plotting Parameters
LW = 1; % line width
hLW = 1.5; % highlight line width for empty trials (to indicate that it is plotting highest occupancy region, not the food quadrant) 
buff = 0.5; % linewidth buffer to match the non null to null width...
dataString = cell([1,num.exp]); % for adding a legend -- to fill with group names
r = 1; % rows
c = nPlots; % columns
xlimits = []; % initialize empty x for finding uniform x axis scaling later

% FIGURE:
fig = getfig('',0);
for i = 1:num.exp
    kolor = grouped(i).color;
    dataString{i} = grouped(i).name; % save name of the group in plot order
    for roi = 1:nPlots % plot each parameter seperately
        subplot(r,c,roi); hold on

        % determine if this is a food or null trial
        if data(i).emptytrial 
            subfield = null_quad; 
            hColor = foreColor;
            highlight = false;
        else 
            subfield = food_quad;
            highlight = false;
        end
        % pull the subfield data structure from the grouped data set
        if quad_regions(roi) % sub regions (requires '.food' or '.low' extension etc)
            yy = grouped(i).(pName{roi}).(subfield);
        else
            yy = grouped(i).(pName{roi}); % no subregions in the metric (e.g., ring)
        end

         switch dType
             case 1 % single trial lines
                for trial = 1:num.trial(i)
                    x = yy.temps;
                    rawY = [grouped(i).(pName{roi}).increasing.raw(:,trial),grouped(i).(pName{roi}).decreasing.raw(:,trial)];
                    y = mean(rawY,2,'omitnan')*scaler(roi);
                    if highlight
                        plot(x,y,'LineWidth',LW+hLW,'Color',hColor)
                    end
                    plot(x,y,'color',kolor,'linewidth',LW + buff)
                end
    
             case 2 % avg lines (combined heating and cooling)
                x = yy.temps;
                rawY = [yy.increasing.raw,yy.decreasing.raw];
                y = mean(rawY,2,'omitnan')*scaler(roi);
                y_err = (std(rawY,0,2,'omitnan')*scaler(roi))./sqrt(num.trial(i));
                if highlight
                    plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle',ls)
                end
                plot(x,y,'color',kolor,'linewidth',LW + buff,'LineStyle',ls)
                plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    
             case 3 % separated heating and cooling
                x = yy.temps;
                YC = yy.decreasing.raw;
                YH = yy.increasing.raw;
                % cooling
                y = mean(YC,2,'omitnan')*scaler(roi);
                y_err = (std(YC,0,2,'omitnan')*scaler(roi))./sqrt(num.trial(i));
                plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
                if highlight
                    plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle','--')
                end
                plot(x,y,'color',kolor,'linewidth',LW + buff,'linestyle', '--')
                % heating
                y = mean(YH,2,'omitnan')*scaler(roi);
                y_err = (std(YH,0,2,'omitnan')*scaler(roi))./sqrt(num.trial(i));
                plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
                if highlight
                    plot(x,y,'LineWidth',LW+hLW,'Color',hColor,'LineStyle','-')
                end
                plot(x,y,'color',kolor,'linewidth',LW + buff,'linestyle', '-')          
         end
         xlimits = [xlimits, xlim];
    end
end

% find the xlimits
xlimits = [min(xlimits), max(xlimits)];

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c]);
for roi = 1:nPlots 
    subplot(r,c,roi); hold on
    h_line(nullD(roi),foreColor,'--',1) % plot the null distribution value if it exists...
    ylabel(y_lab{roi})
    xlabel('temperature (\circC)')
    if autoLim
        xlim(xlimits)
    else 
        xlim(manual_xlims)
    end
end

% adjust the figure size if selected
if narrow_fig
    cur_pos = fig.Position;
    new_pos = [cur_pos(1:2), narrow_fig_size];
    fig.Position = new_pos;
end

% save name of parameters being plotted as well: 
param_names = pName{1};
for i = 2:length(pName)
    param_names = [param_names, '_',  pName{i}];
end

set(findall(fig, 'Type', 'axes'), 'FontSize', 22)        % tick labels
set(findall(fig, 'Type', 'text'), 'FontSize', 24)        % all text/titles/labels

save_figure(fig,[paper_figs 'Multimetric tuning curves ' param_names],fig_type);

%% FIGURE: [temp rate comparison] flies in location at coldest points in the temp ramp
% includes static data for that type

clearvars('-except',initial_vars{:})
[foreColor,~] = formattingColors(blkbgd);
[title_str, pName,y_dir,y_lab, nullD, scaler, dType, dir_end, ext] = PlotParamSelection(false);

% load static temp hold data
temp = load([baseFolder, 'Grouped Data Structures/Berlin Temp Holds Caviar/' pName ' data.mat']);
static_data = temp.y; clear temp;
exp_loc = find(strcmp({static_data(:).name}, 'Berlin hold 17 caviar'));
if ext
    static_y = static_data(exp_loc).(pName).food.all(1:54000,:); % first 5 hours of data avg.
else
    static_y = static_data(exp_loc).(pName).all(1:54000,:); % first 5 hours of data avg.
end

% trial averages for static condition
static_avg = mean(static_y, 1, 'omitnan')';   % trials x 1

% extract data across structures
temp_buff = 1.5;
plotData = [];
[temp_rate, temp_diff] = deal(nan(1, num.exp));
for ii = 1:num.exp
    exp = expOrder(ii);
    tp = getTempTurnPoints(data(exp).T.TempProtocol{1});
    if ext
        temps = grouped(exp).(pName).food.temps;
        roi = (temps >= tp.threshLow) & (temps <= tp.threshLow + temp_buff);
        y = grouped(exp).(pName).decreasing.food.raw(roi,:);
    else
        temps = grouped(exp).(pName).temps;
        roi = (temps >= tp.threshLow) & (temps <= tp.threshLow + temp_buff);
        y = grouped(exp).(pName).decreasing.raw(roi,:);
    end
    y_avg = mean(y, 1, 'omitnan');
    plotData = autoCat(plotData, y_avg', false);
    temp_rate(ii) = abs(max(tp.rates));
    temp_diff(ii) = tp.threshHigh - tp.threshLow - 1;
end

% append static data as additional column
plotData = autoCat(plotData, static_avg, false);
all_rates = [temp_rate, 0];   % 0 = static condition

% --- x axis layout with break ---
% place static group at a visually spaced x location
x_break  = min(temp_rate) - 0.25;   % where the axis break sits
x_static = x_break - 0.25;          % x position of static group on plot

% map real temp_rates to plot x positions
x_pos  = temp_rate;           % dynamic groups plot at their real rate
x_pos_static = x_static;           % static plots at artificial x

% -------------------------------------------------------------------------------
% FIGURE
buff = 0.01;
SZ   = 100;
fig  = getfig('', false, [463 660]);
hold on

% dynamic groups
for ii = 1:num.exp
    kolor = grouped(expOrder(ii)).color;
    y = rmnan(plotData(:, ii));
    x = shuffle_data(temp_rate(ii) + linspace(-buff, buff, length(y)));
    scatter(x, y, SZ, kolor, 'filled', ...
        'XJitter', 'density', 'XJitterWidth', 0.05, 'MarkerFaceAlpha', 0.7)
    m   = mean(y, 'omitnan');
    sem = std(y, 0, 'omitnan') ./ sqrt(length(y));
    errorbar(temp_rate(ii), m, sem, ...
        'Color', foreColor, 'LineWidth', 2, 'CapSize', 6)
    scatter(temp_rate(ii), m, SZ+50, foreColor, 'filled', 'square')
end

% static group
y_s = rmnan(static_avg);
x_s = shuffle_data(x_pos_static + linspace(-buff, buff, length(y_s)));
scatter(x_s, y_s, SZ, foreColor, 'filled', ...
    'XJitter', 'density', 'XJitterWidth', 0.05, 'MarkerFaceAlpha', 0.7)
m_s   = mean(y_s, 'omitnan');
sem_s = std(y_s, 0, 'omitnan') ./ sqrt(length(y_s));
errorbar(x_pos_static, m_s, sem_s, ...
    'Color', foreColor, 'LineWidth', 2, 'CapSize', 6)
scatter(x_pos_static, m_s, SZ+50, foreColor, 'filled', 'square')

% --- axis break: railroad double-hash on x axis ---
hash_h   = 0.06 * diff(ylim);   % height of hash marks
hash_w   = 0.07;                 % horizontal width of each hash
hash_gap = 0.05;                 % gap between the two lines

for hh = [-1, 1]
    x_center = x_break + hh * hash_gap/2;
    plot([x_center - hash_w, x_center + hash_w], ...
         [0 - hash_h/2, 0 + hash_h/2], ...
        'Color', foreColor, 'LineWidth', 1.5, 'Clipping', 'off')
end

% formatting
formatFig(fig, blkbgd);
ylabel(y_lab)
xlabel('temp rate (\circC/min)')
ylim([0, max(ylim)]);
% custom x ticks: real rates + static label
x_tick_pos    = [x_pos_static, sort(temp_rate)];
x_tick_labels = ['0', arrayfun(@num2str, sort(temp_rate), 'UniformOutput', false)];
set(gca, 'XTick', x_tick_pos, 'XTickLabel', x_tick_labels, 'XTickLabelRotation', 90)
set(findall(fig, 'Type', 'axes'), 'FontSize', 20)
set(findall(fig, 'Type', 'text'), 'FontSize', 20)
xlim([x_pos_static-0.25, max(temp_rate)+0.25])

% save_figure(fig, [paper_figs, 'temp rate vs ' pName ' scatter at coldest point'], '-pdf', false, false);


% % -------------------------------------------------------------------------------
% % STATS — anova1 with Bonferroni correction, including static group
% datastats.all = [];
% datastats.id  = [];
% 
% for ii = 1:num.exp
%     y = rmnan(plotData(:, ii));
%     datastats.all = [datastats.all; y];
%     datastats.id  = [datastats.id; repmat(temp_rate(ii), length(y), 1)];
% end
% 
% % add static group with rate = 0
% y_s = rmnan(static_avg);
% datastats.all = [datastats.all; y_s];
% datastats.id  = [datastats.id;  zeros(length(y_s), 1)];
% 
% % one-way ANOVA
% [p, ~, stats] = anova1(datastats.all, datastats.id, 'off');
% fprintf('\n\nTemperature rate statistics — one-way ANOVA: p = %.4f\n\n', p)
% 
% % Bonferroni multiple comparisons
% alpha = 0.05;
% [c, ~, ~, gnames] = multcompare(stats, 'Alpha', alpha, 'Display', 'off');
% m_hyp        = size(c, 1);
% sigThreshold = alpha / m_hyp;
% 
% fprintf('Bonferroni corrected threshold: p < %.4f\n\n', sigThreshold)
% fprintf('%-20s %-20s %10s %10s\n', 'Group 1', 'Group 2', 'p-value', 'sig')
% for idx = 1:m_hyp
%     g1    = gnames{c(idx, 1)};
%     g2    = gnames{c(idx, 2)};
%     pval  = c(idx, 6);
%     isSig = pval <= sigThreshold;
%     fprintf('%-20s %-20s %10.4f %10s\n', g1, g2, pval, string(isSig))
% end

% -------------------------------------------------------------------------------------------------------

% STATS — each dynamic group vs static only
fprintf('\n\nPairwise comparisons vs static (0°C/min)\n\n')
fprintf('%-20s %10s %10s %10s\n', 'Group', 'p-value', 'p-corr', 'sig')

n_comparisons = num.exp;  % one test per dynamic group
alpha_corr  = alpha / n_comparisons;  % Bonferroni corrected threshold

y_static = rmnan(static_avg);

for ii = 1:num.exp
    y_dynamic = rmnan(plotData(:, ii));
    
    % two-sample t-test
    [~, pval] = ttest2(y_dynamic, y_static);
    
    isSig = pval <= alpha_corr;
    fprintf('%-20s %10.4f %10.4f %10s\n', ...
        num2str(temp_rate(ii)), pval, alpha_corr, string(isSig))
end
fprintf('\nBonferroni corrected threshold: p < %.4f\n', alpha_corr)

% add significance stars above dynamic groups
yl = ylim;
star_y = yl(2) + 0.05 * diff(yl);
for ii = 1:num.exp
    y_dynamic = rmnan(plotData(:, ii));
    [~, pval] = ttest2(y_dynamic, y_static);
    if pval <= alpha_corr
        text(temp_rate(ii), star_y, '*', ...
            'Color', foreColor, 'FontSize', 35, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Clipping', 'off')
    end
end

save_figure(fig, [paper_figs, 'temp rate vs ' pName ' scatter at coldest point']);

%% FIGURE: [temp rate tuning curves -- only the cooling data plotted together to show the temp-but-not-rate-dependence]


%% FIGURE: Heating & Cooling separated single parameter tuning curves -- select your metric
clearvars('-except',initial_vars{:})
[foreColor, ~] = formattingColors(blkbgd); % get background colors

autoLim = true;
y_lim = [0 100]; % if manually adjusting the axis limits for the y-axis 
xlim_auto = true; % change the time range for the x axis
nMax =  num.exp; 

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end,ext] = PlotParamSelection(false);
plot_err = true;

if isempty(title_str)
    return
end
fig_dir = [saveDir, 'temp tuning curves/'];
% set figure folder
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end

% % set up figure aligments
r = 1; %rows
c = 2; %columns

if contains(expGroup,'LTS 15-35')
    xlimits = [13, 37];
    xPos = [14.5, 25]; 
end

typeNames = {'cooling', 'warming'};
LW = 1.25;
% sSpan = 180;
dataString = cell([1,num.exp]);

% FIGURE:
fig = getfig('',true);
for i = 1:num.exp % num.exp:-1:1
    kolor = grouped(i).color;
    switch ext
        case true % subregions exist
            yy = grouped(i).(pName).food;
        case false % no subregions
            yy = grouped(i).(pName);
    end

    hold on
    if strcmp(pName, 'dist')
        x = grouped(i).(pName).distavgbytemp(:,1);
        YC = grouped(i).decreasing.all;
        YH = grouped(i).increasing.all;
    else
        x = yy.temps;
        YC = yy.decreasing.raw;
        YH = yy.increasing.raw;
    end
    % cooling
    subplot(r,c,1); hold on
    y = mean(YC,2,'omitnan')*scaler;
    y_err = (std(YC,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'color',kolor,'linewidth',LW)
    % heating
     subplot(r,c,2); hold on
    y = mean(YH,2,'omitnan')*scaler;
    y_err = (std(YH,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'color',kolor,'linewidth',LW)          

     dataString{i} = grouped(i).name;
end

% FORMATTING
for type = 1:2
    subplot(r, c, type)
    % formatting
    if type == 1
        set(gca, 'xdir', 'reverse')
        ylabel(y_lab)
    else
        set(gca, 'ycolor', 'none')
    end
    xlabel('temp (\circC)')
    %xlim(xlimits)
    % ylim(ylimits)
    title(typeNames{type},'color', foreColor)
end
matchAxis(fig, true); % match the y axes between heating and cooling
formatFig(fig, blkbgd,[r,c]);

for type = 1:2
    subplot(r, c, type)
    % ylim([0 90])
    h_line(nullD,'grey',':',2) %36.2
    set(gca,'ydir',y_dir)
    ylimits = y_lim;
    if contains(expGroup,'LTS 15-35')
        pos = [xPos(1,type), ylimits(1), 10, range(ylimits)]; % [lower-left X, lower-left Y, X-width, Y-height]
        h = rectangle('Position', pos, ...
                  'FaceColor', foreColor, ...   % RGB color
                  'FaceAlpha', 0.2, ...
                  'EdgeColor', 'none');
    end
    if type == 2
        set(gca, 'ycolor', 'none')
    end
end
legend(strrep(dataString,'_',' '), 'textcolor', foreColor, 'box', 'off','fontsize', 12,'location', 'northwest')

% add time arrows 
for type = 1:2 
    subplot(r,c,type)
    addTimeArrow(gca, foreColor)
end
set(findall(fig, 'Type', 'axes'), 'FontSize', 20)
set(findall(fig, 'Type', 'text'), 'FontSize', 20)
 
% save figure
save_figure(fig,[paper_figs  title_str ' temp tuning curves separeted H and C'],fig_type);

%% FIGURE: Heating & Cooling separated single parameter tuning curves -- select your metric
clearvars('-except',initial_vars{:})
[foreColor, ~] = formattingColors(blkbgd); % get background colors

autoLim = true;
y_lim = [0 100]; % if manually adjusting the axis limits for the y-axis 
xlim_auto = true; % change the time range for the x axis
nMax =  num.exp; 

% Select the type of information to plot: 
[title_str, pName,y_dir,y_lab,nullD,scaler,dType,dir_end,ext] = PlotParamSelection(false);
plot_err = true;

if isempty(title_str)
    return
end
fig_dir = [saveDir, 'temp tuning curves/'];
% set figure folder
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end

% % set up figure aligments
r = 1; %rows
c = 2; %columns

if contains(expGroup,'LTS 15-35')
    xlimits = [13, 37];
    xPos = [14.5, 25]; 
end

typeNames = {'cooling', 'warming'};
LW = 1.25;
% sSpan = 180;
dataString = cell([1,num.exp]);

% FIGURE:
fig = getfig('',0, [633 680]);
for i = 1:num.exp % num.exp:-1:1
    kolor = grouped(i).color;
    switch ext
        case true % subregions exist
            yy = grouped(i).(pName).food;
        case false % no subregions
            yy = grouped(i).(pName);
    end

    hold on
    if strcmp(pName, 'dist')
        x = grouped(i).(pName).distavgbytemp(:,1);
        YC = grouped(i).decreasing.all;
        YH = grouped(i).increasing.all;
    else
        x = yy.temps;
        YC = yy.decreasing.raw;
        YH = yy.increasing.raw;
    end
    % cooling
    subplot(r,c,1); hold on
    y = mean(YC,2,'omitnan')*scaler;
    y_err = (std(YC,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'color',kolor,'linewidth',LW)
    % heating
     subplot(r,c,2); hold on
    y = mean(YH,2,'omitnan')*scaler;
    y_err = (std(YH,0,2,'omitnan')*scaler)./sqrt(num.trial(i));
    plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.35);
    plot(x,y,'color',kolor,'linewidth',LW)          

     dataString{i} = grouped(i).name;
end

% FORMATTING
for type = 1:2
    subplot(r, c, type)
    % formatting
    if type == 1
        set(gca, 'xdir', 'reverse')
        ylabel(y_lab)
    else
        set(gca, 'ycolor', 'none')
    end
    xlabel('temp (\circC)')
    %xlim(xlimits)
    % ylim(ylimits)
    title(typeNames{type},'color', foreColor)
end
matchAxis(fig, true); % match the y axes between heating and cooling
formatFig(fig, blkbgd,[r,c]);

for type = 1:2
    subplot(r, c, type)
    % ylim([0 90])
    h_line(nullD,'grey',':',2) %36.2
    set(gca,'ydir',y_dir)
    ylimits = y_lim;
    if contains(expGroup,'LTS 15-35')
        pos = [xPos(1,type), ylimits(1), 10, range(ylimits)]; % [lower-left X, lower-left Y, X-width, Y-height]
        h = rectangle('Position', pos, ...
                  'FaceColor', foreColor, ...   % RGB color
                  'FaceAlpha', 0.2, ...
                  'EdgeColor', 'none');
    end
    if type == 2
        set(gca, 'ycolor', 'none')
    end

    if contains(expGroup, 'temp rate')
        ylim([-1.5, 40])
        set(gca, 'ytick', 0:10:40, 'xtick', 17:4:25)
    end
end

set(findall(fig, 'Type', 'axes'), 'FontSize', 20)
set(findall(fig, 'Type', 'text'), 'FontSize', 20)

for type = 1:2
    subplot(r, c, type)
    addTimeArrow(gca, foreColor, -0.08, -0.04, 0.04)
end

legend(strrep(dataString,'_',' '), 'textcolor', foreColor, 'box', 'off','fontsize', 12,'location', 'northwest')

% save figure
save_figure(fig,[paper_figs 'temp tuning curve h and c separated ' title_str figcolor(blkbgd)]);


%% [load from scratch -- sleep debt comparisons]

path_rate = "S:\Evyn\DATA\Grouped Data Structures\Berlin temp rate caviar\Sleep\Sleep duration by thermal stress norm axes-wht.fig";
path_shift = 






























