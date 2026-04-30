
blkbgd = false;
[foreColor,backColor] = formattingColors(blkbgd); %get background colors


%% FIGURE: time course of all the different behaviors ordered by precidence
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

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
saveDir = createFolder([figDir, 'Paper Figures/']);
save_figure(fig, [saveDir, 'Fig 1 timecourse of behavior ', grouped(g).name]);

%% FIGURE: timecourse of 







%%
autoLim = true;
% Y limit ranges
speed_lim = [0,10]; %speed
dist_lim = [5, 30]; %distance
dt_lim = [12, 34];      %distance-temp
% dt_lim = [10, 30];        %distance-temp
nMax =  num.exp;%
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

% set up figure aligments
r = 5; %rows
c = 3; %columns
sb(1).idx = 1:2; %temp timecourse
sb(2).idx = [4,5,7,8]; %distance from food timecourse %TODO: normalize this to something more intuitive?
sb(3).idx = [10,11,13,14]; %speed timecourse
sb(4).idx = 3:3:15; %binned distance alignment

LW = 0.75;
sSpan = 360; % 2 minute smoothing 
dataString = cell([1,num.exp]);

% FIGURE:  
fig = getfig('',true);
for g = 1:nMax
%     i = expOrder(ii);
    x = grouped(g).time;   
    kolor = grouped(g).color;

    %temp
    subplot(r,c,sb(1).idx); hold on
        y = grouped(g).temp;
        plot(x,y,'LineWidth',1,'Color',kolor)

    %distance
    subplot(r,c,sb(2).idx); hold on
        y = smooth(grouped(g).dist.avg,'moving',sSpan);
        % for jj = 1:50
        %     y = smooth(y,'moving',sSpan);
        % end
%         y_err = smooth(grouped(i).dist.err,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor)
        if ~autoLim
            ylim(dist_lim)
        end

%     %speed
    subplot(r,c,sb(3).idx); hold on
        y = smooth(grouped(g).speed.avg,'moving',sSpan);
        % for jj = 1:50
        %     y = smooth(y,'moving',sSpan);
        % end
        %         y_err = smooth(grouped(i).speed.err,'moving',sSpan);
        plot(x,y,'LineWidth',LW,'Color',kolor)
        if ~autoLim
            ylim(speed_lim)
        end

    %temp dependent distance
    subplot(r,c,sb(4).idx); hold on
        x = grouped(g).dist.distavgbytemp(:,1);
        y = grouped(g).dist.distavgbytemp(:,2);
        y_err = grouped(g).dist.distavgbytemp_err(:,2);
        plot_error_fills(plot_err, x, y, y_err, kolor,  fig_type, 0.4);
        plot(x,y,'color',kolor,'linewidth',LW+1)

        dataString{g} = grouped(g).name;
end

% FORMATING AND LABELS
formatFig(fig,blkbgd,[r,c],sb);
% temp
subplot(r,c,sb(1).idx)
ylabel('\circC')
set(gca,"XColor",backColor)
% distance
subplot(r,c,sb(2).idx)
ylabel('proximity to food (mm)')
set(gca,"XColor",backColor)
set(gca,'ydir','reverse')
% speed
subplot(r,c,sb(3).idx)
ylabel('speed (mm/s)')
xlabel('time (min)')
% temp-distance relationship
subplot(r,c,sb(4).idx)
ylabel('proximity to food (mm)')

xlabel('temp (\circC)')
if ~autoLim
    ylim(dt_lim)
end
h_line(18.1,'grey',':',1) %36.2
set(gca,'ydir','reverse')
%
dataString = strrep(dataString, '_', ' ');
legend(dataString,'textcolor', foreColor, 'location', 'northwest', 'box', 'off','fontsize', 5)

% save figure
save_figure(fig,[saveDir expGroup ' timecourse summary'],fig_type);