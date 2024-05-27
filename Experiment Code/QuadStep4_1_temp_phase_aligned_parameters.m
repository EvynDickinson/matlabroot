
%% Setup saving directory
fig_dir = [saveDir 'Phase Aligned/'];
if ~exist(fig_dir, "dir")
    mkdir(fig_dir)
end

initial_vars{end+1} = 'fig_dir';

%% Distance to food
clearvars('-except',initial_vars{:})

% PULL DATA: 
plotData = [];
for i = 1:num.exp
    % Cooling first
    temps_1 = flip(grouped(i).decreasing.temps);
    dist_1 = flip(grouped(i).decreasing.avg);
    err_1 = flip(grouped(i).decreasing.err);   
    loc = isnan(dist_1);
    temps_1(loc) = [];
    dist_1(loc) = [];
    err_1(loc) = [];
    % Warming
    temps_2 = grouped(i).increasing.temps;
    dist_2 = grouped(i).increasing.avg;
    err_2 = flip(grouped(i).decreasing.err);
    loc = isnan(dist_2);
    temps_2(loc) = [];
    dist_2(loc) = [];
    err_2(loc) = [];
    % combine for plotting
    if temps_1(end)==temps_2(1)
        % Average the lowest temp point... 
        low_dist = (dist_1(end)+dist_2(1))/2;
        low_err = (err_1(end)+err_2(1))/2;
        low_temp_idx = length(temps_1);
        plotData(i).X = [temps_1; temps_2(2:end)];
        plotData(i).Y = [dist_1(1:end-1); low_dist; dist_2(2:end)];
        plotData(i).err = [err_1(1:end-1); low_err; err_2(2:end)];
    end
    % save temp label
    x_tick_label = [];
    for tt = 1:length(plotData(i).X)
        x_tick_label{tt}  = num2str(plotData(i).X(tt));
    end
    plotData(i).xlabel = x_tick_label;
end

% PLOT DATA
LW = 2;

fig = getfig('',1);
hold on
for i = 1:num.exp
    x = 1:length(plotData(i).X);
    plot(x, plotData(i).Y,'color', grouped(i).color, 'linewidth', LW)
end
v_line(low_temp_idx,'grey', '--',1.5)
set(gca, 'XTick', 1:2:length(x_tick_label))
set(gca, 'XTickLabel', plotData(i).xlabel(1:2:end))
set(gca, 'YDir','reverse')
formatFig(fig, blkbgd);
xlabel('temperature (\circC)')
ylabel('proximity to food (mm)')

% Save figure
save_figure(fig,[fig_dir 'distance to food'],fig_type);

%% SPEED
clearvars('-except',initial_vars{:})

plotData = [];

for i = 1:num.exp
    x = data(i).G(1).TR.temps(1:end-1);
    [h_speed,c_speed] = deal([]);
    rateIdx = data(i).G(1).TR.rateIdx;

    % heating locations
    idx = find(data(i).G(1).TR.rates>0);
    heatloc = rateIdx==idx;
    h_style = '-';
    % cooling locations
    idx = find(data(i).G(1).TR.rates<0);
    coolloc = rateIdx==idx;
    %pull temp avg data
    tempIdx = data(i).G(1).TR.tempIdx;
    speed_data = grouped(i).speed.all;
    for t = 1:length(x)
        temp_loc = (tempIdx==t);
        % heating
        loc = temp_loc & heatloc;
        raw = mean(speed_data(loc,:),'omitnan');
        h_speed(t,:) =  raw;
        % cooling
        loc = temp_loc & coolloc;
        raw = mean(speed_data(loc,:),'omitnan');
        c_speed(t,:) =  raw;
    end
    
    % Average the lowest temp point... 
    temps_1 = flip(x);
    temps_2 = x;
    speed_1 = flip(mean(c_speed,2,'omitnan'));
    speed_2 = mean(h_speed,2,'omitnan');
    loc = isnan(speed_1);
    temps_1(loc) = [];
    speed_1(loc) = [];
    loc = isnan(speed_2);
    temps_2(loc) = [];
    speed_2(loc) = [];
    % avg low temp point:
    low_speed = (speed_1(end)+speed_2(1))/2;
    low_temp_idx = length(temps_1);
    plotData(i).X = [temps_1'; temps_2(2:end)'];
    plotData(i).Y = [speed_1(1:end-1); low_speed; speed_2(2:end)];

    % save temp label
    x_tick_label = [];
    for tt = 1:length(plotData(i).X)
        x_tick_label{tt}  = num2str(plotData(i).X(tt));
    end
    plotData(i).xlabel = x_tick_label;
end


% PLOT DATA
LW = 2;

fig = getfig('',1);
hold on
for i = 1:num.exp
    x = 1:length(plotData(i).X);
    plot(x, plotData(i).Y,'color', grouped(i).color, 'linewidth', LW)
end
v_line(low_temp_idx,'grey', '--',1.5)
set(gca, 'XTick', 1:2:length(x_tick_label))
set(gca, 'XTickLabel', plotData(i).xlabel(1:2:end))
formatFig(fig, blkbgd);
xlabel('temperature (\circC)')
ylabel('speed (mm/s)')

% Save figure
save_figure(fig,[fig_dir 'speed'],fig_type);











