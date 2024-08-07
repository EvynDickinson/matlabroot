
%% Setup saving directory
fig_dir = [saveDir 'Phase Aligned/'];
if ~exist(fig_dir, "dir")
    mkdir(fig_dir)
end

initial_vars{end+1} = 'fig_dir';

%% Distance to food
clearvars('-except',initial_vars{:})

switch questdlg('Absolute or normalized distance to food?','', 'Absolute', 'Normalized', 'Cancel', 'Absolute')
    case 'Absolute'
        ylab = 'proximity to food (mm)';
        save_str = 'absolute distance to food';
        ydir = 'reverse';
        y_abs = true;
    case 'Normalized'
        ylab = 'proximity to food (z-score)';
        ydir = 'reverse';
        save_str = 'normalized distance to food';
        y_abs = false;
    case 'Cancel'
        return
end

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
        if y_abs
            plotData(i).Y = [dist_1(1:end-1); low_dist; dist_2(2:end)];
        else
            plotData(i).Y = zscore([dist_1(1:end-1); low_dist; dist_2(2:end)]);
        end
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

fig = getfig('',1,[646 684]);
hold on
for i = 1:num.exp
    x = 1:length(plotData(i).X);
    plot(x, plotData(i).Y,'color', grouped(i).color, 'linewidth', LW)
end
v_line(low_temp_idx,'grey', '--',1.5)
idx = [1:4:low_temp_idx, low_temp_idx+1:4:length(x)];
set(gca, 'XTick', idx)
set(gca, 'XTickLabel', plotData(i).xlabel(idx))
set(gca, 'YDir',ydir)
formatFig(fig, blkbgd);
xlabel('temperature (\circC)')
ylabel(ylab)

% Save figure
save_figure(fig,[fig_dir save_str],fig_type);

%% Number of flies on the food
clearvars('-except',initial_vars{:})

switch questdlg('Absolute or normalized number of flies on food?','', 'Absolute', 'Normalized', 'Cancel', 'Absolute')
    case 'Absolute'
        ylab = 'number of flies on food';
        save_str = 'absolute flies on food';
        y_abs = true;
    case 'Normalized'
        ylab = 'flies on the food (z-score)';
        save_str = 'normalized flies on food';
        y_abs = false;
    case 'Cancel'
        return
end

% PULL DATA: 
plotData = [];
for i = 1:num.exp
    % Cooling first
    temps_1 = flip(grouped(i).FoF.temps);
    dist_1 = flip(grouped(i).FoF.decreasing(:,1));
    err_1 = flip(grouped(i).FoF.decreasing(:,2));   
    loc = isnan(dist_1);
    temps_1(loc) = [];
    dist_1(loc) = [];
    err_1(loc) = [];
    % Warming
    temps_2 = grouped(i).FoF.temps;
    dist_2 = grouped(i).FoF.increasing(:,1);
    err_2 = flip(grouped(i).FoF.decreasing(:,2));
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
        plotData(i).X = [temps_1'; temps_2(2:end)'];
        if y_abs
            plotData(i).Y = [dist_1(1:end-1); low_dist; dist_2(2:end)];
        else
            plotData(i).Y = zscore([dist_1(1:end-1); low_dist; dist_2(2:end)]);
        end
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

fig = getfig('',1,[646 684]);
hold on
for i = 1:num.exp
    x = 1:length(plotData(i).X);
    plot(x, plotData(i).Y,'color', grouped(i).color, 'linewidth', LW)
end
v_line(low_temp_idx,'grey', '--',1.5)
idx = [1:4:low_temp_idx, low_temp_idx+1:4:length(x)];
set(gca, 'XTick', idx)
set(gca, 'XTickLabel', plotData(i).xlabel(idx))
formatFig(fig, blkbgd);
xlabel('temperature (\circC)')
ylabel(ylab)

% Save figure
save_figure(fig,[fig_dir save_str],fig_type);

%% SPEED
clearvars('-except',initial_vars{:})

switch questdlg('Absolute or normalized speed?','', 'Absolute', 'Normalized', 'Cancel', 'Absolute')
    case 'Absolute'
        ylab = 'speed (mm/s)';
        save_str = 'absolute speed';
        y_abs = true;
    case 'Normalized'
        ylab = 'speed (z-score)';
        save_str = 'normalized speed';
        y_abs = false;
    case 'Cancel'
        return
end

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
    if y_abs
        plotData(i).Y = [speed_1(1:end-1); low_speed; speed_2(2:end)];
    else
        plotData(i).Y = zscore([speed_1(1:end-1); low_speed; speed_2(2:end)]);
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

fig = getfig('',1,[646 684]);
hold on
for i = 1:num.exp
    x = 1:length(plotData(i).X);
    plot(x, plotData(i).Y,'color', grouped(i).color, 'linewidth', LW)
end
v_line(low_temp_idx,'grey', '--',1.5)
idx = [1:4:low_temp_idx, low_temp_idx+1:4:length(x)];
set(gca, 'XTick', idx)
set(gca, 'XTickLabel', plotData(i).xlabel(idx))
formatFig(fig, blkbgd);
xlabel('temperature (\circC)')
ylabel(ylab)

% Save figure
save_figure(fig,[fig_dir save_str],fig_type);

%% SLEEP
clearvars('-except',initial_vars{:})

switch questdlg('Absolute or normalized sleep?','', 'Absolute', 'Normalized', 'Cancel', 'Absolute')
    case 'Absolute'
        ylab = 'sleep (fraction of flies)';
        save_str = 'absolute sleep';
        y_abs = true;
    case 'Normalized'
        ylab = 'sleep (z-score)';
        save_str = 'normalized sleep';
        y_abs = false;
    case 'Cancel'
        return
end

plotData = [];

for i = 1:num.exp
    % Average the lowest temp point... 
    temps_1 = flip(sleep(i).temps);
    temps_2 = sleep(i).temps;
    sleep_1 = flip(sleep(i).decreasing(:,1));
    sleep_2 = sleep(i).increasing(:,1);
    loc = isnan(sleep_1);
    temps_1(loc) = [];
    sleep_1(loc) = [];
    loc = isnan(sleep_2);
    temps_2(loc) = [];
    sleep_2(loc) = [];
    % avg low temp point:
    low_sleep = (sleep_1(end)+sleep_2(1))/2;
    low_temp_idx = length(temps_1);
    plotData(i).X = [temps_1'; temps_2(2:end)'];
    if y_abs
        plotData(i).Y = [sleep_1(1:end-1); low_sleep; sleep_2(2:end)];
    else
        plotData(i).Y = zscore([sleep_1(1:end-1); low_sleep; sleep_2(2:end)]);
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

fig = getfig('',1,[646 684]);
hold on
for i = 1:num.exp
    x = 1:length(plotData(i).X);
    plot(x, plotData(i).Y,'color', grouped(i).color, 'linewidth', LW)
end
v_line(low_temp_idx,'grey', '--',1.5)
idx = [1:4:low_temp_idx, low_temp_idx+1:4:length(x)];
set(gca, 'XTick', idx)
set(gca, 'XTickLabel', plotData(i).xlabel(idx))
formatFig(fig, blkbgd);
xlabel('temperature (\circC)')
ylabel(ylab)

% Save figure
save_figure(fig,[fig_dir save_str],fig_type);





%% Overlay by z-score
clearvars('-except',initial_vars{:})
LW = 2;

for i = 1:num.exp
    fig = getfig('',1,[482 684]);
        hold on
        % =========================================
        % SLEEP (white)
        % Average the lowest temp point... 
        temps_1 = flip(sleep(i).temps);
        temps_2 = sleep(i).temps;
        sleep_1 = flip(sleep(i).decreasing(:,1));
        sleep_2 = sleep(i).increasing(:,1);
        loc = isnan(sleep_1);
        temps_1(loc) = [];
        sleep_1(loc) = [];
        loc = isnan(sleep_2);
        temps_2(loc) = [];
        sleep_2(loc) = [];
        % avg low temp point:
        low_sleep = (sleep_1(end)+sleep_2(1))/2;
        low_temp_idx = length(temps_1);
        X = [temps_1'; temps_2(2:end)'];
        Y = [sleep_1(1:end-1); low_sleep; sleep_2(2:end)];
        % normalize to 1 and plot
        Y = zscore(Y);
        x = 1:length(X);
        plot(x , Y,'color', 'white', 'linewidth', LW)
        % =========================================
        
        % =========================================
        % SPEED (teal)
        x = data(i).G(1).TR.temps(1:end-1);
        [h_speed,c_speed] = deal([]);
        rateIdx = data(i).G(1).TR.rateIdx;
    
        % heating locations
        idx = find(data(i).G(1).TR.rates>0);
        heatloc = rateIdx==idx;
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
        X = [temps_1'; temps_2(2:end)'];
        Y = [speed_1(1:end-1); low_speed; speed_2(2:end)];
        
        % normalize and plot:
        Y = zscore(Y);
        x = 1:length(X);
        plot(x, Y, 'color', Color('teal'), 'linewidth', LW)
        % =========================================
        
        % =========================================
        % DISTANCE (orange)
        % Cooling first
        temps_1 = flip(grouped(i).decreasing.temps);
        dist_1 = flip(grouped(i).decreasing.avg);
        loc = isnan(dist_1);
        temps_1(loc) = [];
        dist_1(loc) = [];
        % Warming
        temps_2 = grouped(i).increasing.temps;
        dist_2 = grouped(i).increasing.avg;
        loc = isnan(dist_2);
        temps_2(loc) = [];
        dist_2(loc) = [];
        % combine for plotting
        if temps_1(end)==temps_2(1)
            % Average the lowest temp point... 
            low_dist = (dist_1(end)+dist_2(1))/2;
            low_temp_idx = length(temps_1);
            x = [temps_1; temps_2(2:end)];
            y = [dist_1(1:end-1); low_dist; dist_2(2:end)];
        end
        % normalize and plot
        Y = zscore(y);
        X = 1:length(x);
        plot(X, Y, 'color', Color('orangered'), 'linewidth', LW)
        % =========================================
      
        % FLIES ON FOOD (purple) 
        % Cooling first
        temps_1 = flip(grouped(i).FoF.temps);
        dist_1 = flip(grouped(i).FoF.decreasing(:,1));
        loc = isnan(dist_1);
        temps_1(loc) = [];
        dist_1(loc) = [];
        % Warming
        temps_2 = grouped(i).FoF.temps;
        dist_2 = grouped(i).FoF.increasing(:,1);
        loc = isnan(dist_2);
        temps_2(loc) = [];
        dist_2(loc) = [];
        % combine for plotting
        if temps_1(end)==temps_2(1)
            % Average the lowest temp point... 
            low_dist = (dist_1(end)+dist_2(1))/2;
            low_temp_idx = length(temps_1);
            X = [temps_1'; temps_2(2:end)'];
            Y = zscore([dist_1(1:end-1); low_dist; dist_2(2:end)]);
        end
        x = 1:length(X);
        plot(x, Y, 'color', Color('blueviolet'), 'linewidth', LW)

        % FORMATING
        % save temp label
        x_tick_label = [];
        for tt = 1:length(x)
            x_tick_label{tt}  = num2str(X(tt));
        end
        
        v_line(low_temp_idx,'grey', '--',1.5)
        idx = [1:4:low_temp_idx, low_temp_idx+1:4:length(x)];
        set(gca, 'XTick', idx)
        set(gca, 'XTickLabel', x_tick_label(idx))
        formatFig(fig, blkbgd);
        xlabel('temperature (\circC)')
        ylabel('z-score')
        title(grouped(i).name,'color', 'w')
        
    % Save figure
    save_figure(fig,[fig_dir grouped(i).name ' zscore'],fig_type);
    
end


































