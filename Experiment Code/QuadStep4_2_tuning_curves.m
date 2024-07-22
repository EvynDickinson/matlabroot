


%% FIGURE : distance tuning curve single point (must be done post-COM analysis)
clearvars('-except',initial_vars{:})

figDir = [saveDir 'Tuning Curves/'];
if ~isfolder(figDir)
    mkdir(figDir)
end

% Find master temperature list:
threshLow = [];
threshHigh = [];
for exp = 1:num.exp
    tp = getTempTurnPoints(data(exp).temp_protocol);
    threshLow = min([threshLow, tp.threshLow]);
    threshHigh = max([threshHigh, tp.threshHigh]);
end

% Temperature range: 
dynamic_temps = floor(threshLow):2:ceil(threshHigh);
nTemps = length(dynamic_temps);

% ============= Plot the data =============
buff_1 = 0.15; % buffer distance for scatter points
buff_2 = 0.35; % buffer distance for average line
offset_buffer = 0.1;
sz = 50;
LW = 2;
y_limits = [5,35];

fig = getfig('',1); hold on

for exp = 1:num.exp
    
    idx = 0+(offset_buffer*(exp-1));
    
    % Cooling first
    temp_list = flip(dynamic_temps);
    for tt = 1:nTemps
        idx = idx +1; % x-plotting position count so that the location is based on order
    
        % find frames that are at the appropriate temperature:
        cooling_idx = find(grouped(exp).position.temp_rates<0); %this assumes there is only one 'non-zero' temp rate
        temp_idx = find(grouped(exp).position.temp_list==temp_list(tt)); 
        if isempty(temp_idx) || isempty(cooling_idx)
            continue % skip locations without data for the temp or rate of change
        end
        %  frame locations
        frames = grouped(exp).position.loc(cooling_idx,temp_idx).frames;
        if isnan(frames)
            continue %skip if there are no frames for this condition
        end
        x_roi = [idx-buff_1, idx+buff_1];
        x_line = [idx-buff_2, idx+buff_2];
    
        % proximity to food for dynamic trials
        y_all = grouped(exp).dist.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        y_err = std(y_mean,0);
        scatter(idx, y_line, sz, grouped(exp).color, "filled")
        errorbar(idx, y_line, y_err, 'linewidth', LW, 'color', grouped(exp).color,'CapSize', 0)

        % % SCATTER VERSION:
        % x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        % scatter(x,y_mean, sz, grouped(exp).color,'filled')
        % plot(x_line,[y_line,y_line],'color', grouped(exp).color, 'linewidth', LW,'LineStyle','-')

    end 
    
    % Warming second
    temp_list = (dynamic_temps);
    for tt = 1:nTemps
        idx = idx +1;
    
        % find frames that are at the appropriate temperature:
        warming_idx = find(grouped(exp).position.temp_rates>0); %this assumes there is only one 'non-zero' temp rate
        temp_idx = find(grouped(exp).position.temp_list==temp_list(tt)); 
        if isempty(temp_idx) || isempty(warming_idx)
            continue % skip locations without data for the temp or rate of change
        end
        %  frame locations
        frames = grouped(exp).position.loc(warming_idx,temp_idx).frames;
        if all(isnan(frames))
            continue %skip if there are no frames for this condition
        end
        x_roi = [idx-buff_1, idx+buff_1];
        x_line = [idx-buff_2, idx+buff_2];
    
        % proximity to food for dynamic trials
        y_all = grouped(exp).dist.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        y_err = std(y_mean,0);
        scatter(idx, y_line, sz, grouped(exp).color, "filled")
        errorbar(idx, y_line, y_err, 'linewidth', LW, 'color', grouped(exp).color,'CapSize', 0)

        % % SCATTER VERSION:
        % x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        % scatter(x,y_mean, sz,  grouped(exp).color,'filled')
        % plot(x_line,[y_line,y_line],'color',  grouped(exp).color, 'linewidth', LW,'LineStyle','-')
    end
        
end


% Labels and formating
x_tick_label = [];
temp_names_label = [flip(dynamic_temps),dynamic_temps];
for i = 1:length(temp_names_label)
    x_tick_label{i}  = num2str(temp_names_label(i));
end
ylim(y_limits)
ax = gca;
set(ax, 'YDir', 'reverse')
set(ax, 'XTick', 1:nTemps*2,'XTickLabel', x_tick_label,'TickDir', 'out')
fig = formatFig(fig, blkbgd);
xlabel('Temperature (\circC)')
ylabel('proximity to food (mm)')

save_figure(fig,[figDir 'distance to food'],fig_type);

 

%% FIGURE : sleep tuning curve single point 
clearvars('-except',initial_vars{:})

figDir = [saveDir 'Tuning Curves/'];
if ~isfolder(figDir)
    mkdir(figDir)
end

% Find master temperature list:
threshLow = [];
threshHigh = [];
for exp = 1:num.exp
    tp = getTempTurnPoints(data(exp).temp_protocol);
    threshLow = min([threshLow, tp.threshLow]);
    threshHigh = max([threshHigh, tp.threshHigh]);
end

% Temperature range: 
dynamic_temps = floor(threshLow):2:ceil(threshHigh);
nTemps = length(dynamic_temps);

% ============= Plot the data =============
buff_1 = 0.15; % buffer distance for scatter points
buff_2 = 0.35; % buffer distance for average line
offset_buffer = 0.1;
sz = 50;
LW = 2;
% y_limits = [5,35]; %these should change for sleeping....

fig = getfig('',1); hold on

for exp = 1:num.exp
    
    idx = 0+(offset_buffer*(exp-1));
    
    % Cooling first
    temp_list = flip(dynamic_temps);
    for tt = 1:nTemps
        idx = idx +1; % x-plotting position count so that the location is based on order
    
        % find frames that are at the appropriate temperature:
        cooling_idx = find(grouped(exp).position.temp_rates<0); %this assumes there is only one 'non-zero' temp rate
        temp_idx = find(grouped(exp).position.temp_list==temp_list(tt)); 
        if isempty(temp_idx) || isempty(cooling_idx)
            continue % skip locations without data for the temp or rate of change
        end
        %  frame locations
        frames = grouped(exp).position.loc(cooling_idx,temp_idx).frames;
        if all(isnan(frames))
            continue %skip if there are no frames for this condition
        end
        x_roi = [idx-buff_1, idx+buff_1];
        x_line = [idx-buff_2, idx+buff_2];
    
        % proximity to food for dynamic trials
        % y_all = grouped(exp).dist.all(frames,:);
        y_all = sleep(exp).fract_sleep(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        y_err = std(y_mean,0);
        scatter(idx, y_line, sz, grouped(exp).color, "filled")
        errorbar(idx, y_line, y_err, 'linewidth', LW, 'color', grouped(exp).color,'CapSize', 0)

        % % SCATTER VERSION:
        % x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        % scatter(x,y_mean, sz, grouped(exp).color,'filled')
        % plot(x_line,[y_line,y_line],'color', grouped(exp).color, 'linewidth', LW,'LineStyle','-')
    end 
    
    % Warming second
    temp_list = (dynamic_temps);
    for tt = 1:nTemps
        idx = idx +1;
    
        % find frames that are at the appropriate temperature:
        warming_idx = find(grouped(exp).position.temp_rates>0); %this assumes there is only one 'non-zero' temp rate
        temp_idx = find(grouped(exp).position.temp_list==temp_list(tt)); 
        if isempty(temp_idx) || isempty(warming_idx)
            continue % skip locations without data for the temp or rate of change
        end
        %  frame locations
        frames = grouped(exp).position.loc(warming_idx,temp_idx).frames;
        if isnan(frames)
            continue %skip if there are no frames for this condition
        end
        x_roi = [idx-buff_1, idx+buff_1];
        x_line = [idx-buff_2, idx+buff_2];
    
        % proximity to food for dynamic trials
        y_all = sleep(exp).fract_sleep(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        y_err = std(y_mean,0);
        scatter(idx, y_line, sz, grouped(exp).color, "filled")
        errorbar(idx, y_line, y_err, 'linewidth', LW, 'color', grouped(exp).color,'CapSize', 0)

        % % SCATTER VERSION:
        % x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        % scatter(x,y_mean, sz,  grouped(exp).color,'filled')
        % plot(x_line,[y_line,y_line],'color',  grouped(exp).color, 'linewidth', LW,'LineStyle','-')
    end
        
end

% Labels and formating
x_tick_label = [];
temp_names_label = [flip(dynamic_temps),dynamic_temps];
for i = 1:length(temp_names_label)
    x_tick_label{i}  = num2str(temp_names_label(i));
end
% ylim(y_limits)
y_limits = ylim;
ylim([0,y_limits(2)]) % set the bottom of the y axis to zero since there can't be negative sleep

ax = gca;
% set(ax, 'YDir', 'reverse')
set(ax, 'XTick', 1:nTemps*2,'XTickLabel', x_tick_label,'TickDir','out')
fig = formatFig(fig, blkbgd);
xlabel('Temperature (\circC)')
ylabel('Fraction of flies sleeping')

save_figure(fig,[figDir 'fraction of flies sleeping'],fig_type);

%% FIGURE: (1/experiment group) Distance to food as a function of change in temperature from 25C
% right now this is only designed for LTS 15-35 trials
clearvars('-except',initial_vars{:})

figDir = [saveDir 'Tuning Curves/'];
if ~isfolder(figDir)
    mkdir(figDir)
end

PT = 25; % preferred temperature
sz = 15;

for exp = 1:num.exp
    
    % Plot data structure for easy looping
    PD = struct;
    kolors = {'seagreen', 'crimson'};
    conds = {'approach', 'retreat'};
    types = {'increasing', 'decreasing'}; %heating or cooling
    
    for cc = 1:2 % conditions = approach or retreat towards preferred temp
        PD.(conds{cc}) = []; %set empty fill space
        PD.(conds{cc}).color = Color(kolors{cc});
        for t = 1:2 % run both warming data and cooling data
            temps = grouped(exp).(types{t}).temps;
            % pull the locations for pooled data that are approaching or retreating 
            % from the preferred temperature:
            switch types{t}
                case 'increasing'
                    PD.approach.(types{t}).idx = temps<PT;
                    PD.retreat.(types{t}).idx = temps>PT;
                case 'decreasing'
                    PD.approach.(types{t}).idx = temps>PT;
                    PD.retreat.(types{t}).idx = temps<PT;
            end
            PD.(conds{cc}).(types{t}).temps = temps(PD.(conds{cc}).(types{t}).idx);
            PD.(conds{cc}).(types{t}).tempDiff = abs(PD.(conds{cc}).(types{t}).temps-PT);
            PD.(conds{cc}).(types{t}).dist = grouped(exp).(types{t}).all(PD.(conds{cc}).(types{t}).idx,:);
            PD.(conds{cc}).(types{t}).dist_avg = mean(PD.(conds{cc}).(types{t}).dist,2,'omitnan');
            PD.(conds{cc}).(types{t}).dist_err = std(PD.(conds{cc}).(types{t}).dist,0,2,'omitnan');% ./num.trial(exp);
        end
    end
    
    
    % FIGURE: 
    plot_err = true; 
    LW = 2;
    LT = {'-','--'};% warming = solid line | cooling = dashed line
    
    fig = getfig('', 1,[524 539]); hold on
    % plot avg + dev of each condition
    for cc = 1:2
        for t = 1:2 %1=warming, 2=cooling
    
            % pull params to plot
            kolor = PD.(conds{cc}).color;
            x = PD.(conds{cc}).(types{t}).tempDiff;
            y = PD.(conds{cc}).(types{t}).dist_avg;
            y_err = PD.(conds{cc}).(types{t}).dist_err;
            loc = isnan(y)|isnan(y_err);
            x(loc) = [];
            y(loc) = [];
            y_err(loc) = []; 
    
            % plot data
            plot(x,y,'color',kolor,'linewidth',LW,'LineStyle',LT{t})
            plot_error_fills(plot_err, x, y, y_err, kolor,  '-png', 0.35);
        end
    end
    
    % Formatting
    fig = formatFig(fig, blkbgd);
    title(grouped(exp).name)
    xlabel('|\DeltaT| from 25\circC')
    ylabel('Proximity to food (mm)')
    set(gca, 'YDir', 'reverse')
    xlim([0,10])
    ylim([5,35])
    h_line(18.1, 'grey',':',0.5)
    
    save_figure(fig,[figDir grouped(exp).name ' distance as absolute diff from preferred temp'],fig_type);

end


%% FIGURE: (1/experiment group) Sleep as a function of change in temperature from 25C
% right now this is only designed for LTS 15-35 trials
clearvars('-except',initial_vars{:})

figDir = [saveDir 'Tuning Curves/'];
if ~isfolder(figDir)
    mkdir(figDir)
end

PT = 25; % preferred temperature
sz = 15;

for exp = 1:num.exp
    
    % Plot data structure for easy looping
    PD = struct;
    kolors = {'seagreen', 'crimson'};
    conds = {'approach', 'retreat'};
    types = {'increasing', 'decreasing'}; %heating or cooling
    temps = sleep(exp).temps;

    for cc = 1:2 % conditions = approach or retreat towards preferred temp
        PD.(conds{cc}) = []; %set empty fill space
        PD.(conds{cc}).color = Color(kolors{cc});
        for t = 1:2 % run both warming data and cooling data
            % pull the locations for pooled data that are approaching or retreating 
            % from the preferred temperature:
            switch types{t}
                case 'increasing'
                    PD.approach.(types{t}).idx = temps<PT;
                    PD.retreat.(types{t}).idx = temps>PT;
                case 'decreasing'
                    PD.approach.(types{t}).idx = temps>PT;
                    PD.retreat.(types{t}).idx = temps<PT;
            end
            PD.(conds{cc}).(types{t}).temps = temps(PD.(conds{cc}).(types{t}).idx);
            PD.(conds{cc}).(types{t}).tempDiff = abs(PD.(conds{cc}).(types{t}).temps-PT);
            PD.(conds{cc}).(types{t}).dist_avg = sleep(exp).(types{t})(PD.(conds{cc}).(types{t}).idx,1);
            PD.(conds{cc}).(types{t}).dist_err = sleep(exp).(types{t})(PD.(conds{cc}).(types{t}).idx,2); % (STD)
        end
    end
    
    
    % FIGURE: 
    plot_err = true; 
    LW = 2;
    LT = {'-','--'};% warming = solid line | cooling = dashed line
    
    fig = getfig('', 1,[524 539]); hold on
    % plot avg + dev of each condition
    for cc = 1:2
        for t = 1:2 %1=warming, 2=cooling
    
            % pull params to plot
            kolor = PD.(conds{cc}).color;
            x = PD.(conds{cc}).(types{t}).tempDiff;
            y = PD.(conds{cc}).(types{t}).dist_avg;
            y_err = PD.(conds{cc}).(types{t}).dist_err;
            loc = isnan(y)|isnan(y_err);
            x(loc) = [];
            y(loc) = [];
            y_err(loc) = []; 
    
            % plot data
            plot(x,y,'color',kolor,'linewidth',LW,'LineStyle',LT{t})
            plot_error_fills(plot_err, x, y, y_err, kolor,  '-png', 0.35);
        end
    end
    
    % Formatting
    fig = formatFig(fig, blkbgd);
    title(grouped(exp).name)
    xlabel('|\DeltaT| from 25\circC')
    ylabel('Fraction of flies sleeping')
    xlim([0,10])
    y_lim = ylim;
    ylim([0,y_lim(2)])
    
    save_figure(fig,[figDir grouped(exp).name ' sleep as function of absolute diff from preferred temp'],fig_type);

end


%% FIGURE: speed vs eccentricity divided by approach vs retreat from PT
% right now this is only designed for LTS 15-35 trials
clearvars('-except',initial_vars{:})
[foreColor,backColor] = formattingColors(blkbgd);
figDir = [saveDir 'Tuning Curves/'];
if ~isfolder(figDir)
    mkdir(figDir)
end

PT = 25; % preferred temperature

for exp = 1:num.exp
    
    % Plot data structure for easy looping
    PD = struct;
    kolors = {'seagreen', 'crimson'};
    conds = {'approach', 'retreat'};
    types = {'increasing', 'decreasing'}; %heating or cooling
    
    for cc = 1:2 % conditions = approach or retreat towards preferred temp
        PD.(conds{cc}) = []; %set empty fill space
        PD.(conds{cc}).color = Color(kolors{cc});
        for t = 1:2 % run both warming data and cooling data
            temps = grouped(exp).position.temp_list;
            % pull the locations for pooled data that are approaching or retreating 
            % from the preferred temperature:
            switch types{t}
                case 'increasing'
                    PD.approach.(types{t}).idx = temps<PT;
                    PD.retreat.(types{t}).idx = temps>PT;
                    rate_idx = find(grouped(exp).position.temp_rates>0); % location in the position matrix for increasing frame locations
                case 'decreasing'
                    PD.approach.(types{t}).idx = temps>PT;
                    PD.retreat.(types{t}).idx = temps<PT;
                    rate_idx = find(grouped(exp).position.temp_rates<0); % location in the position matrix for decreasing frame locations
            end

            % pull all the frames that correspond to the condition (e.g., increasing & approach)
            temp_locs = find(PD.(conds{cc}).(types{t}).idx);
            speed_mat = [];
            eccent_mat = [];
            for bin = 1:length(temp_locs)
                frames = grouped(exp).position.loc(rate_idx, temp_locs(bin)).frames;
                if all(isnan(frames)) || isempty(frames)
                    % TODO: add nan placeholder first
                    continue
                end
                % avg speed per trial in this bin:
                x = mean(grouped(exp).speed.all(frames,:),1,'omitnan');
                speed_mat = [speed_mat; x'];
                
                % eccentricity for those frames (have to pull from data)
                for trial = 1:num.trial(exp)
                    y = data(exp).data(trial).occupancy.eccentricity(frames,1);
                    y = mean(y, 'omitnan');
                    eccent_mat = [eccent_mat; y];
                end
            end
            
            PD.(conds{cc}).(types{t}).speed = speed_mat;
            PD.(conds{cc}).(types{t}).eccentricity = eccent_mat;
        end
    end
    
    
    % % FIGURE: 
    % sz = 7;    
    % LT = {'-', ':'};
    % LW = 2;
    % 
    % fig = getfig('', 1,[524 539]); hold on
    %     % plot avg + dev of each condition
    %     for cc = 1:2 % approach, retreat
    %         for t = 1:2 
    %             kolor = PD.(conds{cc}).color;
    %             % all data points
    %             x = PD.(conds{cc}).(types{t}).speed;
    %             y = PD.(conds{cc}).(types{t}).eccentricity;
    %             x_err = std(x,0,'omitnan');%./num.trial(exp);
    %             x_avg = mean(x,'omitnan');
    %             y_err = std(y,0,'omitnan');%./num.trial(exp);
    %             y_avg = mean(y,'omitnan');
    %             switch types{t}
    %                 case 'increasing'
    %                     scatter(x,y,sz,kolor,'filled')
    %                 case 'decreasing'
    %                     scatter(x,y,sz,kolor)
    %             end
    %             % avg + error
    %             scatter(x_avg,y_avg,70,'MarkerFaceColor',kolor,'MarkerEdgeColor',foreColor,'LineWidth',1)
    %             errorbar(x_avg,y_avg,x_err, 'horizontal', 'LineStyle',LT{t},'color', kolor, 'linewidth', LW);
    %             errorbar(x_avg,y_avg,y_err, 'vertical', 'LineStyle', LT{t},'color', kolor, 'linewidth', LW);
    %         end
    %     end

% FIGURE:  Combined warming and cooling
    sz = 7;    
    LW = 3;

    fig = getfig('', 1,[524 539]); hold on
        % plot avg + dev of each condition
        for cc = 1:2 % approach, retreat
            kolor = PD.(conds{cc}).color;
            % pull the data points
            x = [PD.(conds{cc}).(types{1}).speed; PD.(conds{cc}).(types{2}).speed];
            y = [PD.(conds{cc}).(types{1}).eccentricity; PD.(conds{cc}).(types{2}).eccentricity];
            x_err = std(x,0,'omitnan');%./num.trial(exp);
            x_avg = mean(x,'omitnan');
            y_err = std(y,0,'omitnan');%./num.trial(exp);
            y_avg = mean(y,'omitnan');
            % plot all the points
            scatter(x,y,sz,kolor,'filled')
            % avg + error
            errorbar(x_avg,y_avg,x_err, 'horizontal','color', kolor, 'linewidth', LW);
            errorbar(x_avg,y_avg,y_err, 'vertical', 'color', kolor, 'linewidth', LW);
            scatter(x_avg,y_avg,70,'MarkerFaceColor',kolor,'MarkerEdgeColor',foreColor,'LineWidth',1)
        end
        
    
    % Formatting
    fig = formatFig(fig, blkbgd);
    title(grouped(exp).name)
    xlabel('speed (mm)')
    ylabel('distance from center of arena (mm)')
    % set(gca, 'YDir', 'reverse')
    xlim([0,15])
    ylim([14,28])
    % h_line(18.1, 'grey',':',0.5)
    
    save_figure(fig,[figDir grouped(exp).name ' speed vs eccentricity by approach and retreat'],fig_type);

end



























