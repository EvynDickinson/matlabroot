 
% need to run prep code 'normalization distance' in COM 
clearvars('-except',initial_vars{:})


%% FIGURE: proximity to food for static and dynamic W & C separated 
buff_1 = 0.25;
buff_2 = 0.5;
sz = 60;
LW = 2;
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

temps_idx = [3, 9, 15, 19];  % 17,20,23,25 [temp locations that match the assigned hold temps]
temp_actual = grouped(4).position.temp_list(temps_idx);

dyn_exp = [4]; %start with just the caviar trial
stat_idx = [1,2,3,6];

fig = getfig('',1); hold on

for tt = 1:length(temps_idx)
    %  frame locations
    frames_up = grouped(4).position.loc(3,temps_idx(tt)).frames;
    frames_down = grouped(4).position.loc(1,temps_idx(tt)).frames;
    x_roi = [temp_actual(tt)-buff_1, temp_actual(tt)+buff_1];
    x_line = [temp_actual(tt)-buff_2, temp_actual(tt)+buff_2];

    %proximity to food for dynamic trials
    for i = dyn_exp
        % warming
        y_all = grouped(i).dist.all(frames_up,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, grouped(i).color,'filled')
        plot(x_line,[y_line,y_line],'color', grouped(i).color, 'linewidth', LW,'LineStyle','-')
        % cooling
        y_all = grouped(i).dist.all(frames_down,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        scatter(x,y_mean, sz, grouped(i).color)
        plot(x_line,[y_line,y_line],'color', grouped(i).color, 'linewidth', LW,'LineStyle','--')
    end

    % plot concordant static trial
    i = stat_idx(tt);
    y_all = grouped(i).dist.all([frames_up; frames_down],:);
    y_mean = mean(y_all, 1,'omitnan');
    y_line = mean(y_mean);
    x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
    scatter(x,y_mean, sz, foreColor,'filled')
    plot(x_line,[y_line,y_line],'color', foreColor, 'linewidth', LW,'LineStyle','-')
end

% Figure formatting
formatFig(fig, blkbgd);
h_line(18.1, 'grey',':',0.5) 
xlim([16,26])
ylabel('proximity to food (mm)')
xlabel('temperature (\circC)')
ax = gca;
set(gca,'ydir','reverse','XTick',temp_actual)
ax.FontSize = 18;

save_figure(fig,[saveDir 'static vs dynamic distance to food warm cool separated'],fig_type);

%% FIGURE: distance from food collapsed across heating/cooling
buff_1 = 0.25;
buff_2 = 0.5;
sz = 60;
LW = 2;
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

temps_idx = [3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
temp_actual = grouped(4).position.temp_list(temps_idx);

dyn_exp = [4]; %start with just the caviar trial
stat_idx = [1,2,3,6];

fig = getfig('',1); hold on

for tt = 1:length(temps_idx)
    %  frame locations
    frames_up = grouped(4).position.loc(3,temps_idx(tt)).frames;
    frames_down = grouped(4).position.loc(1,temps_idx(tt)).frames;
    x_roi = [temp_actual(tt)-buff_1, temp_actual(tt)+buff_1];
    x_line = [temp_actual(tt)-buff_2, temp_actual(tt)+buff_2];

    % %proximity to food for dynamic trials
    % for i = dyn_exp
    %     % warming
    %     y_all = grouped(i).dist.all([frames_up; frames_down],:);
    %     y_mean = mean(y_all, 1,'omitnan');
    %     y_line = mean(y_mean);
    %     x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
    %     scatter(x,y_mean, sz, grouped(i).color,'filled')
    %     plot(x_line,[y_line,y_line],'color', grouped(i).color, 'linewidth', LW,'LineStyle','-')
    % end
    % 
    % plot concordant static trial
    i = stat_idx(tt);
    y_all = grouped(i).dist.all([frames_up; frames_down],:);
    y_mean = mean(y_all, 1,'omitnan');
    y_line = mean(y_mean);
    x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
    scatter(x,y_mean, sz, foreColor,'filled')
    plot(x_line,[y_line,y_line],'color', foreColor, 'linewidth', LW,'LineStyle','-')
end

% Figure formatting
formatFig(fig, blkbgd);
h_line(18.1, 'grey',':',0.5) 
xlim([16,26])
ylim([5,35])
ylabel('proximity to food (mm)')
xlabel('temperature (\circC)')
ax = gca;
set(gca,'ydir','reverse','XTick',temp_actual)
ax.FontSize = 18;

clearvars('-except',initial_vars{:})

%% FIGURE: animation style overview 
%   close all
fig_dir = [saveDir '/PPT FIgures/Distance/'];
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end
plot_caviar = 0;
plot_static = 1;
plot_water = 0;
plot_waxed = 0;
others = 0;

% INDEX LIST
caviar_idx = 4;
waxed_idx = 5;
water_idx = 7;
stat_idx = [1,2,3,6]; % 17, 20, 23, 25
empty_idx = 8;
AN_MP_idx = 9;
arista_idx = 10;

% plot multioptions
plot_list = 1:num.exp;
plot_list(stat_idx) = [];

% COLOR LIST
static_color = Color('gold');
grouped(water_idx).color = Color('white');
grouped(waxed_idx).color = Color('lime');
grouped(caviar_idx).color = Color('dodgerblue');
grouped(empty_idx).color = Color('navajowhite');
grouped(arista_idx).color = Color('DarkOrchid');
grouped(AN_MP_idx).color = Color('red');

% ploting parameters
buff_1 = 0.25;
buff_2 = 0.5;
sz = 40;
LW = 3;
x_limits = [0,19];
y_limits = [5,35];
% [foreColor,backColor] = formattingColors(blkbgd); %get background colors

% experiment indexes:

static_temps = [3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
% static_temps = [19,15,9, 3, 3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
static_temp_actual = grouped(caviar_idx).position.temp_list(static_temps);

dynamic_temps = 3:2:19; % all whole number temps
temp_actual = grouped(caviar_idx).position.temp_list(dynamic_temps);
temp_names_label = [flip(temp_actual),temp_actual];
for i = 1:length(temp_names_label)
    x_tick_label{i}  = num2str(temp_names_label(i));
end

fig = getfig('',1); hold on
    %empty axis:
    formatFig(fig, blkbgd);
    xlim(x_limits)
    ylim(y_limits)
    h_line(18.1, 'grey',':',0.5) 
    
    ylabel('proximity to food (mm)')
    xlabel('temperature (\circC)')
    ax = gca;
    set(ax,'ydir','reverse')
    set(ax,'XTick', 1:length(x_tick_label), "XTickLabel", x_tick_label)
    ax.FontSize = 18;

% save_figure(gcf,[fig_dir 'static vs dynamic distance to food empty frame'],fig_type, false, false);

% ==============================================================================
% 1) FULL EXP WARMING AND COOLING
caviar_color = grouped(caviar_idx).color;
if plot_caviar
    idx = 0;
    % Cooling first
    temp_list = flip(dynamic_temps);
    for tt = 1:length(dynamic_temps)
        idx = idx +1;
    
        %  frame locations
        frames = grouped(caviar_idx).position.loc(1,temp_list(tt)).frames;
        x_roi = [idx-buff_1, idx+buff_1];
        x_line = [idx-buff_2, idx+buff_2];
    
        % proximity to food for dynamic trials
        y_all = grouped(caviar_idx).dist.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, caviar_color,'filled')
        plot(x_line,[y_line,y_line],'color', caviar_color, 'linewidth', LW,'LineStyle','-')
    end
    
    % Warming second
    temp_list = (dynamic_temps);
    for tt = 1:length(dynamic_temps)
        idx = idx +1;
    
        %  frame locations
        frames = grouped(4).position.loc(3,temp_list(tt)).frames;
        x_roi = [idx-buff_1, idx+buff_1];
        x_line = [idx-buff_2, idx+buff_2];
    
        % proximity to food for dynamic trials
        y_all = grouped(caviar_idx).dist.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, caviar_color,'filled')
        plot(x_line,[y_line,y_line],'color', caviar_color, 'linewidth', LW,'LineStyle','-')
    end
        
    save_figure(gcf,[fig_dir 'static vs dynamic distance to food caviar'],fig_type, false, false);
end

% ==============================================================================
% 2) STATIC HOLD TEMPERATURES
if plot_static
    idx_loc = [1,3,6,9,10,13,16,18]; % x positions of the selected temperatures
    idx = 0;
    % Cooling first
    temp_list = flip(static_temps);
    expList = flip(stat_idx);
    for tt = 1:length(static_temps)
        idx = idx +1; 
        %  frame locations
        frames = grouped(4).position.loc(1,temp_list(tt)).frames; % 1=cooling, 3=warming
        x_roi = [idx_loc(idx)-buff_1, idx_loc(idx)+buff_1];
        x_line = [idx_loc(idx)-buff_2, idx_loc(idx)+buff_2];
        % pull data
        i = expList(tt);
        y_all = grouped(i).dist.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, static_color,'filled')
        plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
    end
    
    % warming second
    temp_list = (static_temps);
    expList = (stat_idx);
    for tt = 1:length(static_temps)
        idx = idx +1; 
        %  frame locations
        frames = grouped(4).position.loc(3,temp_list(tt)).frames; % 1=cooling, 3=warming
        x_roi = [idx_loc(idx)-buff_1, idx_loc(idx)+buff_1];
        x_line = [idx_loc(idx)-buff_2, idx_loc(idx)+buff_2];
        % pull data
        i = expList(tt);
        y_all = grouped(i).dist.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, static_color,'filled')
        plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
    end
    
    save_figure(gcf,[fig_dir 'static vs dynamic distance to food static'],fig_type, false, true);
end

% ==============================================================================
% WATER COMPARISON
if plot_water
    water_color = grouped(water_idx).color;
    % 1) FULL EXP WARMING AND COOLING
    idx = 0;
    % Cooling first
    temp_list = flip(dynamic_temps);
    for tt = 1:length(dynamic_temps)
        idx = idx +1;
    
        %  frame locations
        frames = grouped(water_idx).position.loc(1,temp_list(tt)).frames;
        x_roi = [idx-buff_1, idx+buff_1];
        x_line = [idx-buff_2, idx+buff_2];
    
        % proximity to food for dynamic trials
        y_all = grouped(water_idx).dist.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, water_color,'filled')
        plot(x_line,[y_line,y_line],'color', water_color, 'linewidth', LW,'LineStyle','-')
    end
    
    % Warming second
    temp_list = (dynamic_temps);
    for tt = 1:length(dynamic_temps)
        idx = idx +1;
        %  frame locations
        frames = grouped(water_idx).position.loc(3,temp_list(tt)).frames;
        x_roi = [idx-buff_1, idx+buff_1];
        x_line = [idx-buff_2, idx+buff_2];
    
        % proximity to food for dynamic trials
        y_all = grouped(water_idx).dist.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, water_color,'filled')
        plot(x_line,[y_line,y_line],'color', water_color, 'linewidth', LW,'LineStyle','-')
    end
    
    save_figure(gcf,[fig_dir 'static vs dynamic distance to food water'],fig_type, false, false);

end

% ==============================================================================
% WAXED COMPARISON
if plot_waxed
    waxed_color = grouped(waxed_idx).color;
    % 1) FULL EXP WARMING AND COOLING
    idx = 0;
    % Cooling first
    temp_list = flip(dynamic_temps);
    for tt = 1:length(dynamic_temps)
        idx = idx +1;
    
        %  frame locations
        frames = grouped(waxed_idx).position.loc(1,temp_list(tt)).frames;
        x_roi = [idx-buff_1, idx+buff_1];
        x_line = [idx-buff_2, idx+buff_2];
    
        % proximity to food for dynamic trials
        y_all = grouped(waxed_idx).dist.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, waxed_color,'filled')
        plot(x_line,[y_line,y_line],'color', waxed_color, 'linewidth', LW,'LineStyle','-')
    end
    
    % Warming second
    temp_list = (dynamic_temps);
    for tt = 1:length(dynamic_temps)
        idx = idx +1;
        %  frame locations
        frames = grouped(waxed_idx).position.loc(3,temp_list(tt)).frames;
        x_roi = [idx-buff_1, idx+buff_1];
        x_line = [idx-buff_2, idx+buff_2];
    
        % proximity to food for dynamic trials
        y_all = grouped(waxed_idx).dist.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, waxed_color,'filled')
        plot(x_line,[y_line,y_line],'color', waxed_color, 'linewidth', LW,'LineStyle','-')
    end
    
    save_figure(gcf,[fig_dir 'static vs dynamic distance to food waxed'],fig_type, false, false);

end

% ==============================================================================


% ALL OTHERS COMPARISON
if others
   for exp = 1:length(plot_list) 
    fig = getfig('',1); hold on
        %empty axis:
        formatFig(fig, blkbgd);
        xlim(x_limits)
        ylim(y_limits)
        h_line(18.1, 'grey',':',0.5) 
        
        ylabel('proximity to food (mm)')
        xlabel('temperature (\circC)')
        ax = gca;
        set(ax,'ydir','reverse')
        set(ax,'XTick', 1:length(x_tick_label), "XTickLabel", x_tick_label)
        ax.FontSize = 18;
    
        ii = plot_list(exp);
        % 1) FULL EXP WARMING AND COOLING
        idx = 0;
        % Cooling first
        temp_list = flip(dynamic_temps);
        for tt = 1:length(dynamic_temps)
            idx = idx +1;
        
            %  frame locations
            frames = grouped(ii).position.loc(1,temp_list(tt)).frames;
            x_roi = [idx-buff_1, idx+buff_1];
            x_line = [idx-buff_2, idx+buff_2];
        
            % proximity to food for dynamic trials
            y_all = grouped(ii).dist.all(frames,:);
            y_mean = mean(y_all, 1,'omitnan');
            y_line = mean(y_mean);
            x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
            scatter(x,y_mean, sz, grouped(ii).color,'filled')
            plot(x_line,[y_line,y_line],'color', grouped(ii).color, 'linewidth', LW,'LineStyle','-')
        end
        
        % Warming second
        temp_list = (dynamic_temps);
        for tt = 1:length(dynamic_temps)
            idx = idx +1;
            %  frame locations
            frames = grouped(ii).position.loc(3,temp_list(tt)).frames;
            x_roi = [idx-buff_1, idx+buff_1];
            x_line = [idx-buff_2, idx+buff_2];
        
            % proximity to food for dynamic trials
            y_all = grouped(ii).dist.all(frames,:);
            y_mean = mean(y_all, 1,'omitnan');
            y_line = mean(y_mean);
            x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
            scatter(x,y_mean, sz, grouped(ii).color,'filled')
            plot(x_line,[y_line,y_line],'color', grouped(ii).color, 'linewidth', LW,'LineStyle','-')
        end

        save_figure(fig,[fig_dir 'static vs dynamic distance to food ' grouped(ii).name],fig_type);
    
    end
end


%% ANALYSIS: flies on food:
well_radius = 3; % 5 mm diameter of the physical well -- give 0.5mm buffer zone outside well
well_rad = well_radius * pix2mm; %convert mm to pixels

% Calculate number of flies on food for each frame
% plotData = [];
for i = 1:num.exp
    N = [];
    for trial = 1:num.trial(i)
        onFood = [];
        well_loc = data(i).T.foodLoc(trial);
        well_center = data(i).data(trial).data.wellcenters(:,well_loc);
        x = data(i).data(trial).data.x_loc;
        y = data(i).data(trial).data.y_loc;
        for frame = 1:size(x,1)
            distances = sqrt((x(frame,:)-well_center(1)).^2+((y(frame,:)-well_center(2)).^2));
            onFood(frame,1) = sum(distances<well_rad);
        end
        N = autoCat(N,onFood,false);
    end
    grouped(i).fliesonfood = N;
    disp(['Done exp ' num2str(i)])
end

%% FIGURE: Flies on Food
%   close all
clearvars('-except',initial_vars{:})
fig_dir = [saveDir '/PPT FIgures/Flies on Food/'];
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end

plot_static = 1;
others = 1;

% INDEX LIST
caviar_idx = 4;
waxed_idx = 5;
water_idx = 7;
stat_idx = [1,2,3,6]; % 17, 20, 23, 25
empty_idx = 8;
AN_MP_idx = 9;
arista_idx = 10;

% plot multioptions
plot_list = 1:num.exp;
plot_list(stat_idx) = [];

% COLOR LIST
static_color = Color('gold');

% ploting parameters
buff_1 = 0.25;
buff_2 = 0.5;
sz = 40;
LW = 3;
x_limits = [0,19];
y_limits = [0, 4];
% [foreColor,backColor] = formattingColors(blkbgd); %get background colors

% experiment indexes:
static_temps = [3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
% static_temps = [19,15,9, 3, 3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
static_temp_actual = grouped(caviar_idx).position.temp_list(static_temps);

dynamic_temps = 3:2:19; % all whole number temps
temp_actual = grouped(caviar_idx).position.temp_list(dynamic_temps);
temp_names_label = [flip(temp_actual),temp_actual];
for i = 1:length(temp_names_label)
    x_tick_label{i}  = num2str(temp_names_label(i));
end

fig = getfig('',1); hold on
    %empty axis:
    formatFig(fig, blkbgd);
    xlim(x_limits)
    ylim(y_limits)

    ylabel('flies on food (#)')
    xlabel('temperature (\circC)')
    ax = gca;
    set(ax,'XTick', 1:length(x_tick_label), "XTickLabel", x_tick_label)
    ax.FontSize = 18;

% save_figure(gcf,[fig_dir 'empty frame'],fig_type, false, false);

% ==============================================================================
% STATIC HOLD TEMPERATURES
if plot_static
    idx_loc = [1,3,6,9,10,13,16,18]; % x positions of the selected temperatures
    idx = 0;
    % Cooling first
    temp_list = flip(static_temps);
    expList = flip(stat_idx);
    for tt = 1:length(static_temps)
        idx = idx +1; 
        %  frame locations
        frames = grouped(4).position.loc(1,temp_list(tt)).frames; % 1=cooling, 3=warming
        x_roi = [idx_loc(idx)-buff_1, idx_loc(idx)+buff_1];
        x_line = [idx_loc(idx)-buff_2, idx_loc(idx)+buff_2];
        % pull data
        i = expList(tt);
        y_all = grouped(i).fliesonfood(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, static_color,'filled')
        plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
    end
    
    % warming second
    temp_list = (static_temps);
    expList = (stat_idx);
    for tt = 1:length(static_temps)
        idx = idx +1; 
        %  frame locations
        frames = grouped(4).position.loc(3,temp_list(tt)).frames; % 1=cooling, 3=warming
        x_roi = [idx_loc(idx)-buff_1, idx_loc(idx)+buff_1];
        x_line = [idx_loc(idx)-buff_2, idx_loc(idx)+buff_2];
        % pull data
        i = expList(tt);
        y_all = grouped(i).fliesonfood(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, static_color,'filled')
        plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
    end
    
    save_figure(gcf,[fig_dir 'static temp holds'],fig_type, true);
end

% ALL OTHERS COMPARISON
if others
   for exp = 1:length(plot_list) 
    fig = getfig('',1); hold on
        %empty axis:
        formatFig(fig, blkbgd);
        xlim(x_limits)
        ylim(y_limits)
        
        ylabel('flies on food (#)')
        xlabel('temperature (\circC)')
        ax = gca;
        set(ax,'XTick', 1:length(x_tick_label), "XTickLabel", x_tick_label)
        ax.FontSize = 18;
    
        ii = plot_list(exp);
        % 1) FULL EXP WARMING AND COOLING
        idx = 0;
        % Cooling first
        temp_list = flip(dynamic_temps);
        for tt = 1:length(dynamic_temps)
            idx = idx +1;
        
            %  frame locations
            frames = grouped(ii).position.loc(1,temp_list(tt)).frames;
            x_roi = [idx-buff_1, idx+buff_1];
            x_line = [idx-buff_2, idx+buff_2];
        
            % proximity to food for dynamic trials
            y_all = grouped(ii).fliesonfood(frames,:);
            y_mean = mean(y_all, 1,'omitnan');
            y_line = mean(y_mean);
            x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
            scatter(x,y_mean, sz, grouped(ii).color,'filled')
            plot(x_line,[y_line,y_line],'color', grouped(ii).color, 'linewidth', LW,'LineStyle','-')
        end
        
        % Warming second
        temp_list = (dynamic_temps);
        for tt = 1:length(dynamic_temps)
            idx = idx +1;
            %  frame locations
            frames = grouped(ii).position.loc(3,temp_list(tt)).frames;
            x_roi = [idx-buff_1, idx+buff_1];
            x_line = [idx-buff_2, idx+buff_2];
        
            % proximity to food for dynamic trials
            y_all = grouped(ii).fliesonfood(frames,:);
            y_mean = mean(y_all, 1,'omitnan');
            y_line = mean(y_mean);
            x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
            scatter(x,y_mean, sz, grouped(ii).color,'filled')
            plot(x_line,[y_line,y_line],'color', grouped(ii).color, 'linewidth', LW,'LineStyle','-')
        end

        save_figure(fig,[fig_dir grouped(ii).name],fig_type, true);
    
    end
end

%% FIGURE: NOT DONE Speed
%   close all
clearvars('-except',initial_vars{:})
fig_dir = [saveDir '/PPT FIgures/Speed/'];
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end
close_figures = false;
auto_save = false;
plot_static = 1;
others = 1;

% INDEX LIST
caviar_idx = 4;
waxed_idx = 5;
water_idx = 7;
stat_idx = [1,2,3,6]; % 17, 20, 23, 25
empty_idx = 8;
AN_MP_idx = 9;
arista_idx = 10;

% plot multioptions
plot_list = 1:num.exp;
plot_list(stat_idx) = [];

% COLOR LIST
static_color = Color('gold');

% ploting parameters
buff_1 = 0.25;
buff_2 = 0.5;
sz = 40;
LW = 3;
x_limits = [0,19];
y_limits = [0, 4];
% [foreColor,backColor] = formattingColors(blkbgd); %get background colors

% experiment indexes:
static_temps = [3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
% static_temps = [19,15,9, 3, 3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
static_temp_actual = grouped(caviar_idx).position.temp_list(static_temps);

dynamic_temps = 3:2:19; % all whole number temps
temp_actual = grouped(caviar_idx).position.temp_list(dynamic_temps);
temp_names_label = [flip(temp_actual),temp_actual];
for i = 1:length(temp_names_label)
    x_tick_label{i}  = num2str(temp_names_label(i));
end

fig = getfig('',1); hold on
    %empty axis:
    formatFig(fig, blkbgd);
    xlim(x_limits)
    % ylim(y_limits)

    ylabel('proximity to food (mm)')
    xlabel('temperature (\circC)')
    ax = gca;
    set(ax,'XTick', 1:length(x_tick_label), "XTickLabel", x_tick_label)
    ax.FontSize = 18;

% save_figure(gcf,[fig_dir 'empty frame'],fig_type, false, false);

% ==============================================================================
% STATIC HOLD TEMPERATURES
if plot_static
    idx_loc = [1,3,6,9,10,13,16,18]; % x positions of the selected temperatures
    idx = 0;
    % Cooling first
    temp_list = flip(static_temps);
    expList = flip(stat_idx);
    for tt = 1:length(static_temps)
        idx = idx +1; 
        %  frame locations
        frames = grouped(4).position.loc(1,temp_list(tt)).frames; % 1=cooling, 3=warming
        x_roi = [idx_loc(idx)-buff_1, idx_loc(idx)+buff_1];
        x_line = [idx_loc(idx)-buff_2, idx_loc(idx)+buff_2];
        % pull data
        i = expList(tt);
        y_all = grouped(i).speed.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, static_color,'filled')
        plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
    end
    
    % warming second
    temp_list = (static_temps);
    expList = (stat_idx);
    for tt = 1:length(static_temps)
        idx = idx +1; 
        %  frame locations
        frames = grouped(4).position.loc(3,temp_list(tt)).frames; % 1=cooling, 3=warming
        x_roi = [idx_loc(idx)-buff_1, idx_loc(idx)+buff_1];
        x_line = [idx_loc(idx)-buff_2, idx_loc(idx)+buff_2];
        % pull data
        i = expList(tt);
        y_all = grouped(i).speed.all(frames,:);
        y_mean = mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, static_color,'filled')
        plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
    end
    
    save_figure(gcf,[fig_dir 'static temp holds'],fig_type, auto_save, close_figures);
end

% ALL OTHERS COMPARISON
if others
   for exp = 1:length(plot_list) 
    fig = getfig('',1); hold on
        %empty axis:
        formatFig(fig, blkbgd);
        xlim(x_limits)
        ylim(y_limits)

        ylabel('proximity to food (mm)')
        xlabel('temperature (\circC)')
        ax = gca;
        set(ax,'XTick', 1:length(x_tick_label), "XTickLabel", x_tick_label)
        ax.FontSize = 18;
    
        ii = plot_list(exp);
        % 1) FULL EXP WARMING AND COOLING
        idx = 0;
        % Cooling first
        temp_list = flip(dynamic_temps);
        for tt = 1:length(dynamic_temps)
            idx = idx +1;
        
            %  frame locations
            frames = grouped(ii).position.loc(1,temp_list(tt)).frames;
            x_roi = [idx-buff_1, idx+buff_1];
            x_line = [idx-buff_2, idx+buff_2];
        
            % proximity to food for dynamic trials
            y_all = grouped(ii).speed.all(frames,:);
            y_mean = mean(y_all, 1,'omitnan');
            y_line = mean(y_mean);
            x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
            scatter(x,y_mean, sz, grouped(ii).color,'filled')
            plot(x_line,[y_line,y_line],'color', grouped(ii).color, 'linewidth', LW,'LineStyle','-')
        end
        
        % Warming second
        temp_list = (dynamic_temps);
        for tt = 1:length(dynamic_temps)
            idx = idx +1;
            %  frame locations
            frames = grouped(ii).position.loc(3,temp_list(tt)).frames;
            x_roi = [idx-buff_1, idx+buff_1];
            x_line = [idx-buff_2, idx+buff_2];
        
            % proximity to food for dynamic trials
            y_all = grouped(ii).speed.all(frames,:);
            y_mean = mean(y_all, 1,'omitnan');
            y_line = mean(y_mean);
            x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
            scatter(x,y_mean, sz, grouped(ii).color,'filled')
            plot(x_line,[y_line,y_line],'color', grouped(ii).color, 'linewidth', LW,'LineStyle','-')
        end

        save_figure(fig,[fig_dir grouped(ii).name],fig_type, auto_save, close_figures);
    
    end
end








%% FIGURE: Correlation -  flies on food vs distance to food

clearvars('-except',initial_vars{:})
fig_dir = [saveDir '/PPT FIgures/Correlation/'];
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end
close_figures = false;
auto_save = false;
plot_static = 1;
others = 1;

x_limits = [5,35];
y_limits = [0,4];

% INDEX LIST
caviar_idx = 4;
waxed_idx = 5;
water_idx = 7;
stat_idx = [1,2,3,6]; % 17, 20, 23, 25
empty_idx = 8;
AN_MP_idx = 9;
arista_idx = 10;

% plot multioptions
plot_list = 1:num.exp;
plot_list(stat_idx) = [];

% COLOR LIST
static_color = Color('gold');

% ploting parameters
sz = 40;
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

% experiment indexes:
static_temps = [3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
% static_temps = [19,15,9, 3, 3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
static_temp_actual = grouped(caviar_idx).position.temp_list(static_temps);

dynamic_temps = 3:2:19; % all whole number temps
temp_actual = grouped(caviar_idx).position.temp_list(dynamic_temps);
% temp_names_label = [flip(temp_actual),temp_actual];
% for i = 1:length(temp_names_label)
%     x_tick_label{i}  = num2str(temp_names_label(i));
% end

 [rho,pval] = deal([]);

% ALL OTHERS COMPARISON
if others
   for exp = 1:length(plot_list) 
     x = []; y = [];
    fig = getfig('',1); hold on
        %empty axis:
        formatFig(fig, blkbgd);
        xlim(x_limits)
        ylim(y_limits)

        ylabel('flies on food (#)')
        xlabel('proximity to food (mm)')
        ax = gca;
        set(gca,'xdir', 'reverse','FontSize', 18)
            
        ii = plot_list(exp);
        % 1) FULL EXP WARMING AND COOLING
        idx = 0;
        % Cooling first
        temp_list = flip(dynamic_temps);
        for tt = 1:length(dynamic_temps)
            idx = idx +1;
        
            % cooling  frame locations
            frames = grouped(ii).position.loc(1,temp_list(tt)).frames;
            x_all = grouped(ii).dist.all(frames,:);
            x_mean = mean(x_all, 1,'omitnan');
            y_all = grouped(ii).fliesonfood(frames,:);
            y_mean = mean(y_all, 1,'omitnan');
            scatter(x_mean,y_mean, sz, grouped(ii).color)
            x = [x, x_mean]; y = [y, y_mean];
            
             % warming  frame locations
            frames = grouped(ii).position.loc(3,temp_list(tt)).frames;
            x_all = grouped(ii).dist.all(frames,:);
            x_mean = mean(x_all, 1,'omitnan');
            y_all = grouped(ii).fliesonfood(frames,:);
            y_mean = mean(y_all, 1,'omitnan');
            scatter(x_mean,y_mean, sz, grouped(ii).color,'filled')
            x = [x, x_mean]; y = [y, y_mean];
        end

        [rho(exp),pval(exp)] = corr(x',y');

        save_figure(fig,[fig_dir grouped(ii).name],fig_type, auto_save, close_figures);
    
    end
end


% STATS FIGURE:

fig = getfig('',1); hold on

for exp = 1:length(plot_list) 
   i = plot_list(exp);
    scatter(exp, rho(exp), 100, grouped(i).color,'filled')
end
xlim([0,7])
formatFig(fig, blkbgd);
set(gca,'ydir', 'reverse','xcolor', backColor)
xlabel(' ')
ylabel('Spearmans correlation')

save_figure(fig,[fig_dir 'proximity vs flies on food correlation'],fig_type, auto_save, close_figures);





 save_figure(fig,[saveDir grouped(i).name ' flies on food'],'-png',false,false);



%% FIGURE: occupancy of food region for static and dynamic W & C separated 
clearvars('-except',initial_vars{:})
fig_dir = [saveDir '/Poster Figures/'];
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end
buff_1 = 0.5;
buff_2 = 0.75;
sz = 60;
LW = 2;
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

temps_idx = [3, 9, 15, 19];  % 17,20,23,25 [temp locations that match the assigned hold temps]
temp_actual = grouped(4).position.temp_list(temps_idx);

dyn_exp = [4]; %start with just the caviar trial
stat_idx = [1,2,3,6]; % experiment index for the static trials

fig = getfig('',1,[350 689]); hold on

for tt = 1:length(temps_idx)
    %  frame locations
    frames_up = grouped(4).position.loc(3,temps_idx(tt)).frames;
    frames_down = grouped(4).position.loc(1,temps_idx(tt)).frames;
    x_roi = [temp_actual(tt)-buff_1, temp_actual(tt)+buff_1];
    x_line = [temp_actual(tt)-buff_2, temp_actual(tt)+buff_2];

    % %proximity to food for dynamic trials
    % for i = dyn_exp
    %     % warming
    %     y_all = grouped(i).occ.all(frames_up,:)*100;
    %     y_mean = mean(y_all, 1,'omitnan');
    %     y_line = mean(y_mean);
    %     x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
    %     % scatter(x,y_mean, sz, grouped(i).color,'filled')
    %     % plot(x_line,[y_line,y_line],'color', grouped(i).color, 'linewidth', LW,'LineStyle','-')
    %     scatter(x,y_mean, sz, Color('red'),'filled')
    %     plot(x_line,[y_line,y_line],'color', Color('red'), 'linewidth', LW)
    %     % cooling
    %     y_all = grouped(i).occ.all(frames_down,:)*100;
    %     y_mean = mean(y_all, 1,'omitnan');
    %     y_line = mean(y_mean);
    %     % scatter(x,y_mean, sz, grouped(i).color)
    %     % plot(x_line,[y_line,y_line],'color', grouped(i).color, 'linewidth', LW,'LineStyle','--')
    %     scatter(x,y_mean, sz, Color('dodgerblue'), 'filled')
    %     plot(x_line,[y_line,y_line],'color', Color('dodgerblue'), 'linewidth', LW)
    % end
    % 
    % plot concordant static trial
    static_color = Color('teal');
    i = stat_idx(tt);
    y_all = grouped(i).occ.all([frames_up; frames_down],:)*100;
    y_mean = mean(y_all, 1,'omitnan');
    y_line = mean(y_mean);
    x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
    % scatter(x,y_mean, sz, static_color,'filled')
    % plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
    scatter(x,y_mean, sz, Color('black'),'filled')
    plot(x_line,[y_line,y_line],'color', Color('black'), 'linewidth', LW,'LineStyle','-')
end
% 
% Figure formatting
formatFig(fig, blkbgd);
h_line(14.4, 'grey',':',0.5) 
xlim([16,26])
ylim([0,100])
ylabel('food region occupancy (%)')
xlabel('temperature (\circC)')
ax = gca;
set(gca,'tickdir','out','XTick',temp_actual)
ax.FontSize = 18;

% save_figure(fig,[fig_dir 'occupancy ' grouped(i).name],fig_type);
save_figure(fig,[fig_dir 'occupancy static caviar' ],fig_type);

%% FIGURE: speed for static and dynamic W & C separated 
clearvars('-except',initial_vars{:})
buff_1 = 0.5;
buff_2 = 0.75;
sz = 60;
LW = 2;
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

temps_idx = [3, 9, 15, 19];  % 17,20,23,25 [temp locations that match the assigned hold temps]
temp_actual = grouped(4).position.temp_list(temps_idx);

dyn_exp = [5]; %start with just the caviar trial
stat_idx = [1,2,3,6]; % experiment index for the static trials

fig = getfig('',1,[350 689]); hold on
for tt = 1:length(temps_idx)
    %  frame locations
    frames_up = grouped(4).position.loc(3,temps_idx(tt)).frames;
    frames_down = grouped(4).position.loc(1,temps_idx(tt)).frames;
    x_roi = [temp_actual(tt)-buff_1, temp_actual(tt)+buff_1];
    x_line = [temp_actual(tt)-buff_2, temp_actual(tt)+buff_2];

    % %proximity to food for dynamic trials
    % for i = dyn_exp
    %     % warming
    %     y_all = grouped(i).speed.all(frames_up,:);
    %     y_mean = mean(y_all, 1,'omitnan');
    %     y_line = mean(y_mean);
    %     x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
    %     % scatter(x,y_mean, sz, grouped(i).color,'filled')
    %     % plot(x_line,[y_line,y_line],'color', grouped(i).color, 'linewidth', LW,'LineStyle','-')
    %     scatter(x,y_mean, sz, Color('red'),'filled')
    %     plot(x_line,[y_line,y_line],'color', Color('red'), 'linewidth', LW,'LineStyle','-')
    %     % cooling
    %     y_all = grouped(i).speed.all(frames_down,:);
    %     y_mean = mean(y_all, 1,'omitnan');
    %     y_line = mean(y_mean);
    %     % scatter(x,y_mean, sz, grouped(i).color)
    %     % plot(x_line,[y_line,y_line],'color', grouped(i).color, 'linewidth', LW,'LineStyle','--')
    %     scatter(x,y_mean, sz, Color('dodgerblue'),'filled')
    %     plot(x_line,[y_line,y_line],'color', Color('dodgerblue'), 'linewidth', LW)
    % end

    % plot concordant static trial
    static_color = Color('black');
    i = stat_idx(tt);
    y_all = grouped(i).speed.all([frames_up; frames_down],:);
    y_mean = mean(y_all, 1,'omitnan');
    y_line = mean(y_mean);
    x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
    scatter(x,y_mean, sz, static_color,'filled')
    plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
end
% 
% Figure formatting
formatFig(fig, blkbgd);

xlim([16,26])
ylim([0,9])
ylabel('fly movement speed (mm/s)')
xlabel('temperature (\circC)')
ax = gca;
set(gca,'tickdir','out','XTick',temp_actual)
ax.FontSize = 18;

% save_figure(fig,[saveDir 'static vs dynamic occupancy caviar and waxed'],fig_type);
save_figure(fig,[saveDir 'scatter speed static temps'],fig_type);

%% FIGURE: eccentricity for static and dynamic W & C separated 
clearvars('-except',initial_vars{:})
buff_1 = 0.5;
buff_2 = 0.75;
sz = 60;
LW = 2;
[foreColor,backColor] = formattingColors(blkbgd); %get background colors

temps_idx = [3, 9, 15, 19];  % 17,20,23,25 [temp locations that match the assigned hold temps]
temp_actual = grouped(4).position.temp_list(temps_idx);

dyn_exp = [5]; %start with just the caviar trial
stat_idx = [1,2,3,6]; % experiment index for the static trials

arena_radius = (data(1).data(1).data.r/pix2mm); 

fig = getfig('',1,[350 689]); hold on
for tt = 1:length(temps_idx)
    %  frame locations
    frames_up = grouped(4).position.loc(3,temps_idx(tt)).frames;
    frames_down = grouped(4).position.loc(1,temps_idx(tt)).frames;
    x_roi = [temp_actual(tt)-buff_1, temp_actual(tt)+buff_1];
    x_line = [temp_actual(tt)-buff_2, temp_actual(tt)+buff_2];

    %proximity to food for dynamic trials
    for i = dyn_exp
        % warming
        y_all = grouped(i).ecent.all(frames_up,:);
        y_mean = arena_radius - mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        % scatter(x,y_mean, sz, grouped(i).color,'filled')
        % plot(x_line,[y_line,y_line],'color', grouped(i).color, 'linewidth', LW,'LineStyle','-')
        scatter(x,y_mean, sz, Color('red'),'filled')
        plot(x_line,[y_line,y_line],'color', Color('red'), 'linewidth', LW,'LineStyle','-')
        % cooling
        y_all = grouped(i).ecent.all(frames_down,:);
        y_mean = arena_radius - mean(y_all, 1,'omitnan');
        y_line = mean(y_mean);
        % scatter(x,y_mean, sz, grouped(i).color)
        % plot(x_line,[y_line,y_line],'color', grouped(i).color, 'linewidth', LW,'LineStyle','--')
        scatter(x,y_mean, sz, Color('dodgerblue'),'filled')
        plot(x_line,[y_line,y_line],'color', Color('dodgerblue'), 'linewidth', LW)
    end

    % % plot concordant static trial
    % static_color = Color('black');
    % i = stat_idx(tt);
    % y_all = grouped(i).ecent.all([frames_up; frames_down],:);
    % y_mean = arena_radius - mean(y_all, 1,'omitnan');
    % y_line = mean(y_mean);
    % x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
    % scatter(x,y_mean, sz, static_color,'filled')
    % plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
end
% 
% Figure formatting
formatFig(fig, blkbgd);

xlim([16,26])
ylim([6,20])
ylabel('distance to edge (mm)')
xlabel('temperature (\circC)')
ax = gca;
set(gca,'tickdir','out','XTick',temp_actual)
ax.FontSize = 18;

% save_figure(fig,[saveDir 'static vs dynamic occupancy caviar and waxed'],fig_type);
save_figure(fig,[saveDir 'scatter eccentricity ' grouped(i).name],fig_type);

 %% ANALYSIS AND FIGURE: Sleep quantity per fly for ENTIRE experiment
clearvars('-except',initial_vars{:})

tP = getTempTurnPoints('linear_ramp_F_25-17');
ROI = min(tP.HoldROI):max(tP.HoldROI);
SZ = 15;
idx = [4 1 2 3 6 5];

fig = getfig('', 1, [472 577]); hold on
    for i = 1:length(idx)
        exp = idx(i);
        y = sleep(exp).avg_quant;
        y_err = mean(y, 'omitnan');
    
        scatter(i*ones([1,length(y)]), y, SZ, grouped(exp).color, 'filled', 'XJitter', 'density')
        plot([i-0.4,i+0.4],[y_err, y_err],'color', grouped(exp).color,'linewidth', 2)
    end
ylabel('Sleep Quantity (min/fly)')
typeLabel = {grouped(idx).name}';
formatFig(fig, blkbgd);
set(gca, 'xtick', 1:length(idx), 'xticklabel', typeLabel)

save_figure(fig,[saveDir 'avg sleep per fly'],fig_type);


%% FIGURE: probabiltiy of sleep in each temperature regime
% uses Bayes' Formula to calculate the probability of sleep events
% Formula: P(W|S) = P(S|W)P(W) / (P(S|W)P(W) + P(S|C)P(C) + P(S|H)P(H))
% aka: probability of it being 'warming' if a fly is sleeping

clearvars('-except',initial_vars{:})
buffer = 0.15;
[foreColor,~] = formattingColors(blkbgd);

SZ = 50;

for i = 1:num.exp
    plotData = [];
    all_roi = [];
    tPoints = getTempTurnPoints('linear_ramp_F_25-17');
    
    roiC = tPoints.DownROI;
    roiW = tPoints.UpROI;
    roiH = tPoints.HoldROI;
    all_roi =[roiC,roiW,roiH];
            
    % variables needed to know:
    pW = length(roiW)/length(all_roi); % prob of warming
    pC = length(roiC)/length(all_roi); %prob of cooling
    pH = length(roiH)/length(all_roi); %prob of holding
    
    pSgW = sum(sleep(i).num(roiW,:))/length(all_roi); %prob of sleep given warming condition
    pSgC = sum(sleep(i).num(roiC,:))/length(all_roi); %prob of sleep given cooling condition
    pSgH = sum(sleep(i).num(roiH,:))/length(all_roi); %prob of sleep given holding condition

    % Calculate the probabilities
    probSW = pSgW.*pW ./ ( pSgW.*pW + pSgC.*pC + pSgH.*pH);
    probSC = pSgC.*pC ./ ( pSgW.*pW + pSgC.*pC + pSgH.*pH);
    probSH = pSgH.*pH ./ ( pSgW.*pW + pSgC.*pC + pSgH.*pH);
    
    plotData = [probSC',probSW',probSH']; %cool, warm, hold ploting
    
    % plot the probabilities
    fig = getfig('',1,[481 680]); hold on
     for type = 1:3
        switch type
            case 1 %decreasing
                kolor = Color('dodgerblue');
            case 2 % increasing
                kolor = Color('red');
            case 3 % holding
                kolor = Color('grey');
        end
        %average
        bar(type, mean(plotData(:,type),'omitnan'),'FaceColor',kolor,'FaceAlpha',1,'EdgeColor',foreColor)
        %scatterpoints
        x = shuffle_data(linspace(type-buffer,type+buffer,num.trial(i)));
        scatter(x,plotData(:,type),SZ,foreColor)
     end

    % formatting and labels
    set(gca,'XTick',1:3,'XTickLabel',{'cooling','heating','hold'})
    ylabel('Prob of sleep for each temp condition')
    xlim([0,4])
    ylim([0,1])
    formatFig(fig,blkbgd);
    title(grouped(i).name)
    % Save figure
    save_figure(fig,[saveDir 'Sleep\' expNames{i} 'Bayes probabilty of sleeping by temp type'],fig_type);
end

%% FIGURE: Normalized sleeping over time
plot_err = false;
[~,backColor] = formattingColors(blkbgd); %get background colors
LW = 1.5;
sb(1).idx = 1;
sb(2).idx = 2:3;
r = 3;
c = 1;
sSpan = 90; %30 second smoothing
idx = [4 1 2 3 6 5];
% idx = [1 2 3 6];

fig = getfig('',1);
for i = idx %1:num.exp

    subplot(r,c,sb(1).idx); hold on
        time = grouped(i).time;
        temp = grouped(i).temp;
        % kolor = grouped(i).color;
        switch i 
            case 4 % intact & dynamic caviar
                kolor = Color('dodgerblue');
            case 1 % intact & dynamic caviar
                kolor = Color('black');
            case 2 % intact & dynamic caviar
                kolor = Color('darkslategray');
           case 3 % intact & dynamic caviar
                kolor = Color('dimgray');
            case 6 % intact & dynamic caviar
                kolor = Color('silver');
            case 5 % intact & dynamic caviar
                kolor = Color('gold');
        end
        plot(time, temp, 'color',kolor, 'linewidth', LW)
        ylabel('temp (\circC)')
        xlim([0,365])
    subplot(r,c,sb(2).idx); hold on
        y = smooth(sum(sleep(i).num,2),sSpan, 'moving');
        y_norm = zscore(y);
        plot(time,y_norm,'color',kolor,'linewidth',LW)
        ylabel('sleeping (zscore)')
        xlabel('time (min)')
        xlim([0,365])
end

formatFig(fig,blkbgd,[r,c],sb);
subplot(r,c,sb(1).idx);
set(gca,'xcolor',backColor) 

save_figure(fig,[saveDir 'Sleep\' expGroup ' normalized sleep timecourse'],fig_type);


%%

%% FIGURE: animation style overview 
%   close all
fig_dir = [saveDir '/PPT FIgures/Quadrant Occupancy/'];
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end

plot_caviar = 0;
plot_static = 1;
plot_water = 0;
plot_waxed = 0;
plot_empty = 0;
plot_aristaless = 0;
plot_AN_MPless = 0;

plotList = [plot_static,plot_caviar,plot_water,plot_waxed,plot_empty,plot_aristaless,plot_AN_MPless];

% INDEX LIST IN STRUCTURE FOR THE EXP GROUPS
caviar_idx = 4;
waxed_idx = 5;
water_idx = 7;
stat_idx = [1,2,3,6]; % for temps 17, 20, 23, 25 
empty_idx = 8;
AN_MP_idx = 9;
arista_idx = 10;
dynm_idx = 1:num.exp; dynm_idx(stat_idx) = [];
% dynm_idx = [9, 12, 10];

% plot multioptions
plot_list = 1:num.exp;
plot_list(stat_idx) = [];

% COLOR LIST
static_color = Color('gold');
grouped(water_idx).color = Color('white');
grouped(waxed_idx).color = Color('lime');
grouped(caviar_idx).color = Color('dodgerblue');
grouped(empty_idx).color = Color('navajowhite');
grouped(arista_idx).color = Color('DarkOrchid');
grouped(AN_MP_idx).color = Color('red');

% ploting parameters
buff_1 = 0.25;
buff_2 = 0.5;
sz = 40;
LW = 3;
font_size  = 25;
tickDir = 'in';
x_limits = [0,19];
y_limits = [0 100];
% [foreColor,backColor] = formattingColors(blkbgd); %get background colors

% experiment indexes:
static_temps = [3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
% static_temps = [19,15,9, 3, 3, 9, 15, 19]; % 17,20,23,25 [temp locations that match the assigned hold temps]
static_temp_actual = grouped(caviar_idx).position.temp_list(static_temps);

dynamic_temps = [3,9,15,19];%3:2:19; % all whole number temps
temp_actual = grouped(caviar_idx).position.temp_list(dynamic_temps);
temp_names_label = [flip(temp_actual),temp_actual];
xtick_loc = [1 3 6 9 10 13 16 18];
for i = 1:length(temp_names_label)
    x_tick_label{i}  = num2str(temp_names_label(i));
end


% ============= PLOT STATIC ===================
fig = getfig('',1); hold on
    %empty axis:
    formatFig(fig, blkbgd);
    xlim(x_limits)
    ylim(y_limits)
    set(gca, 'ytick', y_limits(1):20:y_limits(2))
    % h_line(25, 'grey',':',0.5) 
    
    ylabel('flies in quadrant (%)')
    xlabel('temperature (\circC)')
    ax = gca;
    % set(ax,'ydir','reverse')
    set(ax,'XTick', xtick_loc, "XTickLabel", x_tick_label)
    ax.FontSize = font_size;
    ax.TickDir = tickDir;
% save_figure(gcf,[fig_dir 'static vs dynamic distance to food empty frame'],fig_type, false, false);

% STATIC HOLD TEMPERATURES

idx_loc = [1,3,6,9,10,13,16,18]; % x positions of the selected temperatures
idx = 0;
% Cooling first
temp_list = flip(static_temps);
expList = flip(stat_idx);
for tt = 1:length(static_temps)
    idx = idx +1; 
    %  frame locations
    x_roi = [idx_loc(idx)-buff_1, idx_loc(idx)+buff_1];
    x_line = [idx_loc(idx)-buff_2, idx_loc(idx)+buff_2];
    % pull data
    i = expList(tt);
    y_mean = grouped(i).quadrant.decreasing.raw(temp_list(tt),:);
    y_line = mean(y_mean);
    x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
    scatter(x,y_mean, sz, static_color,'filled')
    plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
end

% warming second
temp_list = (static_temps);
expList = (stat_idx);
for tt = 1:length(static_temps)
    idx = idx +1; 
    x_roi = [idx_loc(idx)-buff_1, idx_loc(idx)+buff_1];
    x_line = [idx_loc(idx)-buff_2, idx_loc(idx)+buff_2];
    % pull data
    i = expList(tt);
    y_mean = grouped(i).quadrant.increasing.raw(temp_list(tt),:);
    y_line = mean(y_mean);
    x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
    scatter(x,y_mean, sz, static_color,'filled')
    plot(x_line,[y_line,y_line],'color', static_color, 'linewidth', LW,'LineStyle','-')
end

save_figure(gcf,[fig_dir 'static vs dynamic distance to food static'],fig_type, false, true);


% ============= PLOT DYNAMIC ===================

for exp = dynm_idx
    % MAKE FIGURE
    fig = getfig('',1); hold on
    %empty axis:
    formatFig(fig, blkbgd);
    xlim(x_limits)
    ylim(y_limits)
    set(gca, 'ytick', y_limits(1):20:y_limits(2))
    % h_line(25, 'grey',':',0.5) 
    
    ylabel('flies in quadrant (%)')
    xlabel('temperature (\circC)')
    ax = gca;
    % set(ax,'ydir','reverse')
    set(ax,'XTick', xtick_loc, "XTickLabel", x_tick_label)
    ax.FontSize = font_size;
    ax.TickDir = tickDir;

    % PLOT DATA
    kolor = grouped(exp).color;
    idx = 0;
    % Cooling first
    temp_list = flip(dynamic_temps);
    for tt = 1:length(dynamic_temps)
        idx = idx +1;
        x_loc = xtick_loc(idx);
        x_roi = [x_loc-buff_1, x_loc+buff_1];
        x_line = [x_loc-buff_2, x_loc+buff_2];
        % proximity to food for dynamic trials
        y_mean = grouped(exp).quadrant.decreasing.raw(temp_list(tt),:);
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, kolor,'filled')
        plot(x_line,[y_line,y_line],'color', kolor, 'linewidth', LW,'LineStyle','-')
    end    
    % Warming second
    temp_list = (dynamic_temps);
    for tt = 1:length(dynamic_temps)
        idx = idx +1;
        x_loc = xtick_loc(idx);
        x_roi = [x_loc-buff_1, x_loc+buff_1];
        x_line = [x_loc-buff_2, x_loc+buff_2];
        % proximity to food for dynamic trials
        y_mean = grouped(exp).quadrant.increasing.raw(temp_list(tt),:);
        y_line = mean(y_mean);
        x = shuffle_data(linspace(x_roi(1), x_roi(2), length(y_mean)));
        scatter(x,y_mean, sz, kolor,'filled')
        plot(x_line,[y_line,y_line],'color', kolor, 'linewidth', LW,'LineStyle','-')
    end
        
    save_figure(gcf,[fig_dir 'static vs dynamic distance to food ' grouped(exp).name],fig_type);
end









% [9, 12, 10] %



















