 
% need to run prep code 'normalization distance' in COM 
clearvars('-except',initial_vars{:})


%% FIGURE: compare the differences in sleeping  
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

%% FIGURE: animation style overview #imextra
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











