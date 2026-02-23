
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors
clc

%% Speed between regions timecourse

% pull out when flies in outer ring = 1
% pull out speed

clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

% Create lists to cycle through
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color
nRegions = length(region_name);
r = 3; c = 1;
sb(1).idx = 1;
sb(2).idx = 2:3; 

% FIGURE
fig = getfig('',true,[1450, 900]);
subplot(r,c,sb(1).idx)
plot(data.time,data.temperature,'color', foreColor,'LineWidth', 2)
    ylabel('temp (\circC)')
subplot(r,c,sb(2).idx)
hold on
% Pull out fly speed in each region
for ii = 1:nRegions
    current_var = data.(region_name{ii});
    loc = logical(replaceNaN(current_var,0));    
    speed = data.speed;
    % Replace speed values when fly is not in the wanted region with nans
    speed(~loc) = nan;
    % Remove M/F fly dimension - all flies and their speeds together
    speed_all = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
    % Calculate mean
    m = mean(speed_all,2,'omitnan');    
    % Plot average speed in each region over time
    plot(data.time,smooth(m,300,'moving'),'Color',Color(color_list{ii}))
end

% Format fig
formatFig(fig,blkbgd,[r,c], sb);
subplot(r,c,sb(1).idx)
    set(gca, 'xcolor', 'none')
subplot(r,c,sb(2).idx)
legend(region_name,"TextColor",foreColor,"Color",~foreColor,"EdgeColor",~foreColor)
ylabel('average speed (mm/s)')
xlabel('time (min)')

% Save figure
save_figure(fig,[figDir 'speed between regions timecourse'],fig_type);

%% Speed between regions within temp regimes scatter plot

% pull out when flies are in each region
% pull out speed for each temp regime time period

clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colo

% Create lists to cycle through
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
regime_name = {'WT','WS','CT','CS'}; % temp regime
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color

% Establish sizes for loops
nRegions = length(region_name);
nRegimes = length(regime_name);

% Set up x axis placement
h = []; 
x = [];
for region = 1:nRegions % 3
    h = region;
    for regime = 1:nRegimes % 4
        h = [h,h(end)+nRegimes];   
    end  
    x(region,:) = h(:,1:nRegimes);
end

% Create list for x axis labels and locations
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
x_labels = [];
for regime = 1:nRegimes % 4
    x_labels = [x_labels,median(x(:,regime))];
end

[y_avg, y_sem] = deal(nan([4,1]));

% FIGURE
sz = 50;
buff = 0.2;
fig = getfig('',true,[1450, 900]);
hold on
% Pull out fly speed in each region
for region = 1:nRegions % 3
    current_var = data.(region_name{region});
    loc = logical(replaceNaN(current_var,0));    
    speed = data.speed;
    % Replace speed values when fly is not in the wanted region with nans
    speed(~loc) = nan;
    % Remove M/F fly dimension - all flies and their speeds together
    speed_all = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
    % Pull out fly speed in each temp regime
    for regime = 1:nRegimes % 4
        current_reg = data.tempbin.(regime_name{regime});
        loc2 = replaceNaN(current_reg,0); % TODO can be deleted
        speed2 = speed_all; % M and F combined speed matrix copy
        % Replace speed values not during the wanted te mp regime with nans
        speed2(~loc2,:) = nan;
        % Calculate mean
        avgspeed = mean(speed2,1,'omitnan');  
        y_avg = median(avgspeed,2,'omitnan');
        y_sem = std(avgspeed,0,2,'omitnan')./sqrt(num.trials);
        % Plot average speed in each region during each temp regime
        X = x(region,regime) * ones(size(avgspeed));
        scatter(X,avgspeed,sz,Color(color_list{region}),'filled', 'XJitter','density','XJitterWidth',0.5);
        errorbar(x(region,regime), y_avg, y_sem, 'color', foreColor, 'linestyle', 'none', 'LineWidth', 2,'Marker','_','MarkerSize',30);
    end
end

% Format figure
formatFig(fig,blkbgd);
% y axis specs
ylabel('average speed per fly (mm/s)')
ydim = [.4, .37, .34]; %hard coded for data distribution on figure
ylim([0 20])
yticks(0:5:20)
% x axis specs
xticks(x_labels)
xticklabels(labels)
% annotate on legend to designate regions by color (this might be janky but it works)
for region = 1:nRegions
annotation('textbox',[0.75,ydim(region),.5,.5],'String',region_name{region},'FitBoxToText','on','Color',...
    Color(color_list{region}),'FontSize',15,'EdgeColor',~foreColor)
end

% Save figure
save_figure(fig,[figDir 'speed between regions within temp regimes scatter plot'],fig_type);


%% Speed between temp regimes within regions scatter plot

clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colo

% Create lists to cycle through
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
regime_name = {'WT','WS','CT','CS'}; % temp regime
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color

% Establish sizes for loops
nRegions = length(region_name);
nRegimes = length(regime_name);

% Set up x axis placement
h = [];
x = [];
for regime = 1:nRegimes % 3
    h = regime;
    for region = 1:nRegions % 4
        h = [h,h(end)+nRegimes+1];
    end
    x(regime,:) = h(:,1:nRegions);
end

% Create list for x axis labels and locations
x_labels = [];
for region = 1:nRegions% 4
    x_labels = [x_labels;x(:,region)];
end

[y_avg, y_sem] = deal(nan([4,1]));
Rspeed = [];

% FIGURE
sz = 50;
buff = 0.2;
fig = getfig('',true,[1450, 900]);
hold on
% Pull out fly speed in each region
for region = 1:nRegions % 3
    current_var = data.(region_name{region});
    loc = logical(replaceNaN(current_var,0));    
    speed = data.speed;
    % Replace speed values when fly is not in the wanted region with nans
    speed(~loc) = nan;
    % Remove M/F fly dimension - all flies and their speeds together
    speed_all = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
    % Pull out fly speed in each temp regime
    for regime = 1:nRegimes % 4
        current_reg = data.tempbin.(regime_name{regime});
        loc2 = replaceNaN(current_reg,0); % TODO can be deleted
        speed2 = speed_all; % M and F combined speed matrix copy
        % Replace speed values not during the wanted te mp regime with nans
        speed2(~loc2,:) = nan;
        % Calculate mean
        avgspeed = mean(speed2,1,'omitnan');  
        y_avg = median(avgspeed,2,'omitnan');
        y_sem = std(avgspeed,0,2,'omitnan')./sqrt(num.trials);
        % Plot average speed in each region during each temp regime
        X = x(regime,region) * ones(size(avgspeed));
        scatter(X,avgspeed,sz,Color(color_list{region}),'filled', 'XJitter','density','XJitterWidth',0.5);
        errorbar(x(regime,region), y_avg, y_sem, 'color', foreColor, 'linestyle', 'none', 'LineWidth', 2,'Marker','_','MarkerSize',30);
        % plot(x(regime,region), y_avg, 'color', 'k')
        Rspeed = [Rspeed;avgspeed];
    end
    % Plot line connecting average speed for each fly across regimes
    for fly = 1:length(avgspeed)
        % plot(x(:,region),Rspeed(:,fly),'color',Color(color_list{region}))
    end
    % plot(x(:,region),mean(Rspeed),'color','k')
    Rspeed = []; 
end

% Format figure
formatFig(fig,blkbgd);
% y axis specs
ylabel('average speed per fly (mm/s)')
ydim = [0.4, 0.37, 0.34, 0.31]; %hard coded for data distribution on figure
ylim([0 20])
yticks(0:5:20)
% x axis specs
x_axis_labels = [];
for region = 1:nRegions
    x_axis_labels = [x_axis_labels,regime_name];
end
xticks(x_labels)
xticklabels(x_axis_labels)
% annotate on legend to designate regions by color (this might be janky but it works)
for region = 1:nRegions
annotation('textbox',[0.8,ydim(region),.5,.5],'String',region_name{region},'FitBoxToText','on','Color',...
    Color(color_list{region}),'FontSize',15,'EdgeColor',~foreColor)
end

% Save figure
save_figure(fig,[figDir 'speed between temp regimes within regions scatter plot'],fig_type);


%% Histogram of speed between regions

clearvars('-except',initial_var{:})
% Create lists to cycle through
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color
nRegions = length(region_name);
nBins = 25;

% FIGURE
fig = getfig;
hold on
% Pull out fly speed in each region
for ii = 1:nRegions
    current_var = data.(region_name{ii});
    loc = logical(replaceNaN(current_var,0));    
    speed = data.speed;
    % Replace speed values when fly is not in the wanted region with nans
    speed(~loc) = nan;
    % Remove M/F fly dimension - all flies and their speeds together
    speed_all = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
    histogram(speed_all,"FaceColor",Color(color_list{ii}),'FaceAlpha',0.5)
end

% Format figure
formatFig(fig,blkbgd);
% x axis specs
xlim([0 30]) 
xlabel('speed (mm/s)')
% y axis specs
ylim([0 (3*10^5)])
ylabel('frequency')
% annotate on legend to designate regions by color (this might be janky but it works)
ydim = [0.4, 0.37, 0.34, 0.31]; %hard coded for data distribution on figure
for region = 1:nRegions
annotation('textbox',[0.8,ydim(region),.5,.5],'String',region_name{region},'FitBoxToText','on','Color',...
    Color(color_list{region}),'FontSize',15,'EdgeColor',~foreColor)
end

% Save figure
save_figure(fig,[figDir 'speed between regions histogram'],fig_type);

%% Speed between regions within threat vs safe periods (not distinguishing warm vs cold)

% pull out when flies are in each region
% pull out speed for each temp regime time period

clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colo

% Create lists to cycle through
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
regime_name = {'WT','WS','CT','CS'}; % temp regime
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color

% Establish sizes for loops
nRegions = length(region_name);
nRegimes = 2; % hard coded for threat and safe only

% Set up x axis placement
h = []; 
x = [];
for region = 1:nRegions % 3
    h = region;
    for regime = 1:nRegimes % 4
        h = [h,h(end)+4];   % hard coded
    end  
    x(region,:) = h(:,1:nRegimes);
end



% Create list for x axis labels and locations
labels = {'Threat', 'Safe'};
x_labels = [];
for regime = 1:nRegimes % 4
    x_labels = [x_labels,median(x(:,regime))];
end

[y_avg, y_sem] = deal(nan([2,1]));

% FIGURE
sz = 50;
buff = 0.2;
fig = getfig('',true,[1450, 900]);
hold on
% Pull out fly speed in each region
for region = 1:nRegions % 3
    current_var = data.(region_name{region});
    loc = logical(replaceNaN(current_var,0));    
    speed = data.speed;
    % Replace speed values when fly is not in the wanted region with nans
    speed(~loc) = nan;
    % Remove M/F fly dimension - all flies and their speeds together
    speed_all = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
    % Pull out fly speed in each temp regime
    for regime = 1:nRegimes % 4
        switch regime
            case 1 % threat
                current_reg = data.tempbin.WT | data.tempbin.CT;
            case 2 % safe
                current_reg = data.tempbin.WS | data.tempbin.CS;
        end
        loc2 = replaceNaN(current_reg,0); % TODO can be deleted
        speed2 = speed_all; % M and F combined speed matrix copy
        % Replace speed values not during the wanted te mp regime with nans
        speed2(~loc2,:) = nan;
        % Calculate mean
        avgspeed = mean(speed2,1,'omitnan');  
        y_avg = median(avgspeed,2,'omitnan');
        y_sem = std(avgspeed,0,2,'omitnan')./sqrt(num.trials);
        % Plot average speed in each region during each temp regime
        X = x(region,regime) * ones(size(avgspeed));
        scatter(X,avgspeed,sz,Color(color_list{region}),'filled', 'XJitter','density','XJitterWidth',0.5);
        errorbar(x(region,regime), y_avg, y_sem, 'color', foreColor, 'linestyle', 'none', 'LineWidth', 2,'Marker','_','MarkerSize',30);
    end
end

% Format figure
formatFig(fig,blkbgd);
% y axis specs
ylabel('average speed per fly (mm/s)')
ydim = [.4, .37, .34]; %hard coded for data distribution on figure
ylim([0 20])
yticks(0:5:20)
% x axis specs
xticks(x_labels)
xticklabels(labels)
% annotate on legend to designate regions by color (this might be janky but it works)
for region = 1:nRegions
annotation('textbox',[0.75,ydim(region),.5,.5],'String',region_name{region},'FitBoxToText','on','Color',...
    Color(color_list{region}),'FontSize',15,'EdgeColor',~foreColor)
end

% Save figure
save_figure(fig,[figDir 'speed between regions within threat vs safe'],fig_type);


%% Select temp tuning curve for each region -- TODO: currently hard coded for speed only

clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

% Select the type of information to plot: 
[title_str,pName,y_dir,y_lab,nullD,scaler,dType,~,sexSep,ylimits] = ... 
    PlotParamSelectionHR('Heating and Cooling');
plot_err = true;

if isempty(title_str)
    return
end
fig_dir = [figDir, ' temp tuning curves/'];
% set figure folder
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir)
end
 
% ------------- Plotting parameters -------------

% select temp protocol specific plotting features
autoYLim = false; % Y LIMITS
if any(isnan(ylimits)) % in case new data is added without specs for axis limits
    autoYLim = true;
end

% TODO: update this as new protocols are added to the pipeline
if strcmp(groupName,'Berlin LTS caviar') % X LIMITS
    xlimits = [13, 37];
    autoXLim = false;
    xPos = [15, 25];  % for the shaded threat region on the plot
    data_cut_off = [15, 35];
else
    autoXLim = true;
end

region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color
nRegions = length(region_name);
r = 1; % figure rows
c = 2; % heating and cooling separated plots
LW = 2; % plotting line width
FA = 0.35; % SEM shading face alpha level

% FIGURE
fig = getfig('',1,[998 882]); % short and fat: [1230 637]
hold on
% pull larger data type group
for region = 1:nRegions
    current_var = data.(region_name{region});
    loc = logical(replaceNaN(current_var,0));    
    speed = data.(pName);
    % Replace speed values when fly is not in the wanted region with nans
    speed(~loc) = nan;
    % Remove M/F fly dimension - all flies and their speeds together
    yy = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
    
    % if sexSep % data separated per fly or per group of flies
    %     y_all = [squeeze(yy(:,M,:)), squeeze(yy(:,F,:))];
    % else 
    %     y_all = yy;
    % end
    x  = data.tempbin.temps; % temp bins
    nTemps = length(x); % number of temperature bins
    types  = {'cooling', 'warming'};
    
    % Extract and Plot data:
    for ii = 1:2
        subplot(r,c,ii); hold on
            Idx = data.tempbin.(types{ii});
            rawY = nan([nTemps, 2]); % first col = avg, second = sem
    
            % extract the tuning information across the temp bins
            for tt = 1:nTemps 
                % skip the temp bin if it's not one included in the protocol
                % (e.g. for cases where the temp slightly overshoots or
                % undershoots the value and then we get a 'read' for something
                % like the 14.5 bin when really there are a small number of
                % data points due to temp overshoot
                if x(tt)<data_cut_off(1) ||  x(tt)>data_cut_off(2)
                    continue
                end
                % all the cooling data across the flies that fits this temp bin
                y = yy(Idx(:,tt),:); 
    
                % fill the temp bin data into the appropriate structure
                raw = mean(y,'omitnan').*scaler; % find cooling fly data 
                rawY(tt,1) = mean(raw,'omitnan');
                rawY(tt,2) = sem(raw);
            end
    
            % plot the data: 
             plot_error_fills(plot_err, x, rawY(:,1), rawY(:,2), Color(color_list{region}), fig_type, FA);
             plot(x,rawY(:,1),'color', Color(color_list{region}), 'LineWidth', LW)
    end
end
         
% ------------ formatting ------------
formatFig(fig, blkbgd,[r,c]);
matchAxis(fig, true);
for ii = 1:2
    subplot(r,c,ii) 
    title(types{ii},'color', foreColor,'FontAngle','italic')
    xlabel('temperature (\circC)')
    if ~autoXLim;  xlim(xlimits);  end
    if ~autoYLim;  ylim(ylimits);  end
    set(gca, 'ydir', y_dir)
    % subplot specific adjustments
    if ii==1 % cooling
        set(gca, 'XDir','reverse')
        ylabel(y_lab)
    end
    if ii==2 % warming
        set(gca, 'YColor', 'none')
    end
    h_line(nullD, 'gray', '--',1)
end

% add shaded area for 'threat' temp region
for tt = 1:2
    subplot(r, c, tt)
    ylimits = ylim;
    pos = [xPos(1,tt), ylimits(1), range(xPos), range(ylimits)]; % [lower-left X, lower-left Y, X-width, Y-height]
    rectangle('Position', pos, 'FaceColor', foreColor, ... 
              'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

% add a time arrow to the temp axes graphs
arrow_h = 0.05;
arrow_x = [0.15, 0.18];
arrow_y = [arrow_h, arrow_h];
annotation('textarrow', arrow_x, arrow_y, 'String', 'time ','Color',foreColor,'FontSize',12);
arrow_x = [0.60, 0.63];
arrow_y = [arrow_h, arrow_h];
annotation('textarrow', arrow_x, arrow_y, 'String', 'time ','Color',foreColor,'FontSize',12);
% add dashes between the two temp region halves
arrow_x = [0.495, 0.55];
arrow_y = [0.11, 0.11];
annotation('line', arrow_x, arrow_y,'Color',foreColor,...
    'linestyle','--','linewidth', 1.8);


% Save the Figure
save_figure(fig, [fig_dir title_str ' tuning curve between regions']);


%% Instances where female is chasing male

% Are there any instances of female flies chasing male flies?
% If so, when / where does that occur

% Adapted from 5_1 analysis to instances of male chasing female

clearvars('-except',initial_var{:})

x = 1;
y = 2;

ctime = 1; % time restrictor in seconds for a single bout
fps = fly(1).fps;

fchase = [];

for trial = 1:num.trials
    fchase(trial).name = fly(trial).name;
    % positions of M head and F head and center
    P1 = [fly(trial).f.pos(:,body.head,x),fly(trial).f.pos(:,body.head,y)]; % female head
    P2 = [fly(trial).m.pos(:,body.center,x),fly(trial).m.pos(:,body.center,y)]; % male center
    P3 = [fly(trial).m.pos(:,body.head,x),fly(trial).m.pos(:,body.head,y)]; % male head
    
    % 1) Calculate body vectors
    v1 = P3 - P1;  % Nx2 matrix, vector for male head to female head
    v2 = P3 - P2;  % Nx2 matrix, vector for male head to male center
    
    % 2) Calculate the dot product of v1 and v2 for each time step
    dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);
    
    % 3) Compute the magnitudes of the vectors
    mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2); 
    mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2); 
    
    % 4) Calculate the cosine of the angle
    cosTheta = dotProduct ./ (mag_v1 .* mag_v2);
    
    % 5) Compute the angle in radians and convert to degrees
    angleRadians = acos(cosTheta);  % angle in radians
    angleDegrees = rad2deg(angleRadians);  % convert to degrees
    mfpos_angle = angleDegrees;
    
    % Identify when female position angle is less than 60 degrees from male
    pos_angle = abs(mfpos_angle) <= 60;
    
    % Identify when female is facing male
    facing = [];
    f_items = {'L3', 'L4', 'GX3', 'GX4', 'GY3', 'GY4'};
    % Identify if female is in an appropriate position for each wing direction across each item
    for i = 1:length(f_items)
        facing = [facing, fly(trial).position.(f_items{i})];
    end
    facing = any(facing,2);
    
    % Identify when female is behind male AND facing him
    fbehindm = (facing & pos_angle);
    
    % Identify when female is within 7mm of male
    close_dist = data.IFD(:,trial) <= 7; % mm
    
    % Identify when female is moving
    fmoving = data.speed(:,F,trial) >= 0.1; % min speed up for debate
    
    % Identify when male is moving
    mmoving = data.speed(:,M,trial) >= 0.01; % min speed diff than F in order to include true chase bouts
    
    chase = (fbehindm & close_dist & fmoving & mmoving); % requirements for 'chase all'
    a = diff(chase); 
    % Add the first chase value to the list to account for the starting condition
    b = [chase(1); a]; 
    % Locations in chase where chasing period starts/end
    ch_start = find(b == 1); 
    ch_stop = find(b == -1);
    % If chasing doesn't stop by end, add stop location at end of ch_stop (loc = length of experiment value)
    if chase(end)
        ch_stop(end + 1) = length(time);
    end
    % Calculate the length of each chasing bout
    ch_dur = ch_stop - ch_start;
    % Find where chasing lasts longer than 2sec
    dur_loc = find(ch_dur > (ctime*fps));
    
    % Create new courtship matrix with only true chasing bouts longer than 2sec
    mt = false(size(fly(trial).time));
    if isempty(dur_loc)
        fchase(trial).(['roi_' num2str(ctime) 'sec']) = [];
    else
        for i = 1:length(dur_loc)
            ii = dur_loc(i);
            mt(ch_start(ii):ch_stop(ii)) = true;
            fchase.(['roi_' num2str(ctime) 'sec'])(i,:) = [ch_start(ii), ch_stop(ii)];
        end
    end
    fchase(trial).(['sec_' num2str(ctime)]) = mt; % time restriction 2 seconds
    fchase(trial).roi_all = find(chase);
    fchase(trial).all = chase; % NO time limit
end

% Create new Female Chase folder for saving figs
chDir = 'S:\Evyn\DATA\Courtship Videos/grouped/Berlin LTS caviar/Figures/Female Chase Figures/';
if ~exist(chDir, 'dir')
        mkdir(chDir)
end

loc = [];
loc_2 = [];
y = [];
y2 = [];

% Determine y variables for time restriction and no time limit
for trial = 1:num.trials    
    loc = size(fchase(trial).(['roi_' num2str(ctime) 'sec']),1); % with time restriction
    loc_2 = size(fchase(trial).roi_all,1); % no time limit
    y = [y,loc];
    y2 = [y2,loc_2];
end

% X variables
x = 1;
x2 = 2;

% FIGURE
fig = getfig('', 1, [436 620]);
hold on
% Plot bar plot
bar(x,y,'FaceColor', Color('PaleVioletRed'),'EdgeAlpha',0) % time restricted
bar(x2,mean(y2),'FaceColor', Color('PaleVioletRed'),'EdgeAlpha',0) % all
scatter(x,y,50,Color('MediumVioletRed'),'filled')
scatter(x2,y2,50,Color('MediumVioletRed'),'filled')

% Format figure
formatFig(fig,blkbgd);
xticks([1 2])
xlim([0 2.75])
xticklabels({[num2str(ctime) ' sec'], 'all'})
ylabel('Number of frames of female chase')

% Save figure
save_figure(fig,[chDir 'Frame count of female chase instances ' num2str(ctime) ' sec'],fig_type);

%% Courtship attempts while female is sleeping

scourt = [];

% CI
fsleep = [];
mCI= [];

for trial = 1:num.trials
    scourt(trial).name = fly(trial).name;
    fsleep = data.sleep(:, F, trial);
    fsleep(isnan(fsleep)) = 0;

    % CI
    mCI = data.CI(:, trial);
    mCI(isnan(mCI)) = 0;
    loc = fsleep & mCI;
    idx = find(loc);
    if isempty(idx)
        scourt(trial).numframes_CI = 0;
    else
        scourt(trial).numframes_CI = length(idx);
    end

    % Wing ext
    mwing = data.wing_ext(:, trial);
    mwing(isnan(mwing)) = 0;
    loc = fsleep & mwing;
    idx = find(loc);
    if isempty(idx)
        scourt(trial).numframes_WE = 0;
    else
        scourt(trial).numframes_WE = length(idx);
    end

    % Chase
    mchase = data.court_chase(:, trial);
    mchase(isnan(mchase)) = 0;
    loc = fsleep & mchase;
    idx = find(loc);
    if isempty(idx)
        scourt(trial).numframes_chase = 0;
    else
        scourt(trial).numframes_chase= length(idx);
    end

    % Circling
    mcircling = data.circling_1sec(:, trial);
    mcircling(isnan(mcircling)) = 0;
    loc = fsleep & mcircling;
    idx = find(loc);
    if isempty(idx)
        scourt(trial).numframes_circling= 0;
    else
        scourt(trial).numframes_circling= length(idx);
    end
end

%% 

fps = fly(1).fps;
r = 5;
c = 1;
sb(1).idx = 1; % temperature
sb(2).idx = 2:5; % raster plot
lw = 2;
sSpan = 5*fps;
spike_H = 0.25; %height of each raster
trial_space = 1; %gap between trial lines
spike_W = 0.5; % raster line width

% Create y axis plot points for each trial
dummy = 1;
y = [];
for trial = 1:num.trials    
    ystart = dummy;
    yend = ystart + 2;
    a = ystart:yend;
    y = autoCat(y,a,true,true);
    dummy = ystart + 3;
end

% FIGURE
fig =  getfig('',1); 
    % Plot temperature
    subplot(r,c,sb(1).idx); hold on
        plot(data.time, data.temp, 'color', foreColor,'linewidth', lw)
        ylabel('temp (\circC)')
    % Plot raster of female sleep and male court instances for each trial
    subplot(r,c,sb(2).idx); hold on   
        for trial = 1:num.trials    
            % Pull female sleep and male CI 
            fsleep = data.sleep(:, F, trial);
            fsleep(isnan(fsleep)) = 0;
            mcourt = logical(data.CI(:, trial));
            mcourt(isnan(mcourt)) = 0;
            % X values
            fx = fly(trial).time(fsleep);
            fX = [fx';fx']; % duplicate and transpose to match size of y matrix for raster
            mx = fly(trial).time(mcourt);
            mX = [mx';mx'];
            % Y values
            fY = repmat([y(trial,F);y(trial,F)+spike_H],[1,size(fx,2)]);    
            mY = repmat([y(trial,M);y(trial,M)+spike_H],[1,size(mx,2)]); 
            plot(fX,fY,'color',Color(data,color(F,:)),'linewidth',spike_W)
            plot(mX,mY,'color',Color(data.color(M,:)),'linewidth',spike_W)
        end

       
        
