

%% Speed between regions timecourse

% pull out when flies in outer ring = 1
% pull out speed

clearvars('-except',initial_var{:})
% Create lists to cycle through
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color
nRegions = length(region_name);

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
    % Calculate mean
    m = mean(speed_all,2,'omitnan');    
    % Plot average speed in each region over time
    plot(data.time,smooth(m,300,'moving'),'Color',Color(color_list{ii}))
end

% Format fig
formatFig(fig,blkbgd);
legend(region_name)
ylabel('average speed (mm/s)')
xlabel('time (min)')

% Save figure
% save_figure(fig,[fig_dir 'speed between regions timecourse'],fig_type);)

%% Speed between regions within temp regimes scatter plot

% pull out when flies are in each region
% pull out speed for each temp regime time period

clearvars('-except',initial_var{:})
% Create lists to cycle through
regime_name = {'WT','WS','CT','CS'}; % temp regime
labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color

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


sz = 40;
buff = 0.2;
% [xM, xF] = deal(1:4);
% xM = xM-buff;
% xF = xF+buff;

% FIGURE
fig = getfig;
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
        % Replace speed values not during the wanted temp regime with nans
        speed2(~loc2,:) = nan;
        % Calculate mean
        avgspeed = mean(speed2,1,'omitnan');   

        % Plot average speed in each region during each temp regime
        scatter(x(region,regime),avgspeed,sz,Color(color_list{region}),'filled');
    end
end

formatFig(fig,blkbgd);
legend(region_name)
ylabel('average speed per fly (mm/s)')


%% Speed between temp regimes within regions scatter plot

% pull out when flies are in each region
% pull out speed for each temp regime time period

clearvars('-except',initial_var{:})
% Create lists to cycle through
regime_name = {'WT','WS','CT','CS'}; % temp regime
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color

% Set up x axis placement
h = [];
x = [];
for regime = 1:nRegimes % 3
    h = regime;
    for region = 1:nRegions % 4
        h = [h,h(end)+5];
    end
    x(regime,:) = h(:,1:nRegions);
end

labels = {'Outer Ring', 'Inner Food Quad', 'Inner Empty Quad'};
sz = 40;
buff = 0.2;
% [xM, xF] = deal(1:4);
% xM = xM-buff;
% xF = xF+buff;

Rspeed = [];
% FIGURE
fig = getfig;
hold on
% Pull out fly speed in each region
for region = 1:nRegions % 3
    current_var = data.(region_name{region});
    loc = logical(replaceNaN(current_var,0));    
    speed = data.speed;
    % Replace speed values when fly is not in the wanted region with nans
    speed(~loc) = nan;
    % Pull out fly speed in each temp regime
    for regime = 1:nRegimes % 4
        current_reg = data.tempbin.(regime_name{regime});
        loc2 = replaceNaN(current_reg,0);
        speed2 = speed;
        % Replace speed values not during the wanted temp regime with nans
        speed2(~loc2) = nan;
        % Remove M/F fly dimension - all flies and their speeds together
        speed2 = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
        % Calculate mean
        avgspeed = mean(speed2,1,'omitnan');   
        % Plot average speed in each region during each temp regime
        scatter(x(regime,region),avgspeed,sz,Color(color_list{region}),'filled')
        Rspeed = [Rspeed;avgspeed];
    end
    for fly = 1:length(avgspeed)
        plot(x(:,region),Rspeed(:,fly))
    end
        Rspeed = []; 
end


formatFig(fig,blkbgd);
legend(region_name)
ylabel('average speed per fly (mm/s)')


%% Histogram of speed between regions

clearvars('-except',initial_var{:})
% Create lists to cycle through
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color

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
    histogram(speed_all)
end










