

%% Speed between regions timecourse

% pull out when flies in outer ring = 1
% pull out speed

clearvars('-except',initial_var{:})
% Create lists to cycle through
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color

% FIGURE
fig = getfig;
hold on
% Pull out fly speed in each region
for ii = 1:length(region_name)
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
    plot(data.time,smooth(m,300,'moving'),"Color",Color(color_list{ii}))
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

% Create lists to cycle through
regime_name = {'WT','WS','CT','CS','SS'}; % temp regime
region_name = {'OutterRing','innerFoodQuad','innerEmptyQuad'}; % region
color_list = {'MetroPurple','MetroRed','MetroOrange'}; % color


labels = {'Warm Threat', 'Warm Safe', 'Cool Threat', 'Cool Safe'};
sz = 40;
sex_buff = 0.2;
[xM, xF] = deal(1:4);
xM = xM-sex_buff;
xF = xF+sex_buff;

h = [];
x = NaN(3,4);

% FIGURE
fig = getfig;
hold on
% Pull out fly speed in each regime
for ii = 1:length(regime_name)
    current_var = data.tempbin.(regime_name{ii});
    loc = logical(replaceNaN(current_var,0));    
    speed = data.speed;
    % Replace speed values not during the wanted temp regime with nans
    speed(~loc) = nan;
    % Remove M/F fly dimension - all flies and their speeds together
    speed_all = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
    h = [h,ii];
    % Pull out fly speed in each temp region  
    for region = 1:length(region_name)
        current_reg = data.(region_name{region});
        loc2 = replaceNaN(current_reg,0);
        speed2 = speed_all;
        % Replace speed values not during the wanted temp regime with nans
        speed2(~loc2) = nan;
        % Calculate mean
        mspeed = mean(speed2,2,'omitnan');   
        h = [h,h(end)+4];
        x(ii,:) = h;
        % Plot average speed in each region during each temp regime
        
    end
end






% % FIGURE
% fig = getfig;
% hold on
% % Pull out fly speed in each region
% for ii = 1:length(region_name)
%     current_var = data.(region_name{ii});
%     loc = logical(replaceNaN(current_var,0));    
%     speed = data.speed;
%     % Replace speed values when fly is not in the wanted region with nans
%     speed(~loc) = nan;
%     % Remove M/F fly dimension - all flies and their speeds together
%     speed_all = [squeeze(speed(:,M,:)),squeeze(speed(:,F,:))];
%     % Pull out fly speed in each temp regime
%     for regime = 1:length(regime_name)
%         current_reg = data.tempbin.(regime_name{regime});
%         loc2 = replaceNaN(current_reg,0);
%         speed2 = speed_all;
%         % Replace speed values not during the wanted temp regime with nans
%         speed2(~loc2) = nan;
%         % Calculate mean
%         mspeed = mean(speed2,2,'omitnan');   
%         x = 
%         % Plot average speed in each region during each temp regime
% 
%     end
% end

