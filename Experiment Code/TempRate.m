

% Run this script after initial data processing with:
% QuadStep1.m
% QuadStep2.m
% GroupDataGUI

%% Load data file
% Load data for now from the first step of QuadStep3.m

clearvars('-except',initial_vars{:})

%% Trial by trial figures: 

[threshHigh, threshLow] = getTempThresholds;

% Get the temp rate, temp, amd distance from food
trial = 4;
% temperature
temp = data(trial).occupancy.temp;
plotData(:,1) = temp(1:end-1); 
% distance from food
for well = 1:4
    [kolor,num] = pullFoodColor(data(trial).wellLabels{well});
    if num==1
        y = data(trial).occupancy.dist2wells(well).N(:,1)./pix2mm;
        plotData(:,3) = y(1:end-1);
        break
    end
end

% rate of temperature change
dT = diff(smooth(temp,180)); %180 = 1 minute smoothing kernal
dt = diff(data(trial).occupancy.time); 
plotData(:,2) = dT./dt'; 

% find the mean temp rate during each ramp period:
tPoints = getTempTurnPoints(T.TempProtocol{trial}); %accomodates multiple temp protocols within the data group
for ii = 1:length(tPoints.transitions)-1
    roi = tPoints.transitions(ii):tPoints.transitions(ii+1);
    distD = median(plotData(roi,2));
    plotData(roi,4) = round(distD*ones(range(roi)+1,1),2);
end

% Demo plot for single trial:
time = data(trial).occupancy.time(1:end-1);
fig = figure; set(fig, 'pos', [689 48 1628 960])
% temp
subplot(3,1,1)
plot(time, plotData(:,1),'color','w')
v_line(time(tPoints.transitions),'y',':')
ylabel('Temp (\circC)')
title_str = [T.Date{trial} ' Arena ' T.Arena{trial}];
title(title_str)
% temp rate
subplot(3,1,2); hold on
plot(time, plotData(:,2),'color','w')
plot(time, plotData(:,4),'color','r','linewidth', 2)
v_line(time(tPoints.transitions),'y',':')
ylabel('dT/dt (\circC/min)')
% distance from food
subplot(3,1,3)
plot(time, plotData(:,3),'color',kolor)
v_line(time(tPoints.transitions),'y',':')
xlabel('time (min)')
ylabel('distance (mm)')
yyaxis right
plot(time, smooth(data(trial).occupancy.occ(1:end-1,well),180),'color', 'w')
ylabel('occupancy (prob)')
fig = formatFig(fig, true, [3,1]);
subplot(3,1,3)
yyaxis left 
set(gca,'YColor', 'w')
save_figure(fig, [figDir ExpGroup ' temp temp_rate dist ' title_str], '-png');


% single trial heat map test: 

% sort data into ramp regions: 

% Temp-rate identification and sorting: 
T_rates = sort(unique(plotData(:,4)));
buff = diff(T_rates)/2;
buff = [buff(1); buff];
bin_edges = T_rates - buff;
bin_edges = [bin_edges; T_rates(end)+buff(end)];
rateIdx = discretize(plotData(:,4), bin_edges);
nRates = length(T_rates);


% Temperature range and bin size formatting
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 
t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
if t_roi(end)<ceil(threshHigh)
    t_roi(end+1) = ceil(threshHigh) + binSpace;
end
nTemps = length(t_roi);
tempIdx = discretize(plotData(:,1), t_roi);

% assign distance data to the appropriate place on the heatmap

% turn coordinates of heatmap into a vector to fill in 2D
heatMapData = [];
for col = 1:nTemps
    for row = 1:nRates
    % Pull the data that matches this category
    loc = (rateIdx==row) & (tempIdx==col);
    heatMapData(row,col) = mean(plotData(loc,3));
    end
end

% PLOT DISTANCE AS FUNCTION OF TEMP AND RATE OF TEMP
colormap default

fig = figure; set(fig, 'pos', [560 127 983 417]);
hold on
imAlpha=ones(size(heatMapData));
imAlpha(isnan(heatMapData))=0;
h = imagesc(heatMapData,'AlphaData',imAlpha);
set(gca,'color',0*[1 1 1]);
axis tight
title(title_str)
% Axes formatting
ax = gca;
fig = formatFig(fig, true);
XtickNum = ax.XTick;
ax.XTickLabel = t_roi(XtickNum);
YtickNum = ax.YTick;
ax.YTickLabel = T_rates(YtickNum);
ylabel('\DeltaT/dt (\circC/min)')
xlabel('Temp (\circC)')
% Colorbar formatting
cbh = colorbar(); 
cbh.Label.String = 'Distance from food (mm)';
cbh.Color = Color('white');
% flip colormap around to make yellow closer to food
cmap = colormap;
set(gca, 'colormap', flip(cmap))

save_figure(fig, [figDir ExpGroup ' temp temp_rate dist heatmap ' title_str], '-png');


clearvars('-except',initial_vars{:})



%% Grouped hysteresis heatmap output
% Fly by fly average:

[threshHigh, threshLow] = getTempThresholds;
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 
foodNum = 1; % search for plant food currently
% identify the food well:
for trial = 1:ntrials
    % distance from food
    for well = 1:4
        [kolor,num] = pullFoodColor(data(trial).wellLabels{well});
        if num==foodNum
            wellIdx(trial) = well;
            break
        end
    end
end


% Get the temp rate, temp, amd distance from food
for trial = 1:ntrials

    % temperature
    temp = data(trial).occupancy.temp;
    G(trial).plotData(:,1) = temp(1:end-1); 
    % distance from food
    well = wellIdx(trial);
    y = data(trial).occupancy.dist2wells(well).N(:,1)./pix2mm;
    G(trial).plotData(:,3) = y(1:end-1);


    % rate of temperature change
    dT = diff(smooth(temp,180)); %180 = 1 minute smoothing kernal
    dt = diff(data(trial).occupancy.time); 
    G(trial).plotData(:,2) = dT./dt'; 

    % find the mean temp rate during each ramp period:
    tPoints = getTempTurnPoints(T.TempProtocol{trial}); %accomodates multiple temp protocols within the data group
    for ii = 1:length(tPoints.transitions)-1
        roi = tPoints.transitions(ii):tPoints.transitions(ii+1);
        distD = median(G(trial).plotData(roi,2));
       G(trial).plotData(roi,4) = round(distD*ones(range(roi)+1,1),2);
    end

    % Temp-rate identification and sorting: 
    T_rates = sort(unique(G(trial).plotData(:,4)));
    buff = diff(T_rates)/2;
    buff = [buff(1); buff];
    bin_edges = T_rates - buff;
    bin_edges = [bin_edges; T_rates(end)+buff(end)];
    rateIdx = discretize(G(trial).plotData(:,4), bin_edges);
    G(trial).rateIdx = rateIdx;
    nRates = length(T_rates);


    % Temperature range and bin size formatting
    t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
    if t_roi(end)<ceil(threshHigh)
        t_roi(end+1) = ceil(threshHigh) + binSpace;
    end
    nTemps = length(t_roi);
    tempIdx = discretize(G(trial).plotData(:,1), t_roi);
    G(trial).tempIdx = tempIdx;

    % assign distance data to the appropriate place on the heatmap
    heatMapData = [];
    for col = 1:nTemps
        for row = 1:nRates
        % Pull the data that matches this category
        loc = (rateIdx==row) & (tempIdx==col);
        heatMapData(row,col) = mean(G(trial).plotData(loc,3));
        end
    end
    G(trial).hMap = heatMapData;
end

% Group the distance data (average of trials -- not combined all mashed up)
heatMapData = [];
Hmap = nan(nRates,nTemps,ntrials);
for trial = 1:ntrials
    Hmap(:,:,trial) = G(trial).hMap;
end
heatMapData = mean(Hmap,3);
loc = isnan(heatMapData);

% PLOT DISTANCE AS FUNCTION OF TEMP AND RATE OF TEMP
title_str = [ExpGroup];
colormap default

fig = figure; set(fig, 'pos', [560 127 983 417]);
hold on
imAlpha=ones(size(heatMapData));
imAlpha(isnan(heatMapData))=0;
h = imagesc(heatMapData,'AlphaData',imAlpha);
set(gca,'color',0*[1 1 1]);
axis tight
title(title_str)
% Axes formatting
ax = gca;
fig = formatFig(fig, true);
XtickNum = ax.XTick;
ax.XTickLabel = t_roi(XtickNum);
YtickNum = ax.YTick;
ax.YTickLabel = T_rates(YtickNum);
ylabel('\DeltaT/dt (\circC/min)')
xlabel('Temp (\circC)')
% Colorbar formatting
cbh = colorbar(); 
cbh.Label.String = 'Distance from food (mm)';
cbh.Color = Color('white');
% flip colormap around to make yellow closer to food
cmap = colormap;
set(gca, 'colormap', flip(cmap))

save_figure(fig, [figDir ExpGroup ' temp temp_rate dist heatmap ' title_str], '-png');

% For quick look:
% TODO: more here!
% heatMapData...


%% 



































