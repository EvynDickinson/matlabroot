% This script looks at data files that had no food and used the step ramp
% temperature protocols

%% Initial ANALYSIS: General data organization (forward and backwards compatability
if ~exist('pix2mm','var')
    pix2mm = 12.8;
end
% fill out the data structure
for trial = 1:ntrials
    if isempty(data(trial).occupancy)
        data(trial).occupancy = data(trial).data.occupancy;
    end
    if isempty(data(trial).nflies)
        data(trial).nflies = data(trial).data.nflies;
    end
    if isempty(data(trial).wellLabels)
        data(trial).wellLabels = data(trial).data.wellLabels;
    end
    if isempty(data(trial).dist2wells)
        try
            data(trial).dist2wells = data(trial).data.occupancy.dist2wells;
        catch
            for well = 1:4
                data(trial).dist2wells(:,well) = data(trial).occupancy.dist2wells(well).N(:,1)./pix2mm;
            end
        end
    elseif isfield(data(trial).dist2wells, 'N')
        temp = data(trial).dist2wells;
        data(trial).dist2wells = [];
        for well = 1:4
            data(trial).dist2wells(:,well) = temp(well).N(:,1)./pix2mm;
        end
    end
end

%% ANALYSIS: Average speed for binned temperatures
sSpan = 180;
binSpace = 1;
binbuffer = 0.5;

[tempAll, plotData] = deal([]); 
% Average speed for each temperature
for trial = 1:ntrials
    % pull data for use:
    raw = data(trial).speed.avg;
    temp = (data(trial).occupancy.temp);
    % ROI of useable data:
    tPoints = getTempTurnPoints(T.TempProtocol{trial});
    ROI = [tPoints.DownROI';tPoints.UpROI';tPoints.HoldROI']; 
    tempdata = temp(ROI);
    speeddata = raw(ROI);

    % determine temperature bins
    tLow = tPoints.threshLow-binbuffer;
    tHigh = tPoints.threshHigh+binbuffer;    
    t_roi = tLow:binSpace:tHigh; 
%     t_roi(1) = tLow; t_roi(end) = tHigh;
    % sort data by temperature:
    binID = discretize(tempdata,t_roi);
    for tt = 1:length(t_roi)-1
        loc = binID==tt;
        plotData(trial,tt) = mean(speeddata(loc),1,'omitnan');
        tempAll(trial,tt) = mean(t_roi(tt:tt+1));
    end
end

speed.avg = plotData;
speed.temp = tempAll;

initial_vars{end+1} = 'speed';
clearvars('-except',initial_vars{:})

%% FIGURE: plot avg speed for each temperature

fig = figure;
hold on
for trial = 1:ntrials
    x = speed.temp(trial,:);
    y = speed.avg(trial,:);
    x(isnan(y)) = [];
    y(isnan(y)) = [];
    
    plot(x,y,'Marker','o','LineWidth',1,'color','w')
end
xlabel('Temperature (\circC)')
ylabel('Avg Speed (mm/s)')
formatFig(fig,true);

save_figure(fig, [figDir ExpGroup ' avg speed vs temperature'], '-png');

%% FIGURE: Compare speed for each step to it's control over time

% TODO...










figure; plot(binID)
    

% Pull plotting data:
for ii = 1:nfoods+1 %each food type
    % cut off the high and low ends of data (to clean):
    loc = food(ii).N(:,1)>threshHigh | food(ii).N(:,1)<threshLow;
    food(ii).N(loc,:) = [];
    % sort all the data by temperature:
    food(ii).N(:,3) = discretize(food(ii).N(:,1),t_roi);
    for tt = 1:length(t_roi)
        loc = food(ii).N(:,3)==tt;
        food(ii).avg(tt) = mean(food(ii).N(loc,2),1,'omitnan');
        food(ii).err(tt) = std(food(ii).N(loc,2),0,1,'omitnan')/sqrt(ntrials);
    end
end
    


row = 3;
col = 3;
sb(1).idx = [1,2];
sb(2).idx = [4,5,7,8];
sb(3).idx = [3,6,9];


fig = figure;
% temperature
subplot(row, col, sb(1).idx)
    x = data(1).occupancy.time;
    y = data(1).occupancy.time;
    plot(x,y,'color', 'w', 'LineWidth',1)
    ylabel('\circC')

% speed
subplot(row, col, sb(2).idx)
    hold on
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    