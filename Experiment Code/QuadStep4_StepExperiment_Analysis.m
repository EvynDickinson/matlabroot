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

%% Compare movement and speed



% Average speed for each temperature
for trial = 1:ntrials
    raw = data(trial).speed.avg;
    temp = data(trial).occupancy.temp;
    %ROI of useable data:
    tPoints = getTempTurnPoints(T.TempProtocol{trial});
    %sort into temperature bins
    
    
    
end

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



% group speed over temperature
for trial = 1:ntrials
    % TODO...
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
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    