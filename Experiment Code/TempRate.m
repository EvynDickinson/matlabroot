

% Run this script after initial data processing with:
% QuadStep1.m
% QuadStep2.m
% GroupDataGUI
clear 

%% Load data file
% Load data for now from the first step of QuadStep3.m

clearvars('-except',initial_vars{:})

%% Trial by trial figures: 
disp_figs = true;
auto_save_figs = true;
ind_fig_loc = [figDir 'Trial by trial\'];
if ~isfolder(ind_fig_loc); mkdir(ind_fig_loc); end
vars = [initial_vars(:)', 'trial', 'threshHigh', 'threshLow', 'binSpace',...
        'G', 'disp_figs','auto_save_figs','ind_fig_loc'];

[threshHigh, threshLow] = getTempThresholds;
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 

G = struct;

% Get the temp rate, temp, and distance from food
for trial = 12:ntrials
    T_rates = [];
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

    % === Demo plot for single trial =====:
    title_str = [T.Date{trial} ' Arena ' T.Arena{trial}];
    if disp_figs 
        fig_file = [ind_fig_loc ExpGroup ' temp temp_rate dist ' title_str];
        if isfile([fig_file '.png']) %skip already generated figures
            continue
        end
        time = data(trial).occupancy.time(1:end-1);
        fig = figure; set(fig, 'pos', [689 48 1628 960])
        % temp
        subplot(3,1,1)
        plot(time, plotData(:,1),'color','w')
        v_line(time(tPoints.transitions),'y',':')
        ylabel('Temp (\circC)')
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
        save_figure(fig, fig_file, '-png',auto_save_figs);
    end
    % Temp-rate identification and sorting: 
    buffSize = 0.05;
    for ii = 1:tPoints.nRates
        edges(ii,:) = [tPoints.rates(ii)-buffSize, tPoints.rates(ii)+buffSize];
    end  
    rateData = plotData(:,4);
    nRates = tPoints.nRates;
    rateIdx = discretize(rateData,nRates);
    for ii = 1:nRates
        TD = rateData(rateIdx==ii);
        T_rates(ii) = round(mean(TD),2);
        % do these match the assumed rates?
        idx = find(T_rates(ii) > edges(:,1) & T_rates(ii) < edges(:,2));
        if isempty(idx)
            warndlg('Temp rate not within expected range')
            return
        end
        % set the rate value to the uniform cross-group value:
        rateData(rateIdx==ii) = tPoints.rates(ii);
        T_rates(ii) = tPoints.rates(ii);
%        figure; hold on
%        plot(rateData)
%        plot(plotData(:,4))

    end
    plotData(:,4) = rateData;
    plotData(:,5) = rateIdx; %rate bin
    G(trial).rateIdx = rateIdx;
    G(trial).rates = T_rates;
    G(trial).data = plotData;
    
    % Temperature range and bin size formatting
    t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
    if t_roi(end)<ceil(threshHigh)
        t_roi(end+1) = ceil(threshHigh) + binSpace;
    end
    nTemps = length(t_roi);
    tempIdx = discretize(plotData(:,1), t_roi);
    G(trial).nTemps = nTemps;
    G(trial).tempIdx = tempIdx;
    G(trial).temps = t_roi;

    % turn coordinates of heatmap into a vector to fill in 2D
    heatMapData = [];
    for col = 1:nTemps
        for row = 1:nRates
        % Pull the data that matches this category
        loc = (rateIdx==row) & (tempIdx==col);
        heatMapData(row,col) = mean(plotData(loc,3));
        end
    end
    G(trial).heatmap = heatMapData;

    % PLOT DISTANCE AS FUNCTION OF TEMP AND RATE OF TEMP
    if disp_figs
        fig_file = [ind_fig_loc ExpGroup ' temp temp_rate dist heatmap ' title_str];
        if isfile([fig_file '.png']) %skip already generated figures
            continue
        end
%         colormap default
        fig = figure; set(fig, 'pos', [560 127 983 417]);
        hold on
        imAlpha=ones(size(heatMapData));
        imAlpha(isnan(heatMapData))=0;
        imagesc(heatMapData,'AlphaData',imAlpha);
        set(gca,'color',0*[1 1 1]);
        axis tight
        title(title_str)
        % Axes formatting
        ax = gca;
        fig = formatFig(fig, true);
        XtickNum = ax.XTick;
        ax.XTickLabel = t_roi(XtickNum);
        YtickNum = ax.YTick;
        set(gca, 'ytick', 1:nRates,'YTickLabel',T_rates)
        ylabel('\DeltaT/dt (\circC/min)')
        xlabel('Temp (\circC)')
        % Colorbar formatting
        cbh = colorbar(); 
        cbh.Label.String = 'Distance from food (mm)';
        cbh.Color = Color('white');
        % flip colormap around to make yellow closer to food
        cmap = colormap;
        set(gca, 'colormap', flip(cmap))
        save_figure(fig, fig_file, '-png', auto_save_figs);
    end
    
clearvars('-except',vars{:}) 
end

% Save Data structure:
if questdlg('Save loaded data?')
    save([figDir ExpGroup ' temp rate'])
end


%% Grouped hysteresis heatmap output 
% TODO: reformat this to take advantage of the previously processed data
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

clearvars('-except',initial_vars{:})



%% How do the temp-distance relationships compare between heating and cooling for each temp rate?

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
All_tempRates = [];
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
    All_tempRates = [All_tempRates;T_rates]; %save the rates into group structure
    G(trial).t_rates = T_rates;
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


% Find temperature rates that align (increasing and decreasing)
groupRates = unique(abs(All_tempRates));
nRates = length(groupRates);

% 
% for tt = 1:nRates
%     speed = groupRates(tt);
%     % find opposing temperature rate:
%     loc = find(groupRates==speed);
    
    


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

clearvars('-except',initial_vars{:})





% Organize data for plotting - bin by temp (deg)
emptyData = false(1,3);
plotData = [];
t_roi = floor(threshLow):ceil(threshHigh); 
for K = 1:3 %food type
    for type = 1:2
        % pull appropriate data:
        switch K
            case 1 %plant
                UpData = plant.up;
                DownData = plant.down;
            case 2 %yeast
                UpData = yeast.up;
                DownData = yeast.down;
            case 3 %empty
                UpData = empty.up;
                DownData = empty.down;
        end
        if isempty(UpData) || isempty(DownData)
            emptyData(K) = true;
            continue 
        end
        if type == 1 
            inputData = UpData;
        else 
            inputData = DownData;
        end
%         [loc,idx,cnt_unique,unique_a,len,mt] = deal([]);
        % cut off the high and low ends of data (to clean):
        loc = inputData(:,1)>threshHigh | inputData(:,1)<threshLow;
        inputData(loc,:) = [];
        
        % sort all the data by temperature:
        idx = discretize(inputData(:,1),t_roi);
        [cnt_unique, unique_a] = hist(idx,unique(idx));
        len = max(cnt_unique);
        mt = nan(len,length(unique_a));
        for tt = 1:length(unique_a)
            cue = unique_a(tt); %index number
            loc = idx==cue;
            mt(1:sum(loc),tt) = inputData(loc,2);
            y_err(tt) = std(inputData(loc,2));
            plotData(K,type).y_avg(tt) = mean(inputData(loc,2));
        end
        plotData(K,type).y_err = y_err./sqrt(ntrials);
        plotData(K,type).xdata = t_roi(unique_a);
        inputData = [];
    end
end
 

% LINE GRAPH PLOT CODE
% PLOT the grouped & binned data points :
nPlots = sum(~emptyData);
nrows = 1; ncols = nPlots;
titleList = {'Plant', 'Yeast', 'Empty'};
CList = {'red', 'deepskyblue'}; %heating and cooling colors

ii = 0;
fig = figure; set(fig, 'pos', [132 83 365*nPlots 693]);
for K = 1:3
    % skips absent food types
    if emptyData(K)
        continue
    else
        ii = ii+1;
    end
    subplot(nrows, ncols, ii)
    hold on
    for type = 1:2
        kolor = Color(CList{type});
        x = plotData(K,type).xdata;
        y = plotData(K,type).y_avg;
        yerr = plotData(K,type).y_err;
        fill_data = error_fill(x, y, yerr);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
          set(h, 'facealpha', 0.2)
        plot(x,y,'color', kolor, 'linewidth', 2)
    end
    ylimits(ii,:) = ylim;
    xlabel('temperature (\circC)')
    ylabel('distance from well (mm)')
    title(titleList{K})
end
fig = formatFig(fig,true, [nrows,ncols]);
%set uniform y axis
for ii = 1:nPlots
subplot(nrows,ncols,ii)
ylim([min(ylimits(:,1)),max(ylimits(:,2))])
end
l = legend({'SEM','Heating', 'SEM','Cooling'});
set(l, 'textcolor', 'w','position', [0.5997 0.1556 0.1781 0.1176])
save_figure(fig, [figDir ExpGroup ' temp hysteresis'], '-png');






























