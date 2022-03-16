
%    GroupDataGUI

%% LOAD: multiple trials that are grouped into a structure
clear; warning off

if strcmpi('Yes', questdlg('Use excel named structure?','','Yes','No', 'Cancel', 'No'))
    % baseFolder = getCloudPath;
    [excelfile, Excel, xlFile] = load_QuadBowlExperiments;
    baseFolder = getCloudPath;
    
    % Select structure to load:
    [~,~,structInfo] = getExcelStructureNames(true);
    ExpGroup = structInfo.StructName;
    ntrials = structInfo.numTrials;
    
    % Make a data folder for structure images
    figDir = [baseFolder 'Data structures\' ExpGroup '\'];
    if ~exist(figDir, 'dir'); mkdir(figDir); end
    
    % Load data from each trial in the structure
    data = [];
    fprintf('\nLoading trials: \n')
    for trial = 1:ntrials
        % print the experiments as they are loaded
        trialExpID = excelfile{structInfo.rowNum(trial), Excel.expID};
        trialDate = excelfile{structInfo.rowNum(trial), Excel.date};
        trialArena = excelfile{structInfo.rowNum(trial), Excel.arena};
        trialName = [trialExpID ' Arena ' trialArena ' ' trialDate];
        disp(trialName)
        
        % build the path for the trial data
        dirc = [baseFolder, trialDate, '\Arena ' trialArena '\analysis\' trialExpID trialArena  ' timecourse data.mat'];
        
        % load data
        todel = load(dirc);
        variList = fieldnames(todel);
    %     todel = load(dirc, varList{:});
        data(trial).trialName = trialName;
        for ii = 1:length(variList)
            data(trial).(variList{ii}) = todel.(variList{ii});
        end
    end

else %select data structure from folder names out of GroupDataGUI:
    [baseFolder, folder] = getCloudPath(3);
    list_dirs = dir(folder); 
    list_dirs = {list_dirs(:).name};
    list_dirs(1:2) = [];
    idx = listdlg('ListString', list_dirs,'ListSize', [250, 400]);
    
    figDir = [folder list_dirs{idx} '/'];
    ExpGroup = list_dirs{idx};
    
    % load the directory list:
    load([figDir 'fileList.mat'])
    
    % extract information
    ntrials = size(T,1);
    dates = T.Date;
    arenas = T.Arena;
    expID = T.ExperimentID;
    
    % load data
    n = cell(1,ntrials);
    data = struct('dist2wells', n, 'wellLabels', n, 'occupancy', n,'nflies',n,'data', n);  
    var2load = fieldnames(data);
    
    for trial = 1:ntrials
        filePath = [baseFolder, dates{trial}, '/Arena ' arenas{trial} '/analysis/'];
        if ~isfolder(filePath)
            filePath = [baseFolder, dates{trial}, '/Arena ' arenas{trial} '/'];
        end
        
        temp = load([filePath expID{trial} arenas{trial} ' timecourse data.mat'], var2load{:});
        for ii = 1:length(var2load)
            try data(trial).(var2load{ii}) = temp.(var2load{ii});
            catch; data(trial).(var2load{ii}) = [];
            end
        end
        disp([expID{trial} arenas{trial}])
    end
end

%Pull and reorganize data within the structure:

pix2mm = 12.8; %conversion from pixels to mm for these videos
    
initial_vars = {'ExpGroup','baseFolder', 'T', 'data', 'figDir', 'filePath',...
                'initial_vars', 'folder', 'ntrials', 'pix2mm'};
clearvars('-except',initial_vars{:})
if questdlg('Save loaded data?')
    save([figDir ExpGroup ' raw'])
end
fprintf('Data loaded\n')

%% General data organization (forward and backwards compatability

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
                data(trial).dist2wells = data(trial).occupancy.dist2wells(well).N(:,1)./pix2mm;
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

% find the types of food available:
names = [];
for trial = 1:ntrials
    names = unique([data(trial).wellLabels; names]);
end
nfoods = length(names)-1;
foodNames = names;
foodNames(strcmp(names, 'Empty')) = [];

initial_vars{end+1} = 'nfoods';
initial_vars{end+1} = 'foodNames';
clearvars('-except',initial_vars{:})

%% FIGURE: visual check of temperature alignment across the experiments:
fig = figure; hold on
for trial = 1:ntrials
    X = data(trial).occupancy.time;
    Y = data(trial).occupancy.temp;
    plot(X, Y, 'linewidth', 1)
end
xlabel('Time (min)')
ylabel('Temp (\circ)')
title({'temperature alignment across experiments';...
      ['N = ' num2str(ntrials)]})
formatFig(fig, true);

save_figure(fig, [figDir ExpGroup ' temperature alignment'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% FIGURE: Distance to food vs temperature -- ALL trial lines
sSpan = 360;
fig = figure; hold on
for trial = 1:ntrials
    X = data(trial).occupancy.temp;
    for well = 1:4
        Y = data(trial).dist2wells(:,well);
        % reorder the data by temperature:
        [plotX, idx] = sort(X);
        plot(plotX, smoothdata(Y(idx),'movmean', sSpan), 'linewidth', 1,'color', pullFoodColor(data(trial).wellLabels{well}))
    end
end
ylabel('Distance to well (mm)')
xlabel('Temp (\circC)')
title({'Location from food sources by temperature';...
      ['N = ' num2str(ntrials)]})
formatFig(fig, true);
save_figure(fig, [figDir ExpGroup ' all temp vs distance'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% FIGURE: Average distance | occupancy vs temp across trials:

inputVar =  questdlg('Which data type to compare?','','distance','occupation probability','Cancel','distance');

food = struct;
switch inputVar
    case 'distance'
        ylab = 'Distance from well (mm)';
        L_loc = 'southwest';
    case 'occupation probability'
        ylab = inputVar;
        L_loc = 'northwest';
    case 'Cancel'
        return
end

%allocate empty data structure & set params
for ii = 1:nfoods+1
    food(ii).N = [];
    if ii <= nfoods
        food(ii).name = foodNames{ii};
    else 
        food(ii).name = 'Empty';
    end
end
for trial = 1:ntrials
    x = data(trial).occupancy.temp;
    for well = 1:4
        switch inputVar
            case 'distance'
                y = data(trial).dist2wells(:,well);
            case 'occupation probability'
                y = data(trial).occupancy.occ(:,well);
        end
        % empty food well:
        wellID = data(trial).wellLabels{well};
        if contains(wellID,'Empty')
            food(nfoods+1).N = [food(nfoods+1).N; x, y];
        else 
            loc = find(strcmp(wellID,foodNames));
            food(loc).N = [food(loc).N; x, y];
        end
    end
end 

% Temperature range and bin size formatting
[threshHigh, threshLow] = getTempThresholds;
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 
t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
if t_roi(end)<ceil(threshHigh)
    t_roi(end+1) = ceil(threshHigh) + binSpace;
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
    
% FIGURE: plot the avg. temp vs. distance data
CList = Color('SteelBlue', 'Navy', nfoods);
CList(nfoods+1,:) = Color('White');

fig = figure; set(fig,'pos', [67 82 675 692]);
hold on
for ii = 1:nfoods+1
    kolor = CList(ii,:);
    x = t_roi;
    y = food(ii).avg;
%     y_err = food(ii).err;
%     fill_data = error_fill(x,y, y_err);
%     h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
%     set(h, 'facealpha', 0.35)
    plot(x, y, 'Color', kolor, 'LineWidth', 2)
%     plot(x, y+y_err, 'Color', kolor, 'LineWidth', 0.25)
%     plot(x, y-y_err, 'Color', kolor, 'LineWidth', 0.25)
end
%Labels:
% xlim([7,20])
xlabel('temperature (\circC)')
ylabel(ylab)
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);
ax = gca;
set(ax, 'FontSize', 18)
str = strrep([foodNames; 'Empty'],'_',' ');
legend(str,'textcolor', 'w', 'location', L_loc, 'box', 'off','fontsize',12)

save_figure(fig, [figDir ExpGroup ' temp vs ' inputVar ' bin size ' num2str(binSpace)], '-png');
clearvars('-except',initial_vars{:})

%% FIGURE: Average movement | cluserting | eccentricity vs temp across trials:

inputVar =  questdlg('Which data type to compare?','','movement','clustering','eccentricity','movement');
switch inputVar
    case 'movement'
        ylab = 'movement (au)';
        L_loc = 'southeast';
    case 'clustering'
        ylab = 'Inter-fly-distance (mm)';
        L_loc = 'northwest';
    case 'eccentricity'
        ylab = 'eccentricity (mm)';
        L_loc = 'northwest';
    case ''
        return
end

%allocate empty data structure & set params
food = struct;
food.N = [];
for trial = 1:ntrials
    switch inputVar
        case 'movement'
            x = data(trial).occupancy.temp(1:end-1);
            y = data(trial).occupancy.movement;
        case 'clustering'
            y = data(trial).occupancy.IFD';
            x = data(trial).occupancy.temp;
        case 'eccentricity'
            y = data(trial).occupancy.eccentricity;
            x = data(trial).occupancy.temp;
    end 
    food.N = [food.N; x,y];
end 

% Temperature range and bin size formatting
[threshHigh, threshLow] = getTempThresholds;
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 
t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
if t_roi(end)<ceil(threshHigh)
    t_roi(end+1) = ceil(threshHigh) + binSpace;
end
    
% cut off the high and low ends of data (to clean):
loc = food.N(:,1)>threshHigh | food.N(:,1)<threshLow;
food.N(loc,:) = [];
% sort all the data by temperature:
food.N(:,3) = discretize(food.N(:,1),t_roi);
for tt = 1:length(t_roi)
    loc = food.N(:,3)==tt;
    food.avg(tt) = mean(food.N(loc,2),1,'omitnan');
    food.err(tt) = std(food.N(loc,2),0,1,'omitnan')/sqrt(ntrials);
end

    
% FIGURE: plot the avg. temp vs. distance data
kolor = Color('Teal');

fig = figure; set(fig,'pos', [67 82 675 692]);
hold on
x = t_roi;
y = food.avg;
y_err = food.err;
plot(x, y, 'Color', kolor, 'LineWidth', 2)
plot(x, y+y_err, 'Color', kolor, 'LineWidth', 0.25)
plot(x, y-y_err, 'Color', kolor, 'LineWidth', 0.25)

%Labels:
% xlim([7,20])
xlabel('temperature (\circC)')
ylabel(ylab)
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);
ax = gca;
set(ax, 'FontSize', 18)

save_figure(fig, [figDir ExpGroup ' temp vs ' inputVar ' bin size ' num2str(binSpace)], '-png');
clearvars('-except',initial_vars{:})








%% ANALYSIS: Get temp rate information from all trials
G = struct;
[threshHigh, threshLow] = getTempThresholds;
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 
vars = [initial_vars(:)', 'trial', 'threshHigh', 'threshLow', 'binSpace','G', 'vars'];

% Get the temp rate, temp, and distance from food
for trial = 1:ntrials
%     if isfield(data(trial).data, 'TR')
%         G(trial).TR = data(trial).data.TR;
%         continue
%     end
    [TR,T_rates,plotData] = deal([]); % temp rate data will go here then get saved into Arena Data

    % temperature
    temp = data(trial).occupancy.temp; temp = reshape(temp,numel(temp),1);
    time = data(trial).occupancy.time; time = reshape(time,numel(time),1);
    plotData(:,1) = temp(1:end-1);
    dT = diff(smooth(temp,180)); %180 = 1 minute smoothing kernal
    dt = diff(time); 
    plotData(:,2) = dT./dt; 

    % distance from food
    wellLabels = data(trial).wellLabels;
    for well = 1:4
        [kolor,num] = pullFoodColor(wellLabels{well});
        if ~strcmp(wellLabels{well},'Empty')
            y = data(trial).dist2wells(:,well);
            plotData(:,3) = y(1:end-1);
            foodWell = well;
            break
        end
    end

    % find the mean temp rate during each ramp period:
    tPoints = getTempTurnPoints(T.TempProtocol{trial}); %accomodates multiple temp protocols within the data group
    keepLoc = false(1,size(plotData,1));
    for ii = 1:length(tPoints.transitions)-1
        roi = tPoints.transitions(ii):tPoints.transitions(ii+1);
        distD = median(plotData(roi,2));
        plotData(roi,4) = round(distD*ones(range(roi)+1,1),2);
        keepLoc(roi) = true;
    end
    plotData(~keepLoc,:) = nan; % exclude data outside prescribed ROI regions
%     figure; plot(time(1:end-1),plotData(:,2)); hold on;scatter(time(keepLoc),plotData(keepLoc,2))

    % Temp-rate identification and sorting: 
    buffSize = 0.05;
    if any(ismember(strsplit(T.TempProtocol{trial},'_'),'sweeps')) %account for high freq jitter
        buffSize = 0.1;
    end
    for ii = 1:tPoints.nRates
        edges(ii,:) = [tPoints.rates(ii)-buffSize, tPoints.rates(ii)+buffSize];
    end  
    rateData = plotData(:,4);
    nRates = tPoints.nRates;
    rateIdx = discretize(rateData,nRates);
    for ii = 1:nRates
        TD = rateData(rateIdx==ii);
        T_rates(ii) = round(mean(TD,'omitnan'),2);
        % do these match the assumed rates?
        idx = find(T_rates(ii) > edges(:,1) & T_rates(ii) < edges(:,2));
        if isempty(idx)
            warndlg('Temp rate not within expected range')
            return
        end
        % set the rate value to the uniform cross-group value:
        rateData(rateIdx==ii) = tPoints.rates(ii);
        T_rates(ii) = tPoints.rates(ii);
%        figure; hold on; plot(rateData); plot(plotData(:,4))     
    end
    plotData(:,4) = rateData;
    plotData(:,5) = rateIdx; %rate bin
    TR.rateIdx = rateIdx;
    TR.rates = T_rates;
    TR.data = plotData;

    % Temperature range and bin size formatting
    t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
    if t_roi(end)<ceil(threshHigh)
        t_roi(end+1) = ceil(threshHigh) + binSpace;
    end
    nTemps = length(t_roi);
    tempIdx = discretize(plotData(:,1), t_roi);
    TR.nTemps = nTemps;
    TR.tempIdx = tempIdx;
    TR.temps = t_roi;
    
    % turn coordinates of heatmap into a vector to fill in 2D
    [HM,FC,FD,dist_mat,FC_mat] = deal([]);
    flyCount = data(trial).occupancy.flyCount;
    if size(flyCount,2)>1 %aka new data with white arena
        flyCount = flyCount(:,Alphabet(T.Arena{trial}));
    end
    for col = 1:nTemps
        for row = 1:nRates
        % Pull the data that matches this category
        loc = (rateIdx==row) & (tempIdx==col);

        % fly count data: (for tracking assessment)
        fc = flyCount(loc);
        FC(row,col).data = fc;
        FC(row,col).mean = mean(fc,'omitnan');
        FC(row,col).err = std(fc,0,'omitnan');
        FC_mat.avg(row,col) = FC(row,col).mean;
        FC_mat.err(row,col) = FC(row,col).err;

        % distance data:
        dist = plotData(loc,3);
        HM(row,col) = mean(dist,'omitnan'); %heat map
        FD(row,col).data = dist;
        FD(row,col).mean = mean(dist,'omitnan');
        FD(row,col).err = std(dist,0,'omitnan');
        dist_mat.avg(row,col) = FD(row,col).mean;
        dist_mat.err(row,col) = FD(row,col).err;
        end
    end
    TR.heatmap = HM;
    TR.dist_mat = dist_mat;
    TR.FC_mat = FC_mat;
    TR.nRates = length(TR.rates);
    G(trial).TR = TR;

clearvars('-except',vars{:}) 
end






%% FIGURE: temp rate hysteresis 

%% FIGURE: grouped -- temp distance hysteresis figures 
% TODO: reformat this to take advantage of the previously processed data
% Fly by fly average:

% Find the total number and id of temp rates:
allRates=[];
for trial = 1:ntrials
    allRates = [allRates,G(trial).TR.rates];
end
tRates = sort(unique(allRates));
nRates = length(tRates);
nTemps = G(1).TR.nTemps; %these should all be the same since they're held constant above

% Group the data for each temp rate
tempData = nan(nRates,nTemps,ntrials);
for trial = 1:ntrials
    for rr = 1:nRates
        idx = find(G(trial).TR.rates==tRates(rr));
        if isempty(idx)
            continue
        end
        tempData(rr,:,trial) = G(trial).TR.heatmap(idx,:);
    end
end

% Find the mean and err of each temp bin:
plotData.avg = mean(tempData,3,'omitnan');
plotData.err = std(tempData,0,3,'omitnan')./sqrt(ntrials);

heatMapData = plotData.avg;
t_roi = G(trial).TR.temps;
title_str = [ExpGroup ' (n = ' num2str(ntrials) ')'];

% ========= HeatMap of dT/dt vs T =============
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
    try ax.YTickLabel = tRates(YtickNum);
        
    catch
        ax.YTick = 1:length(tRates);
        ax.YTickLabel = tRates;
    end
    ylabel('\DeltaT/dt (\circC/min)')
    xlabel('Temp (\circC)')
    % Colorbar formatting
    cbh = colorbar(); 
    cbh.Label.String = 'Distance from food (mm)';
    cbh.Color = Color('white');
    % flip colormap around to make yellow closer to food
    cmap = colormap;
    set(gca, 'colormap', flip(cmap))

save_figure(fig, [figDir ExpGroup ' temp temp_rate dist heatmap ' ExpGroup], '-png');

% ========== Line plots of each rate comparison ==============
LS = {'--','-.','-'}; %cooling|stationary|heating

fig = figure;
hold on
for rr = 1:nRates
    if tRates(rr)>0
        lstyle = LS{3};
    elseif tRates(rr)<0
        lstyle = LS{1};
    elseif tRates(rr)==0
        continue
        lstyle = LS{2};
    end
    kolor = pullFoodColor(tRates(rr));
    x = t_roi;
    y = plotData.avg(rr,:);
    y_err = plotData.err(rr,:);
    plot(x,y,'color', kolor, 'linewidth', 2.5, 'linestyle', lstyle)%Color(cList{rr})
    if nRates<4
        plot(x,y+y_err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
        plot(x,y-y_err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
    end
end
title(ExpGroup)
ylabel('Distance to food (mm)')
xlabel('Temperature (\circC)')
formatFig(fig, true);
% legend:
idx = 0; str = [];
for ii = 1:nRates
    idx = idx+1;
    str{idx} = [num2str(tRates(ii)) '\circC/min'];
    if nRates<4
        str{idx+1} = ''; 
        str{idx+2} = ''; 
        idx = idx+2;
    end
end
legend(str,'textcolor', 'w', 'location', 'southwest', 'box', 'off')
%Save figure
save_figure(fig, [figDir ExpGroup ' temp temp_rate dist all rates ' ExpGroup], '-png');



% ========== Line plots of separated by rate MANUAL ADJUST FOR MORE RATES ==============
if nRates > 3

    LS = {'--','-.','-'}; %cooling|stationary|heating
    % legStr = {'SEM','Cooling','','','SEM','Heating','',''};
    legStr = {'Cooling','','','Heating','',''};
    ncol = 2;
    nrow = 2;

    fig = figure; set(fig, 'pos',[296 35 1224 956])
    % All trials
        subplot(nrow, ncol, 1); hold on
        for rr = 1:nRates
            kolor = pullFoodColor(tRates(rr));
            if tRates(rr)>0
                lstyle = LS{3};
            elseif tRates(rr)<0
                lstyle = LS{1};
            elseif tRates(rr)==0
                lstyle = LS{2};
                continue
            end
            x = t_roi;
            y = plotData.avg(rr,:);
            plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
        end
        title([ExpGroup ' (n = ' num2str(ntrials) ')'])
        ylabel('Distance to food (mm)')
        xlabel('Temperature (\circC)')
    % Fast
        subplot(nrow, ncol, 2); hold on
        for rr = [1,7]
            kolor = pullFoodColor(tRates(rr));
            if tRates(rr)>0
                lstyle = LS{3};
            elseif tRates(rr)<0
                lstyle = LS{1};
            elseif tRates(rr)==0
                lstyle = LS{2};
                continue
            end
            x = t_roi(1:end-1);
            y = plotData.avg(rr,1:end-1);
    %         kolor = Color(cList{rr});
            err = plotData.err(rr,1:end-1);
    %         fill_data = error_fill(x, y, err);
    %         h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
    %         set(h, 'facealpha', 0.2)
            plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
            plot(x,y+err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
            plot(x,y-err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
        end
        title(['dT/dt = ' num2str(tRates(rr)) ' (\circC/min)'])
        ylabel('Distance to food (mm)')
        xlabel('Temperature (\circC)')
        l = legend(legStr);
        set(l,'textcolor', 'w','box', 'off')
    % Medium
        subplot(nrow, ncol, 3); hold on
        for rr = [2,6]
            if tRates(rr)>0
                lstyle = LS{3};
            elseif tRates(rr)<0
                lstyle = LS{1};
            elseif tRates(rr)==0
                lstyle = LS{2};
                continue
            end
            x = t_roi(1:end-1);
            y = plotData.avg(rr,1:end-1);
    %         kolor = Color(cList{rr});
            kolor = pullFoodColor(tRates(rr));
            err = plotData.err(rr,1:end-1);
    %         fill_data = error_fill(x, y, err);
    %         h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
    %         set(h, 'facealpha', 0.2)
            plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
            plot(x,y+err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
            plot(x,y-err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
        end
        title(['dT/dt = ' num2str(tRates(rr)) ' (\circC/min)'])
        ylabel('Distance to food (mm)')
        xlabel('Temperature (\circC)')
        l = legend(legStr);
        set(l,'textcolor', 'w','box', 'off')
    % Slow
        subplot(nrow, ncol, 4); hold on
        for rr = [3,5]
            if tRates(rr)>0
                lstyle = LS{3};
            elseif tRates(rr)<0
                lstyle = LS{1};
            elseif tRates(rr)==0
                lstyle = LS{2};
                continue
            end
            x = t_roi(1:end-1);
            y = plotData.avg(rr,1:end-1);
    %         kolor = Color(cList{rr});
            kolor = pullFoodColor(tRates(rr));
            err = plotData.err(rr,1:end-1);
    %         fill_data = error_fill(x, y, err);
    %         h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
    %         set(h, 'facealpha', 0.2)
            plot(x,y,'color', kolor, 'linewidth', 2, 'linestyle', lstyle)
            plot(x,y+err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
            plot(x,y-err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
        end
        title(['dT/dt = ' num2str(tRates(rr)) ' (\circC/min)'])
        ylabel('Distance to food (mm)')
        xlabel('Temperature (\circC)')
        l = legend(legStr);
        set(l,'textcolor', 'w','box', 'off')

    fig = formatFig(fig, true,[nrow,ncol]);
    save_figure(fig, [figDir ExpGroup ' temp rate food preference dependence ' ExpGroup], '-png');
end
clearvars('-except',vars{:})

%% FIGURE: hysteresis across different foods and summed for all rates

% find the 'categories' of food:
[foodIdx,foodBin] = deal([]);
bin_foods = {'Plant', 'Molasses', 'Glucose', 'Water', 'ACV', 'Sugar','Merlot'};
for ii = 1:nfoods
    for jj = 1:length(bin_foods)
        if strfind(foodNames{ii}, bin_foods{jj})
            foodIdx(ii,1) = jj;
            foodBin{ii,1} = bin_foods{jj};
        end
    end
end
foodT = table(foodNames,foodIdx,foodBin);
disp(foodT)

% Find the total number and id of temp rates:
allRates=[];
for trial = 1:ntrials
    allRates = [allRates,G(trial).TR.rates];
end
tRates = sort(unique(allRates));
nRates = length(tRates);
nTemps = G(1).TR.nTemps; %these should all be the same since they're held constant above
abs_rates = unique(abs(tRates));

%allocate empty data structure & set params
HM = struct;
for ii = 1:size(foodT,1)
    % Group the data for each temp rate & food type:
    tempData = nan(nRates,nTemps,1);
    loc = 0;
    for trial = 1:ntrials
        % is this trial within the food bin?
        food_ID = data(trial).wellLabels{~strcmp(data(trial).wellLabels,'Empty')};
        N = foodT.foodIdx(strcmp(food_ID,foodNames));
        % Add data to structure if it is
        if N==ii
            loc = loc + 1;
            for rr = 1:nRates
                idx = find(G(trial).TR.rates==tRates(rr));
                if isempty(idx)
                    continue
                end
                tempData(rr,:,loc) = G(trial).TR.heatmap(idx,:);
            end
        end
    end
    % Find the mean and err of each temp bin for this food:
    HM(ii).avg = mean(tempData,3,'omitnan');
    HM(ii).err = std(tempData,0,3,'omitnan')./sqrt(loc);
    HM(ii).data = plotData.avg;
    HM(ii).name = foodT.foodBin{foodT.foodIdx==ii};
end

% FIGURE: same temp rates compared across foods:
LS = {'--','-.','-'}; %cooling|stationary|heating
x = G(trial).TR.temps;

for rr = 1:length(abs_rates)
  fig = figure; hold on; set(fig, 'pos', [57 118 598 642]);
    idx = 1; str = [];
    for ii = 1:size(foodT,1) % cycle through each food group
        rateList = find(abs(tRates)==abs_rates(rr));
        for jj = rateList
            plotRate = tRates(jj);
            if plotRate>0
                lstyle = LS{3};
                str_tag = 'heating';
            elseif plotRate<0
                lstyle = LS{1};
                str_tag = 'cooling';
            elseif plotRate==0
                lstyle = LS{2};
            end
            kolor = pullFoodColor(foodT.foodBin{ii});
            y = HM(ii).avg(jj,:);
            y_err = HM(ii).err(jj,:);
            plot(x,y,'color', kolor, 'linewidth', 2.5, 'linestyle', lstyle);
            str{idx} = [foodT.foodBin{ii} ' ' str_tag]; idx = idx+1;
        end
    end
    % legends, labels, keys
    title(['\pm' num2str(abs_rates(rr)) '\circC/min temp ramps'])
    ylabel('Distance to food (mm)')
    xlabel('Temperature (\circC)')
    formatFig(fig, true);
    legend(str,'textcolor', 'w', 'location', 'northeast', 'box', 'off','fontsize', 10)
%Save figure
save_figure(fig, [figDir 'temp hysteresis ' foodT.foodBin{ii} ' food ' num2str(abs_rates(rr))], '-png');
end



%%
    
    

title(ExpGroup)
ylabel('Distance to food (mm)')
xlabel('Temperature (\circC)')
formatFig(fig, true);
% legend:
idx = 0; str = [];
for ii = 1:nRates
    idx = idx+1;
    str{idx} = [num2str(tRates(ii)) '\circC/min'];
    if nRates<4
        str{idx+1} = ''; 
        str{idx+2} = ''; 
        idx = idx+2;
    end
end
legend(str,'textcolor', 'w', 'location', 'southwest', 'box', 'off')










heatMapData = plotData.avg;
t_roi = G(trial).TR.temps;
title_str = [ExpGroup ' (n = ' num2str(ntrials) ')'];

% ========= HeatMap of dT/dt vs T =============
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
    try ax.YTickLabel = tRates(YtickNum);
        
    catch
        ax.YTick = 1:length(tRates);
        ax.YTickLabel = tRates;
    end
    ylabel('\DeltaT/dt (\circC/min)')
    xlabel('Temp (\circC)')
    % Colorbar formatting
    cbh = colorbar(); 
    cbh.Label.String = 'Distance from food (mm)';
    cbh.Color = Color('white');
    % flip colormap around to make yellow closer to food
    cmap = colormap;
    set(gca, 'colormap', flip(cmap))

save_figure(fig, [figDir ExpGroup ' temp temp_rate dist heatmap ' ExpGroup], '-png');

% ========== Line plots of each rate comparison ==============
LS = {'--','-.','-'}; %cooling|stationary|heating

fig = figure;
hold on
for rr = 1:nRates
    if tRates(rr)>0
        lstyle = LS{3};
    elseif tRates(rr)<0
        lstyle = LS{1};
    elseif tRates(rr)==0
        continue
        lstyle = LS{2};
    end
    kolor = pullFoodColor(tRates(rr));
    x = t_roi;
    y = plotData.avg(rr,:);
    y_err = plotData.err(rr,:);
    plot(x,y,'color', kolor, 'linewidth', 2.5, 'linestyle', lstyle)%Color(cList{rr})
    if nRates<4
        plot(x,y+y_err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
        plot(x,y-y_err,'color', kolor, 'linewidth', 0.5, 'linestyle', lstyle)
    end
end
title(ExpGroup)
ylabel('Distance to food (mm)')
xlabel('Temperature (\circC)')
formatFig(fig, true);
% legend:
idx = 0; str = [];
for ii = 1:nRates
    idx = idx+1;
    str{idx} = [num2str(tRates(ii)) '\circC/min'];
    if nRates<4
        str{idx+1} = ''; 
        str{idx+2} = ''; 
        idx = idx+2;
    end
end
legend(str,'textcolor', 'w', 'location', 'southwest', 'box', 'off')
%Save figure
save_figure(fig, [figDir ExpGroup ' temp temp_rate dist all rates ' ExpGroup], '-png');






inputVar =  questdlg('Which data type to compare?','','distance','occupation probability','Cancel','distance');

food = struct;
switch inputVar
    case 'distance'
        ylab = 'Distance from well (mm)';
        L_loc = 'southwest';
    case 'occupation probability'
        ylab = inputVar;
        L_loc = 'northwest';
    case 'Cancel'
        return
end

%allocate empty data structure & set params
for ii = 1:nfoods+1
    food(ii).N = [];
    if ii <= nfoods
        food(ii).name = foodNames{ii};
    else 
        food(ii).name = 'Empty';
    end
end
for trial = 1:ntrials
    x = data(trial).occupancy.temp;
    for well = 1:4
        switch inputVar
            case 'distance'
                y = data(trial).dist2wells(:,well);
            case 'occupation probability'
                y = data(trial).occupancy.occ(:,well);
        end
        % empty food well:
        wellID = data(trial).wellLabels{well};
        if contains(wellID,'Empty')
            food(nfoods+1).N = [food(nfoods+1).N; x, y];
        else 
            loc = find(strcmp(wellID,foodNames));
            food(loc).N = [food(loc).N; x, y];
        end
    end
end 

% Temperature range and bin size formatting
[threshHigh, threshLow] = getTempThresholds;
binSpace = str2double(cell2mat(inputdlg('Bin size for temperature?','',[1,35],{'1'}))); 
t_roi = floor(threshLow):binSpace:ceil(threshHigh); 
if t_roi(end)<ceil(threshHigh)
    t_roi(end+1) = ceil(threshHigh) + binSpace;
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
    
% FIGURE: plot the avg. temp vs. distance data
CList = Color('SteelBlue', 'Navy', nfoods);
CList(nfoods+1,:) = Color('White');

fig = figure; set(fig,'pos', [67 82 675 692]);
hold on
for ii = 1:nfoods+1
    kolor = CList(ii,:);
    x = t_roi;
    y = food(ii).avg;
%     y_err = food(ii).err;
%     fill_data = error_fill(x,y, y_err);
%     h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
%     set(h, 'facealpha', 0.35)
    plot(x, y, 'Color', kolor, 'LineWidth', 2)
%     plot(x, y+y_err, 'Color', kolor, 'LineWidth', 0.25)
%     plot(x, y-y_err, 'Color', kolor, 'LineWidth', 0.25)
end
%Labels:
% xlim([7,20])
xlabel('temperature (\circC)')
ylabel(ylab)
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);
ax = gca;
set(ax, 'FontSize', 18)
str = strrep([foodNames; 'Empty'],'_',' ');
legend(str,'textcolor', 'w', 'location', L_loc, 'box', 'off','fontsize',12)

save_figure(fig, [figDir ExpGroup ' temp vs ' inputVar ' bin size ' num2str(binSpace)], '-png');
clearvars('-except',initial_vars{:})














































