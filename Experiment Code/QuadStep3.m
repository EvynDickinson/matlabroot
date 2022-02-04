
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

else %select data structure(s) from folder names out of GroupDataGUI:
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
    data = struct('dist2wells', n, 'wellLabels', n, 'occupancy', n,'nflies',n); %'data', n, 
    var2load = fieldnames(data);
    
    for trial = 1:ntrials
        filePath = [baseFolder, dates{trial}, '/Arena ' arenas{trial} '/analysis/'];
        temp = load([filePath expID{trial} arenas{trial} ' timecourse data.mat'], var2load{:});
        for ii = 1:length(var2load)
            try data(trial).(var2load{ii}) = temp.(var2load{ii});
            catch; data(trial).(var2load{ii}) = [];
            end
        end
        disp([expID{trial} arenas{trial}])
    end
end

pix2mm = 12.8; %conversion from pixels to mm for these videos
    
initial_vars = {'ExpGroup','baseFolder', 'T', 'data', 'figDir', 'filePath',...
                'initial_vars', 'folder', 'ntrials', 'pix2mm'};
clearvars('-except',initial_vars{:})
if questdlg('Save loaded data?')
    save([figDir ExpGroup ' raw'])
end
fprintf('Data loaded\n')

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
formatFig(fig, true)

save_figure(fig, [figDir ExpGroup ' temperature alignment'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% FIGURE: Distance to food vs temperature -- ALL trial lines
sSpan = 240;
fig = figure; hold on
for trial = 1:ntrials
    X = data(trial).occupancy.temp;
    for well = 1:4
        try Y = (data(trial).dist2wells(well).N(:,1));
        catch 
            Y = data(trial).occupancy.dist2wells(well).N(:,1);
            data(trial).dist2wells = data(trial).occupancy.dist2wells;
        end
        Y = Y./pix2mm; % convert the pixel values to mm
        % reorder the data by temperature:
        [plotX, idx] = sort(X);
        plot(plotX, smooth(Y(idx),sSpan), 'linewidth', 1,'color', pullFoodColor(data(trial).wellLabels{well}))
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

%% FIGURE: Average Distance vs temp across trials:

sSpan = 240;
[plant, yeast, empty] = deal([]);
for trial = 1:ntrials
    for well = 1:4
        x = data(trial).occupancy.temp';
        try y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        catch; y = (data(trial).occupancy.dist2wells(well).N(:,1))./pix2mm;
        end
        % plant food well:
        if strfind(data(trial).wellLabels{well},'Plant')
            plant = [plant; x, y];
        end
        % yeast food well:
        if strfind(data(trial).wellLabels{well},'Yeast')
            yeast = [yeast; x, y];
        end
        % empty food well:
        if strfind(data(trial).wellLabels{well},'Empty')
            empty = [empty; x, y];
        end 
    end
end 

FitType = questdlg('Select data for linear fit:', '', 'Smoothed data', 'All data','Cancel', 'All data');
numFitPoints = 1000;
[threshHigh, threshLow] = getTempThresholds;

emptyFood = false(1,3);
plotData = [];
% Fit the data:
for ii = 1:3
    [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); 
    switch ii
        case 1 % plant food
            inputData = plant;
            kolor = pullFoodColor('Plant');
        case 2 % yeast food
            inputData = yeast;
            kolor = pullFoodColor('Yeast');
        case 3 % empty
            inputData = empty;
            kolor = pullFoodColor('Empty');
    end
    %skip empty data sets (e.g. no yeast in trials)
    if isempty(inputData)
        emptyFood(ii) = true;
        continue; 
    end
    % cut off the high and low ends of data (to clean):
    loc = inputData(:,1)>threshHigh | inputData(:,1)<threshLow;
    inputData(loc,:) = [];
    % sort all the data by temperature:
    [~,idx] = sort(inputData(:,1));
    sortedData = inputData(idx,:);
    smoothData(:,1) = smooth(sortedData(:,1),sSpan);
    smoothData(:,2) = smooth(sortedData(:,2),sSpan);
    switch FitType
        case 'Smoothed data'
            fitdata = smoothData;
        case 'All data'
            fitdata = sortedData;
    end
    % find the line of best fit:
    coefficients = polyfit(fitdata(:,1), fitdata(:,2),1);
    xFit = linspace(fitdata(1,1), fitdata(end,1), numFitPoints);
    yFit = polyval(coefficients, xFit);
    yfit = polyval(coefficients, fitdata(:,1));
    yResid = fitdata(:,2)-yfit;
    SSresid = sum(yResid.^2);
    SStotal = (length(fitdata(:,2))-1) * var(fitdata(:,2));
    rsq = 1 - SSresid/SStotal; % get the R-square value
    % coefficient of correlation
    R = corrcoef(fitdata(:,1),fitdata(:,2));
    R = R(2,1);
    % save data for plotting:
    plotData(ii).smoothed = smoothData;
    plotData(ii).bestfit = [fitdata(:,1),yfit];
    plotData(ii).color = kolor;
    plotData(ii).rsq = round(rsq,2);
    plotData(ii).R = round(R,2);
end

% FIGURE: 
LW = 2;
fig = getfig; hold on
for ii = 1:3
    if emptyFood(ii); continue; end
    scatter(plotData(ii).smoothed(:,1), plotData(ii).smoothed(:,2), 30, plotData(ii).color)
end
for ii = 1:3
    if emptyFood(ii); continue; end
%     plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'r', 'linewidth', LW+1)
    plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'w', 'linewidth', LW)
end
ylim([10,40])

% Labels:
xlabel('temperature (\circC)')
ylabel('distance from well (mm)')
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);
ax = gca;
set(ax, 'FontSize', 25)

% label key:
foodLab = {'Plant', 'Yeast', 'Empty'};
printLoc = [8.5,16; 8.5,15; 8.5,14];
for ii = 1:3
  if ~emptyFood(ii)
    str = [foodLab{ii} ' :   R = ' num2str(plotData(ii).R) '  R^2 = ' num2str(plotData(ii).rsq)];
    text(printLoc(ii,1),printLoc(ii,2), str, 'Color', plotData(ii).color, 'FontSize', 12);
  end
end

save_figure(fig, [figDir ExpGroup ' avg temp vs distance'], '-png');
clearvars('-except',initial_vars{:})

%% FIGURE: timecourse of well distance % CAN'T DISTINGUISH WITHIN FOOD TYPE (E.G. PLANT VS PLANT)
% ADD IN THE TEMPERATURE LINE AS WELL
sSpan = 240;
nrows = 3;
ncols = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;


% extract data to plot:
timeLen = zeros(1,ntrials);
for trial = 1:ntrials
    timeLen(trial) = length(data(trial).occupancy.flycount);
end
len = max(timeLen); %max number of frames (or time) within the selected experiments
uni_time = linspace(1,(len/3)/60,len);

idx = 0;
pdata = [];
nEmpty = sum(ismember(data(1).wellLabels,'Empty'));
if nEmpty==3
    [pdata(1).Y,pdata(2).Y] = deal(nan([len,ntrials]));
    pdata(3).Y = deal(nan([len,ntrials*3])); %ASSUMES 3X EMPTY TO TEST
else 
    [pdata(1).Y,pdata(2).Y] = deal(nan([len,ntrials]));
    pdata(3).Y = deal(nan([len,ntrials*2])); %ASSUMES 2X EMPTY TO TEST
end

for well = 1:4
    for trial = 1:ntrials
        try Y = data(trial).dist2wells(well).N(:,1);
        catch; Y = data(trial).occupancy.dist2wells(well).N(:,1); 
        end
        Y = Y./pix2mm; % convert the pixel values to mm
        Y = smooth(Y,sSpan);
        % sort based on contents not well#
        if contains(data(trial).wellLabels{well},'Yeast','IgnoreCase', true) %yeast:
            pdata(1).Y(1:length(Y),trial) = Y;
            pdata(1).kolor = pullFoodColor('Yeast');
        elseif contains(data(trial).wellLabels{well},'Plant','IgnoreCase', true) %yeast:
            pdata(2).Y(1:length(Y),trial) = Y;
            pdata(2).kolor = pullFoodColor('Plant');
        elseif contains(data(trial).wellLabels{well},'Empty','IgnoreCase', true) %yeast:
            idx = idx+1;
            pdata(3).Y(1:length(Y),idx) = Y;
            pdata(3).kolor = pullFoodColor('Empty');
        end
    end
end
for ii = 1:3 %well contents type
    % average + error across flies: (excludes data points not covered in all trials)
    Y_all = pdata(ii).Y;
    Y_avg = nanmean(Y_all,2);
    Y_err = std(Y_all,0,2);
    nanloc = isnan(Y_err);
    Y_err(nanloc) = [];
    Y_avg(nanloc) = [];
    pdata(ii).y_err = Y_err;
    pdata(ii).y_avg = Y_avg;
end

% FIGURE WITH TEMP VS TIME | DISTANCE VS TIME
fig = figure; set(fig, 'pos', [96 444 660 482]);
subplot(nrows, ncols, sb(1).idx)
    hold on
    for trial = 1:ntrials
        X = data(trial).occupancy.time;
        Y = data(trial).occupancy.temp;
        plot(X, Y, 'linewidth', 1)
    end
    xlabel('Time (min)')
    ylabel('Temp (\circ)')
    title({'Distance from food sources';...
          ['N = ' num2str(ntrials)]})

subplot(nrows, ncols, sb(2).idx)
    hold on
    % error fills
    for ii = 1:3
        kolor = pdata(ii).kolor;
        y_avg = pdata(ii).y_avg;
        y_err = pdata(ii).y_err;
        time = uni_time(1:length(y_avg));
        fill_data = error_fill(time, y_avg, y_err);
        h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
        set(h, 'facealpha', 0.4)
    end
    % plot avg line
    for ii = 1:3
        if isempty(pdata(ii).kolor)
            continue
        end
        plot(uni_time(1:length(y_avg)), pdata(ii).y_avg, 'linewidth', 2,'color', pdata(ii).kolor)
    end
    ylabel('Distance (mm)')
    xlabel('Time (min)')   
formatFig(fig, true,[nrows,ncols],sb);
    
save_figure(fig, [figDir ExpGroup ' avg distance vs time'], '-png');    


clearvars('-except',initial_vars{:})
fprintf('Next\n')


%% FIGURE: Avg movement vs. temperature:
sSpan = 360;
timeLen = zeros(1,ntrials);
for trial = 1:ntrials
    timeLen(trial) = length(data(trial).occupancy.flycount);
end
len = max(timeLen); %max number of frames (or time) within the selected experiments
uni_time = linspace(1,(len/3)/60,len); % universal time component...

[mov, temp] = deal(nan([len,ntrials]));
all = [];
for trial = 1:ntrials
    x = data(trial).occupancy.temp(1:end-1)';
    y = data(trial).occupancy.movement; 
    % movement data
    mov(1:length(y),trial) = y;
    % temperature data
    temp(1:length(x),trial) = x;
    % both:
    all = [all; x,y];
end
% sort the 'all' data by temperature dependence:
[~, idx] = sort(all(:,1));
sortedData = all(idx,:);
% avg and err of all the trials
Y_avg = nanmean(mov,2);
mov_avg = nanmean(Y_avg);
Y_err = std(mov,0,2);
nanloc = isnan(Y_err);
Y_err(nanloc) = [];
Y_avg(nanloc) = [];

[threshHigh, threshLow] = getTempThresholds;
xlimits = [floor(threshLow), ceil(threshHigh)];

 %change this if the temp range changes drastically!
nrows = 1;
ncols = 4;
sb(1).idx = 1:2;
sb(2).idx = 3:4;
kolor = Color('steelblue');

% FIGURE
fig = figure; set(fig,'pos',[49 200 1491 541]);
% TIMECOURSE OF MOVEMENT
subplot(nrows, ncols, sb(1).idx)
    hold on
    time = uni_time(1:length(Y_avg));
    fill_data = error_fill(time, Y_avg, Y_err);
    h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
    set(h, 'facealpha', 0.4)
    plot(time, Y_avg, 'linewidth',0.5,'color', kolor)
    hline(mov_avg, 'w:')
    ylabel('Movement (au)')
    xlabel('Time (min)')
    title(['Avg movement: ' num2str(mov_avg)])
    ylim([-5,20])
% TEMP DEPENDENCE
subplot(nrows, ncols, sb(2).idx) 
    hold on
    scatter(smooth(sortedData(:,1),sSpan),smooth(sortedData(:,2),sSpan),10,kolor,'filled')
    xlim(xlimits)
    xlabel('Temperature (\circC)')
    ylabel('Movement (au)')
fig = formatFig(fig, true, [nrows,ncols],sb);

save_figure(fig, [figDir ExpGroup ' movement vs time+temp'], '-png');    


clearvars('-except',initial_vars{:})
fprintf('Next\n')

% TODO: 
% add an analysis that looks at how different the line is from zero slope

%% FIGURE: Clustering vs temperature 
ROI = 1:40490; % take the first 225 minutes only | [] =  all data

if ntrials > 15
    warndlg('There are too many data points for this analysis'); return
end
sSpan = 360;
timeLen = zeros(1,ntrials);
for trial = 1:ntrials
    timeLen(trial) = length(data(trial).occupancy.flycount);
end
len = max(timeLen); %max number of frames (or time) within the selected experiments
uni_time = linspace(1,(len/3)/60,len); % universal time component...

sSpan = 360;
[clust, temp] = deal(nan([len,ntrials]));
all = [];
for trial = 1:ntrials
    
    x = data(trial).occupancy.temp';
    y = data(trial).occupancy.IFD./pix2mm;  % interfly distance
    
    if ~isempty(ROI)
        x = x(ROI);
        y = y(ROI);
    end

    %normalize to the minimum distance (1=most clustered)...this helps deal with the
    %fly density variability
    minY = min(y);
    maxY = max(y);
    z = abs(1-((y-minY)./(maxY-minY))); %1=most clustered 0 = no clustering
    
    % clustering data
    clust(1:length(z),trial) = z;
    % temperature data
    temp(1:length(x),trial) = x;
    % both:
    all = [all; x,z'];
end
% sort the 'all' data by temperature dependence:
[~, idx] = sort(all(:,1));
sortedData = all(idx,:);
% avg and err of all the trials
Y_avg = nanmean(clust,2);
Y_err = std(clust,0,2);
nanloc = isnan(Y_err);
Y_err(nanloc) = [];
Y_avg(nanloc) = [];

xlimits = [8,20]; %TODO: change this if the temp range changes drastically!
ylimits = [0,1];
nrows = 1;
ncols = 4;
sb(1).idx = 1:2;
sb(2).idx = 3:4;
kolor = Color('slateblue');


% FIGURE
fig = figure; set(fig,'pos',[49 200 1491 541]);
% TIMECOURSE OF MOVEMENT
subplot(nrows, ncols, sb(1).idx)
    hold on
    time = uni_time(1:length(Y_avg));
    fill_data = error_fill(time, Y_avg, Y_err);
    h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
    set(h, 'facealpha', 0.4)
    plot(time, Y_avg, 'linewidth',0.5,'color', kolor)
    ylabel('Clustering Index')
    xlabel('Time (min)')
    ylim(ylimits)
% TEMP DEPENDENCE
subplot(nrows, ncols, sb(2).idx) 
    hold on
    scatter(smooth(sortedData(:,1),sSpan),smooth(sortedData(:,2),sSpan),10,kolor,'filled')
    xlim(xlimits)
    ylim(ylimits)
    xlabel('Temperature (\circC)')
    ylabel('Clustering Index')
fig = formatFig(fig, true, [nrows,ncols],sb);


save_figure(fig, [figDir ExpGroup ' clustering vs time+temp'], '-png');    


clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% FIGURE: Fly distance binned by Temperature
emptyFood = false(1,3);
[plant, yeast, empty] = deal([]);
for trial = 1:ntrials
    for well = 1:4
        x = data(trial).occupancy.temp';
        try y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        catch; y = (data(trial).occupancy.dist2wells(well).N(:,1))./pix2mm;
        end
        % plant food well:
        if strfind(data(trial).wellLabels{well},'Plant')
            plant = [plant; x, y];
        end
        % yeast food well:
        if strfind(data(trial).wellLabels{well},'Yeast')
            yeast = [yeast; x, y];
        end
        % empty food well:
        if strfind(data(trial).wellLabels{well},'Empty')
            empty = [empty; x, y];
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
plotData = [];
for ii = 1:3 %each food type
    switch ii
        case 1 % plant food
            inputData = plant;
            plotData(ii).kolor = pullFoodColor('Plant');
        case 2 % yeast food
            inputData = yeast;
            plotData(ii).kolor = pullFoodColor('Yeast');
        case 3 % empty
            inputData = empty;
            plotData(ii).kolor = pullFoodColor('Empty');
    end
    if isempty(inputData)
        emptyFood(ii) = true;
        continue
    end
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
        plotData(ii).y_avg(tt) = mean(inputData(loc,2));
    end
    plotData(ii).y_err = y_err./sqrt(ntrials);
    plotData(ii).xdata = t_roi(unique_a);
end
    
% PLOT
fig = figure; set(fig,'pos', [260 196 413 692]);
hold on
for ii = 1:3
    if emptyFood(ii)
        continue
    end
    fill_data = error_fill(plotData(ii).xdata, plotData(ii).y_avg, plotData(ii).y_err);
    h = fill(fill_data.X, fill_data.Y, plotData(ii).kolor, 'EdgeColor','none');
    set(h, 'facealpha', 0.35)
    plot(plotData(ii).xdata, plotData(ii).y_avg, 'Color', plotData(ii).kolor, 'LineWidth', 2)
    plot(plotData(ii).xdata, plotData(ii).y_avg+plotData(ii).y_err, 'Color', plotData(ii).kolor, 'LineWidth', 0.25)
    plot(plotData(ii).xdata, plotData(ii).y_avg-plotData(ii).y_err, 'Color', plotData(ii).kolor, 'LineWidth', 0.25)
end
%Labels:
% xlim([7,20])
xlabel('temperature (\circC)')
ylabel('distance from well (mm)')
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);
ax = gca;
set(ax, 'FontSize', 18)

save_figure(fig, [figDir ExpGroup ' temp_err vs well distance bin size ' num2str(binSpace)], '-png');
clearvars('-except',initial_vars{:})

%% FIGURE: Number of Flies Effect Size
sSpan = 900;
% determine the number of flies in the experiment population:
nflies = [];
for trial = 1:ntrials
    nflies(1,trial) = data(trial).nflies;
end
flyNums = unique(nflies);
inputdata = struct;
% edges = [1,8,13,18,26]; % bin edges
edges = [3,7,10:3:19,26]; % bin edges

% fig = figure; h = histogram(nflies);
% set(h, 'facecolor', Color('deeppink'),'edgecolor', 'k','facealpha', 1)
% xlabel('Number of flies in arena')
% ylabel('Count')
% v_line(edges-0.5,'deepskyblue','-',1)
% fig = formatFig(fig, true);
% set(gca, 'fontsize', 18)
% save_figure(fig, [figDir ExpGroup ' num flies histogram ' num2str(ngroup) ' bins'], '-png');

% binned fly numbers:
bIdx = discretize(nflies,edges);
ngroup = length(edges)-1;
nflies(2,:) = bIdx;

% all fly numbers:
for gg = 1:ngroup
    inputdata(gg).nflies = edges(gg):edges(gg+1);
    [inputdata(gg).plant,inputdata(gg).yeast,inputdata(gg).empty] = deal([]);
    leg{gg} = [num2str(edges(gg)) '-' num2str(edges(gg+1)-1) ' flies']; % legend for figure
end

% sort data by fly numbers
for trial = 1:ntrials
    gg = nflies(2,trial); %grouping number
    for well = 1:4
        % pull data
        x = data(trial).occupancy.temp';
        try y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        catch; y = (data(trial).occupancy.dist2wells(well).N(:,1))./pix2mm;
        end
        [~,num] = pullFoodColor(data(trial).wellLabels{well});
        switch num
            case 1 % PLANT
            inputdata(gg).plant = [inputdata(gg).plant; x,y];
            case 2 % YEAST
            inputdata(gg).yeast = [inputdata(gg).yeast; x,y];
            case 3 % EMPTY
            inputdata(gg).empty = [inputdata(gg).empty; x,y];
        end
    end
end

FitType = questdlg('Select data for linear fit:', '', 'Smoothed data', 'All data','Cancel', 'All data');

switch questdlg('Select linear fit cutoffs:', '', 'Low (8-20)', 'High (6-25)','Cancel', 'Low (8-20)')
    case 'Low (8-20)'
        threshHigh = 19.88;
        threshLow = 8.02;
    case 'High (6-25)'
        threshHigh = 24.5;
        threshLow = 6.5;
    case 'Cancel'
end
numFitPoints = 100;

% Fit the data:
plotData = [];
for gg = 1:ngroup
  for ii = 1:3
    [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); 
    switch ii
        case 1 % plant food
            raw = inputdata(gg).plant;
%             kolor = pullFoodColor('Plant');
        case 2 % yeast food
            raw = inputdata(gg).yeast;
%             kolor = pullFoodColor('Yeast');
        case 3 % empty
            raw = inputdata(gg).empty;
%             kolor = pullFoodColor('Empty');
    end
    % cut off the high and low ends of data (to clean):
    loc = raw(:,1)>threshHigh | raw(:,1)<threshLow;
    raw(loc,:) = [];
    % sort all the data by temperature:
    [~,idx] = sort(raw(:,1));
    sortedData = raw(idx,:);
    smoothData(:,1) = smooth(sortedData(:,1),sSpan);
    smoothData(:,2) = smooth(sortedData(:,2),sSpan);
    % find the line of best fit:
    switch FitType
        case 'Smoothed data'
            fitData = smoothData;
        case 'All data'
            fitData = sortedData;
    end
    coefficients = polyfit(fitData(:,1), fitData(:,2),1);
    xFit = linspace(fitData(1,1), fitData(end,1), numFitPoints);
    yFit = polyval(coefficients, xFit);
    yfit = polyval(coefficients, fitData(:,1));
    yResid = fitData(:,2)-yfit;
    SSresid = sum(yResid.^2);
    SStotal = (length(fitData(:,2))-1) * var(fitData(:,2));
    rsq = 1 - SSresid/SStotal; % get the R-square value
    % coefficient of correlation
    R = corrcoef(fitData(:,1),fitData(:,2));
    R = R(2,1);
    % save data for plotting:
%     plotData(gg,ii).rawData = sortedData;
    plotData(gg,ii).smoothed = smoothData;
    plotData(gg,ii).bestfit = [fitData(:,1),yfit];
%     plotData(ii).color = kolor;
    plotData(gg,ii).rsq = round(rsq,2);
    plotData(gg,ii).R = round(R,2);
    
    
  end
end
        
CList = {'Lavender','Thistle', 'Plum', 'MediumPurple', 'BlueViolet', 'Indigo'};

% PLOT the data:
LW = 3;
nrows = 1;
ncols = 3;
fig = figure; set(fig,'pos', [274 156 973 680]);
for ii = 1:3
    subplot(nrows,ncols,ii)
    hold on
    switch ii
        case 1
            color_mat = Color('PaleGreen', 'DarkGreen', ngroup);
            title('Plant Food')
        case 2
            color_mat = Color('LemonChiffon', 'DarkOrange', ngroup);
            title('Yeast Food')
        case 3
            color_mat = Color('White', 'Gray', ngroup);
            title('Empty')
    end
    for gg = 1:ngroup
        scatter(plotData(gg,ii).smoothed(:,1), plotData(gg,ii).smoothed(:,2), 30, Color(CList{gg}),'filled')
    %     scatter(plotData(gg,ii).rawData(:,1), plotData(gg,ii).rawData(:,2), 30, color_mat(gg,:))
    end
    for gg = 1:ngroup
        plot(plotData(gg,ii).bestfit(:,1),plotData(gg,ii).bestfit(:,2),...
            'color', 'w', 'linewidth', LW+2)
        plot(plotData(gg,ii).bestfit(:,1),plotData(gg,ii).bestfit(:,2),...
            'color', Color(CList{gg}), 'linewidth', LW)
    end    
    % Labels:
    ylim([15,40])
    xlabel('temperature (\circC)')
    if ii==1; ylabel('distance from well (mm)'); end
end
fig = formatFig(fig, true, [nrows,ncols]);
for ii = 3:3
    subplot(nrows,ncols,ii)
    l = legend(leg); set(l, 'TextColor', 'w');
end

save_figure(fig, [figDir ExpGroup ' temp-dist by num flies - ' num2str(ngroup) ' bins'], '-png');

% SLOPE OF BEST FIT LINES
% Find the slope of the line of best fit: 
for gg = 1:ngroup %number of flies
    for ii = 1:3 %food types
        x1 = plotData(gg,ii).bestfit(1,1);
        x2 = plotData(gg,ii).bestfit(end,1);
        y1 = plotData(gg,ii).bestfit(1,2);
        y2 = plotData(gg,ii).bestfit(end,2);
        plotData(gg,ii).m = (y2-y1)/(x2-x1);
        slopes(gg,ii) = plotData(gg,ii).m;
    end
end

% PLOT the slopes of the lines:
fig = figure; hold on
for gg = 1:ngroup
   p = plot(1:3, slopes(gg,:), 'color', Color(CList{gg}), 'linewidth', 2, 'linestyle','--');
   set(p, 'Marker','d','MarkerSize', 10,'MarkerFaceColor', Color(CList{gg}));
end
xlim([0.75,3.25])
ylim([-0.81,0.21])
set(gca, 'xtick', 1:3,'xticklabel', {'Plant', 'Yeast', 'Empty'})
xlabel({' '; 'Foods'})
ylabel('Slope of best-fit line')
fig = formatFig(fig, true);
l = legend(leg); set(l, 'TextColor', 'w');
set(gca, 'fontsize', 15)

save_figure(fig, [figDir ExpGroup ' slope of temp-dist by num flies with ' num2str(ngroup) ' bins'], '-png');


clearvars('-except',initial_vars{:})

%% FIGURE: Well identity vs temp & distance meta analysis
kolors = {'DeepPink', 'Orange', 'Lime', 'DodgerBlue'};
sSpan = 900;
[well_1, well_2, well_3, well_4] = deal([]);
for trial = 1:ntrials
    for well = 1:4
        x = data(trial).occupancy.temp';
        try y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        catch; y = (data(trial).occupancy.dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        end
        % well 1:
        switch well
            case 1
                well_1 = [well_1; x, y];
            case 2
                well_2 = [well_2; x, y];
            case 3
                well_3 = [well_3; x, y];
            case 4
                well_4 = [well_4; x, y];
        end 
    end
end 


% -------------------------- temperature : distance relationship --------------------

numFitPoints = 100;
switch questdlg('Select linear fit cutoffs:', '', 'Low (8-20)', 'High (6-25)','Cancel', 'Low (8-20)')
    case 'Low (8-20)'
        threshHigh = 19.88;
        threshLow = 8.02;
    case 'High (6-25)'
        threshHigh = 24.5;
        threshLow = 6.5;
    case 'Cancel'
end

FitType = questdlg('Select data for linear fit:', '', 'Smoothed data', 'All data','Cancel', 'All data');
plotData = [];
for ii = 1:4
    [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); 
    switch ii
        case 1 
            inputData = well_1;
        case 2 
            inputData = well_2;
        case 3 
            inputData = well_3;
        case 4 
            inputData = well_4;
    end
    kolor = Color(kolors{ii});
    % cut off the high and low ends of data (to clean):
    loc = inputData(:,1)>threshHigh | inputData(:,1)<threshLow;
    inputData(loc,:) = [];
    % sort all the data by temperature:
    [~,idx] = sort(inputData(:,1));
    sortedData = inputData(idx,:);
    smoothData(:,1) = smooth(sortedData(:,1),sSpan);
    smoothData(:,2) = smooth(sortedData(:,2),sSpan);
    % find the line of best fit:
    switch FitType
        case 'Smoothed data'
            rawData = smoothData;
        case 'All data'
            rawData = sortedData;
        case 'Cancel'
            return
    end
    coefficients = polyfit(rawData(:,1), rawData(:,2),1);
    xFit = linspace(rawData(1,1), rawData(end,1), numFitPoints);
    yFit = polyval(coefficients, xFit);
    yfit = polyval(coefficients, rawData(:,1));
    yResid = rawData(:,2)-yfit;
    SSresid = sum(yResid.^2);
    SStotal = (length(rawData(:,2))-1) * var(rawData(:,2));
    rsq = 1 - SSresid/SStotal; % get the R-square value
    % coefficient of correlation
    R = corrcoef(rawData(:,1),rawData(:,2));
    R = R(2,1);
    % save data for plotting:
    plotData(ii).smoothed = smoothData;
    plotData(ii).bestfit = [rawData(:,1),yfit];
    plotData(ii).color = kolor;
    plotData(ii).rsq = round(rsq,2);
    plotData(ii).R = round(R,2);
end

% FIGURE: 
LW = 2;
fig = getfig; hold on
for ii = 1:4
    scatter(plotData(ii).smoothed(:,1), plotData(ii).smoothed(:,2), 30, plotData(ii).color,'filled')
end
for ii = 1:4
    plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'w', 'linewidth', LW+1)
    plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', plotData(ii).color,...
         'linewidth', LW, 'LineStyle', '-.')
end
ylim([10,40])

% Labels:
xlabel('temperature (\circC)')
ylabel('distance from well (mm)')
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);
ax = gca;
set(ax, 'FontSize', 25)
l = legend({'Well 1', 'Well 2', 'Well 3', 'Well 4'});
set(l, 'TextColor', 'w')

% % label key:
% str = ['Plant :   R = ' num2str(plotData(1).R) '  R^2 = ' num2str(plotData(1).rsq)];
% text(8.5,16, str, 'Color', plotData(1).color, 'FontSize', 12);
% str = ['Yeast :  R = ' num2str(plotData(2).R) '  R^2 = ' num2str(plotData(2).rsq)];
% text(8.5,15, str, 'Color', plotData(2).color, 'FontSize', 12);
% str = ['Empty : R = ' num2str(plotData(3).R) '  R^2 = ' num2str(plotData(3).rsq)];
% text(8.5,14, str, 'Color', plotData(3).color, 'FontSize', 12);

save_figure(fig, [figDir ExpGroup ' avg temp vs well location'], '-png');
clearvars('-except',initial_vars{:})

%% FIGURE: Well identity histogram 
kolors = {'DeepPink', 'Orange', 'Lime', 'DodgerBlue'};
sSpan = 900;
[plant, yeast, empty] = deal([]);
for trial = 1:ntrials
    for well = 1:4
        labs = data(trial).wellLabels{well};
        if strfind(labs,'Plant')
            plant = [plant,well];
        end
        if strfind(labs, 'Yeast')
            yeast = [yeast,well];
        end
        if strfind(labs, 'Empty')
            empty = [empty,well];
        end
    end
end 

for well = 1:4
    well_ID(1,well) = sum(plant==well);
    well_ID(2,well) = sum(yeast==well);
    well_ID(3,well) = sum(empty==well);
end
%normalize
kolors = {'Plant', 'Yeast', 'Empty'};
wellID_per = round((well_ID./sum(well_ID,2))*100);
% PLOT
nrows = 1;
ncols = 2;

fig = figure;
% absolute count
subplot(nrows,ncols,1)
hold on
for ii = 1:3
    plot(1:4,well_ID(ii,:),'Marker', '*', 'Color', pullFoodColor(kolors{ii}), 'LineWidth', 2)
end
xlim([0,5])
ylim([0,45])
set(gca, 'XTick', 1:4)
xlabel('Well')
ylabel('Instances (#)')
% percent of instances
subplot(nrows,ncols,2)
hold on
for ii = 1:3
    plot(1:4,wellID_per(ii,:),'Marker', '*', 'Color', pullFoodColor(kolors{ii}), 'LineWidth', 2)
end
xlim([0,5])
ylim([10,40])
set(gca, 'XTick', 1:4)
xlabel('Well')
ylabel('Instances (%)')
% Labels and formatting
fig = formatFig(fig, true,[nrows,ncols]);
l = legend(kolors);
set(l, 'TextColor', 'w')
set(gca, 'FontSize', 15)

save_figure(fig, [figDir ExpGroup ' food identity per well'], '-png');

% % ---------------------------------
% % Sum of each within wells 1+4 vs 2+3:
% pair_1 = [1,4];
% pair_2 = [2,3];
% P_1 = sum(well_ID(:,pair_1),2);
% P_2 = sum(well_ID(:,pair_2),2);
% 
% k_list = {'green', 'gold', 'grey'};
% fig = figure;
% h = bar([P_1';P_2']);
% for ii = 1:3
%     h(ii).FaceColor = Color(k_list{ii});
% end
% ylabel('Instances')
% set(gca, 'XTickLabels', {'Wells 1 & 4', 'Wells 2 & 3'},'FontSize', 20)
% fig = formatFig(fig, true);
% 
% save_figure(fig, [figDir ExpGroup ' well counts paired'], '-png');

clearvars('-except',initial_vars{:})

%% FIGURE: Well identity by arena -- trends

kolors = {'Thistle', 'Orchid', 'MediumPurple', 'Indigo'};
sSpan = 900;
%build empty arena to stuff with data
G = struct('wells',cell(1,4));
for arena = 1:4
    G(arena).wells.well_1 = [];
    G(arena).wells.well_2 = [];
    G(arena).wells.well_3 = [];
    G(arena).wells.well_4 = [];
end

arenaList = T.Arena;
for trial = 1:ntrials
    % pull out data
    switch arenaList{trial}
        case 'A'
            inputdata = G(1).wells;
        case 'B'
            inputdata = G(2).wells;
        case 'C'
            inputdata = G(3).wells;
        case 'D'
            inputdata = G(4).wells;
    end
    % sort data
    for well = 1:4
        x = data(trial).occupancy.temp';
        try y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        catch; y = (data(trial).occupancy.dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        end
        switch well
            case 1
                inputdata.well_1 = [inputdata.well_1; x, y];
            case 2
                inputdata.well_2 = [inputdata.well_2; x, y];
            case 3
                inputdata.well_3 = [inputdata.well_3; x, y];
            case 4
                inputdata.well_4 = [inputdata.well_4; x, y];
        end 
    end
    % put data back into struct
    switch arenaList{trial}
        case 'A'
            G(1).wells = inputdata;
        case 'B'
            G(2).wells = inputdata;
        case 'C'
            G(3).wells = inputdata;
        case 'D'
            G(4).wells = inputdata;
    end
end 


% -------------------------- temperature : distance relationship --------------------

numFitPoints = 100;
switch questdlg('Select linear fit cutoffs:', '', 'Low (8-20)', 'High (6-25)','Cancel', 'Low (8-20)')
    case 'Low (8-20)'
        threshHigh = 19.88;
        threshLow = 8.02;
    case 'High (6-25)'
        threshHigh = 24.5;
        threshLow = 6.5;
    case 'Cancel'
end

FitType = questdlg('Select data for linear fit:', '', 'Smoothed data', 'All data','Cancel', 'All data');
plotData = [];
for arena = 1:4
  for well = 1:4
    [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); 
    switch well
        case 1 
            inputData = G(arena).wells.well_1;
        case 2 
            inputData = G(arena).wells.well_2;
        case 3 
            inputData = G(arena).wells.well_3;
        case 4 
            inputData = G(arena).wells.well_4;
    end
    % cut off the high and low ends of data (to clean):
    loc = inputData(:,1)>threshHigh | inputData(:,1)<threshLow;
    inputData(loc,:) = [];
    % sort all the data by temperature:
    [~,idx] = sort(inputData(:,1));
    sortedData = inputData(idx,:);
    smoothData(:,1) = smooth(sortedData(:,1),sSpan);
    smoothData(:,2) = smooth(sortedData(:,2),sSpan);
    % find the line of best fit:
    switch FitType
        case 'Smoothed data'
            rawData = smoothData;
        case 'All data'
            rawData = sortedData;
        case 'Cancel'
            return
    end
    coefficients = polyfit(rawData(:,1), rawData(:,2),1);
    xFit = linspace(rawData(1,1), rawData(end,1), numFitPoints);
    yFit = polyval(coefficients, xFit);
    yfit = polyval(coefficients, rawData(:,1));
    yResid = rawData(:,2)-yfit;
    SSresid = sum(yResid.^2);
    SStotal = (length(rawData(:,2))-1) * var(rawData(:,2));
    rsq = 1 - SSresid/SStotal; % get the R-square value
    % coefficient of correlation
    R = corrcoef(rawData(:,1),rawData(:,2));
    R = R(2,1);
    % save data for plotting:
    plotData(arena,well).smoothed = smoothData;
    plotData(arena,well).bestfit = [rawData(:,1),yfit];
    plotData(arena,well).rsq = round(rsq,2);
    plotData(arena,well).R = round(R,2);
    % slope of line
    x1 = plotData(arena,well).bestfit(1,1);
    x2 = plotData(arena,well).bestfit(end,1);
    y1 = plotData(arena,well).bestfit(1,2);
    y2 = plotData(arena,well).bestfit(end,2);
    plotData(arena,well).m = (y2-y1)/(x2-x1);
    slopes(arena,well) = plotData(arena,well).m;
  end
end

% -------------------------- PLOTS ---------------------------
nrows = 2;
ncols = 2;
sz = 10;
fig = figure;
for well = 1:4
    subplot(nrows, ncols, well)
    hold on
    for arena = 1:4
        x = plotData(arena,well).smoothed(:,1);
        y = plotData(arena,well).smoothed(:,2);
        scatter(x,y,sz,Color(kolors{arena}),'filled')
    end
    ylim([10,40])
    xlim([7,21])
    set(gca, 'xtick', [8,20])
    title(['Well ' num2str(well)])
end
% labels and formatting
subplot(nrows, ncols, 1)
    ylabel('Distance (mm)')
subplot(nrows, ncols, 3)
    ylabel('Distance (mm)')   
subplot(nrows, ncols, 3)
    xlabel('Temp (\circC)')
subplot(nrows, ncols, 4)
    xlabel('Temp (\circC)')     
fig = formatFig(fig, true, [nrows,ncols]);

save_figure(fig, [figDir ExpGroup ' well distance by arena'], '-png');
    

% PLOT the slopes of the lines:
fig = figure; hold on
for arena = 1:4
   p = plot(1:4, slopes(arena,:), 'color', Color(kolors{arena}), 'linewidth', 2, 'linestyle','--');
   set(p, 'Marker','d','MarkerSize', 10,'MarkerFaceColor', Color(kolors{arena}));
end
xlim([0.75,4.25])
ylim([-0.61,0.11])
set(gca, 'xtick', 1:4)
xlabel('Wells')
ylabel('Slope of best-fit line')
fig = formatFig(fig, true);
l = legend({'Arena A', 'Arena B', 'Arena C', 'Arena D'}); set(l, 'TextColor', 'w');

save_figure(fig, [figDir ExpGroup ' slope of well distance by arena'], '-png');

clearvars('-except',initial_vars{:})

%% FIGURE: individual distance line fit by food
% TODO: update this to work for variable food options (e.g. 1 plant, 3 empty)
% Currently only works for 
sSpan = 120;
numFitPoints = 20;
% Fit each trial individually and then plot them all...
FitType = questdlg('Select data for linear fit:', '', 'Smoothed data', 'All data','Cancel', 'All data');
[threshHigh, threshLow] = getTempThresholds;

%build empty arena to stuff with data
food = struct;
for K = 1:3 %foods
    food(K).slope = [];
    food(K).X = [];
    food(K).Y = [];
end

%Pull data for each food type: [the empty pull only works for 2 empty
%wells...]
for trial = 1:ntrials
    mt = 1;
    mt_data = [];
    nEmpty = sum(strcmp(data(trial).wellLabels,'Empty'));
    for well = 1:4
        smoothData = [];
        % pull data
        [~,K] = pullFoodColor(data(trial).wellLabels{well});
        x = data(trial).occupancy.temp';
        try y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        catch; y = (data(trial).occupancy.dist2wells(well).N(:,1))./pix2mm;
        end
        loc = (x<threshHigh&x>threshLow); % points within temp bounds
        x = x(loc);
        y = y(loc);
        raw = [x,y];
        
        if K==3 %empty
            mt = mt+1;
            if mt<=nEmpty %second round of empty
                mt_data = [mt_data; raw];
                continue %skip to the next loop
            end
            if mt==nEmpty
                raw = mt_data; %pull the data back out to use...
            end
        end

        % fit the data
        [~,idx] = sort(raw(:,1));
        sortedData = raw(idx,:);
        smoothData(:,1) = smooth(sortedData(:,1),sSpan);
        smoothData(:,2) = smooth(sortedData(:,2),sSpan); 
        switch FitType
            case 'Smoothed data'
                fitData = smoothData;
            case 'All data'
                fitData = sortedData;
            case 'Cancel'
                return
        end
        coefficients = polyfit(fitData(:,1), fitData(:,2),1);
        yfit = polyval(coefficients, fitData(:,1));
        bestFit = [fitData(:,1),yfit];
        % slope of line
        x1 = bestFit(1,1);
        x2 = bestFit(end,1);
        y1 = bestFit(1,2);
        y2 = bestFit(end,2);
        slope = (y2-y1)/(x2-x1);
        
        % save data for plotting:
        food(K).slope(end+1,1) = slope;
        food(K).X(end+1,:) = [x1,x2];
        food(K).Y(end+1,:) = [y1,y2];
    end
end

% --------- PLOT --------------
foodOpt = {'Plant', 'Yeast', 'Empty'};
sz = 50; width = 0.35;

fig = figure; set(fig, 'pos', [815 350 442 624])
hold on
% try to set up for just plant vs. yeast right now
N = length(food(1).slope);
% color code by sex
pdata = [food(1).slope,food(2).slope,food(3).slope];
for ii = 1:N
    kolor = pullSexColor(T.Sex{ii});
    p = plot([1,2,3],[pdata(ii,1),pdata(ii,2),pdata(ii,3)],'color', kolor, 'linewidth', 1.5);
    set(p, 'marker', 'd', 'markerfacecolor', kolor,'linestyle', ':');
end
% plot the avg line:
buff = 0.15;
for ii = 1:3 % for each food type:
    plot([ii-buff,ii+buff], [mean(pdata(:,ii)),mean(pdata(:,ii))],...
        'color', 'w', 'linewidth', 3);
end
xlim([0.75,3.25])
% ylim([-2,1.5])
hline(0,'w:')
set(gca,'xtick',[1,2,3],'xticklabels',foodOpt)
xlabel('Food')
ylabel('slope of temp-distance line')
title(ExpGroup)
fig = formatFig(fig, true);

save_figure(fig, [figDir ExpGroup ' temp distance across food'], '-png');

clearvars('-except',initial_vars{:})

%% FIGURE: effects of sex

sSpan = 120;
% determine the sex of flies in the experiment population:
sexList = T.Sex;
sexLab = unique(sexList); %sex labels
nsex = length(sexLab);

% ---- number of each sex group ----
for gg = 1:nsex
    loc = strcmpi(sexList,sexLab{gg});
    sex(gg).n = sum(loc);
    sex(gg).loc = loc;
    N(gg) = sex(gg).n;
end
fig = figure;
bar(N)
xlabel('Sex')
ylabel('Count')
set(gca, 'xtick', 1:nsex, 'xticklabels', sexLab)
fig = formatFig(fig, true);
save_figure(fig, [figDir ExpGroup ' fly sex histogram'], '-png');
% ----------------------------------

% Fit each trial individually and then plot them all...
numFitPoints = 20;
FitType = questdlg('Select data for linear fit:', '', 'Smoothed data', 'All data','Cancel', 'All data');
[threshHigh, threshLow] = getTempThresholds;

%build empty arena to stuff with data
G = struct('food',cell(1,3));
for sex = 1:nsex
  for K = 1:3 %foods
    G(sex).food(K).slope = [];
    G(sex).food(K).X = [];
    G(sex).food(K).Y = [];
  end
end
%Pull data for each sex group:
for trial = 1:ntrials
    sex = find(strcmp(sexList{trial},sexLab)); %sex group identity
    for well = 1:4
        smoothData = [];
        % pull data
        x = data(trial).occupancy.temp';
        try y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
        catch; y = (data(trial).occupancy.dist2wells(well).N(:,1))./pix2mm;
        end
        loc = (x<threshHigh&x>threshLow); % points within temp bounds
        x = x(loc);
        y = y(loc);
        raw = [x,y];

        % fit the data
        [~,idx] = sort(raw(:,1));
        sortedData = raw(idx,:);
        smoothData(:,1) = smooth(sortedData(:,1),sSpan);
        smoothData(:,2) = smooth(sortedData(:,2),sSpan); 
        switch FitType
            case 'Smoothed data'
                fitData = smoothData;
            case 'All data'
                fitData = sortedData;
            case 'Cancel'
                return
        end
        coefficients = polyfit(fitData(:,1), fitData(:,2),1);
        yfit = polyval(coefficients, fitData(:,1));
        bestFit = [fitData(:,1),yfit];
        % slope of line
        x1 = bestFit(1,1);
        x2 = bestFit(end,1);
        y1 = bestFit(1,2);
        y2 = bestFit(end,2);
        slope = (y2-y1)/(x2-x1);
        
        % save data for plotting:
        [~,K] = pullFoodColor(data(trial).wellLabels{well});
        G(sex).food(K).slope(end+1,1) = slope;
        G(sex).food(K).X(end+1,:) = [x1,x2];
        G(sex).food(K).Y(end+1,:) = [y1,y2];
    end
end
           
        
% --------- PLOT BY SEX --------------
foodOpt = {'Plant', 'Yeast', 'Empty'};
sz = 50; width = 0.35;

for K = 1:3
    fig = figure; set(fig, 'pos', [815 350 442 624])
    hold on
    for sex = 1:nsex
        % skip absent food groups
        a = G(sex).food.slope;
        if isempty(a)
            continue
        end
        kolor = pullFoodColor(foodOpt{K});
        plotdata = G(sex).food(K).slope;
        x = linspace(sex,sex+width,length(plotdata));
        scatter(x, plotdata, sz, kolor, 'filled')
        plot([sex-0.1,sex+width+0.1],[mean(plotdata),mean(plotdata)],'color',kolor,'linewidth',3)
    end
    xlim([0,nsex+1])
    hline(0,'w:')
    set(gca,'xtick',1:nsex,'xticklabels',sexLab,'XTickLabelRotation',45)
    xlabel('Sex')
    ylabel('slope of temp-distance line')
    title(foodOpt{K})
    fig = formatFig(fig, true);

    save_figure(fig, [figDir ExpGroup ' temp distance by sex for ' foodOpt{K}], '-png');
end
      
% --------- PLOT BY FOOD TYPE --------------
foodOpt = {'Plant', 'Yeast', 'Empty'};
sz = 50; width = 0.35;
for sex = 1:3 % only cycle through M,F,Mix
    fig = figure; set(fig, 'pos', [815 350 442 624])
    hold on
    % try to set up for just plant vs. yeast right now
    plotdata1 = G(sex).food(1).slope; %plant
    plotdata2 = G(sex).food(2).slope; %yeast
    N = length(plotdata1);
    kolors = Color('white', 'blue',N); %'powderblue', 'blue'
    % for fun:
    temp = plotdata2 - plotdata1;
    [~,idx] = sort(temp);
    pdata = [plotdata1(idx),plotdata2(idx)];
    F(sex).diff = pdata(:,1)-pdata(:,2);
    for ii = 1:N
        p = plot([1,2],[pdata(ii,1),pdata(ii,2)],'color', kolors(ii,:), 'linewidth', 2);
        set(p, 'marker', 'd', 'markerfacecolor', kolors(ii,:));
    end
    xlim([0.75,2.25])
    ylim([-2,1.5])
    hline(0,'w:')

    set(gca,'xtick',[1,2],'xticklabels',{'Plant', 'Yeast'})
    xlabel('Food')
    ylabel('slope of temp-distance line')
    title(sexLab{sex})
    fig = formatFig(fig, true);

    save_figure(fig, [figDir ExpGroup ' temp distance across food for ' sexLab{sex}], '-png');
end

clearvars('-except',initial_vars{:},'F','sexLab')
% ----- Plot relationship between plant and yeast preference ------
CList = {'deeppink', 'dodgerblue', 'white'};
fig = figure; hold on
for sex = 1:3
    pdata = F(sex).diff;
    h = histogram(pdata);
    set(h, 'facecolor',Color(CList{sex}),'facealpha', 0.75)
end
fig = formatFig(fig, true);
l = legend(sexLab{1:3});
set(l,'textcolor', 'w')
ylabel('Count'); xlabel('Diff in plant-yeast temperature slope')
save_figure(fig, [figDir ExpGroup ' temp slope comparison by sex'], '-png');
clearvars('-except',initial_vars{:})

%% FIGURE: Temperature hysteresis -- does temp ramp direction matter?

[threshHigh, threshLow] = getTempThresholds;

% %%-----plot the alignment of the increase|decreasing areas for a single trial------
% trial = 1;
% tPoints = getTempTurnPoints(T.TempProtocol{trial});
% fig = figure; hold on
% plot(data(trial).occupancy.temp,'color', 'w')
% y = rangeLine(fig);
% for ii = 1:tPoints.nUp
%     plot(tPoints.up(ii,:),[y,y],'color', Color('teal'),'linewidth',1)
% end
% v_line(reshape(tPoints.up,numel(tPoints.up),1),'w')
% for ii = 1:tPoints.nUp
%     plot(tPoints.down(ii,:),[y,y],'color', Color('orangered'),'linewidth',1)
% end
% % scatter(tPoints.UpROI,27*ones(1,length(tPoints.UpROI)),10,'g','filled')
% % scatter(tPoints.DownROI,27*ones(1,length(tPoints.DownROI)),10,'r','filled')
% fig = formatFig(fig, true);
% ylabel('temp (\circC)')
% xlabel('time (data points)')
% save_figure(fig,...
%     [figDir ExpGroup ' temp ramp delineation for ' T.TempProtocol{trial}], '-png');
% %-------------------------------------------------------------------------

% Pull the data from segments 
[plant, yeast, empty] = deal(struct('up', [], 'down',[]));
for trial = 1:ntrials
    tPoints = getTempTurnPoints(T.TempProtocol{trial}); %accomodates multiple temp protocols within the data group
    X = data(trial).occupancy.temp; % temperature
    for well = 1:4
        [~,K] = pullFoodColor(data(trial).wellLabels{well});
        % pull the dist2well data for each trial
        try Y = data(trial).occupancy.dist2wells(well).N(:,1)./pix2mm; %distance to wells
        catch; Y = data(trial).dist2wells(well).N(:,1)./pix2mm;
        end
        raw = [X',Y];
        % sort the data by food type
        switch K
            case 1
                plant.up = [plant.up; raw(tPoints.UpROI,:)];
                plant.down = [plant.down; raw(tPoints.DownROI,:)];
            case 2
                yeast.up = [yeast.up; raw(tPoints.UpROI,:)];
                yeast.down = [yeast.down; raw(tPoints.DownROI,:)];
            case 3
                empty.up = [empty.up; raw(tPoints.UpROI,:)];
                empty.down = [empty.down; raw(tPoints.DownROI,:)];
        end
    end
end

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

clearvars('-except',initial_vars{:})


%% PROCESS: auto generate increasing and decreasing temp ramp regions
% step 1: smooth temp ramp and label temp points as increasing or decreasing
X = data(1).occupancy.temp; % temperature
Y = data(1).occupancy.dist2wells(1).N(:,1); %distance to wells

temp = smooth(X,180); % this is the avg. temp for the previous 180 data points (except at start and end)
direct = diff(temp); % this tells me if the previous point was decreasing or increasing
increase = smooth(direct,360)>0;
decrease = smooth(direct,360)<0;


% 10 minute bins for determining temp history
dur = 10*60*3; %min * sec/min * frames/sec
a = smooth(increase,dur)>0.75; %increasing temp points
b = smooth(decrease,dur)>0.75;

upPoints = (a==~b);


figure; 
subplot(3,1,1) 
plot(temp)
title('temp')
subplot(3,1,2); hold on
plot(upPoints, 'color', 'g')
plot(a)
title('increasing')
subplot(3,1,3);
plot(b)
title('decreasing')



figure; 
subplot(3,1,1) 
plot(temp)
subplot(3,1,2)
plot(direct)
subplot(3,1,3);
plot(increase)


figure;
plot(smooth(increase,90))







%% vvv :SPECIALITY SECTIONS BELOW: vvvv
%% Plot female vs. grouped|male trials:

if strcmpi(ExpGroup,'cooling_warming_ramp')
    female_exps = [7,8,9,12];
elseif strcmpi(ExpGroup,'warming_cooling_ramp')
    female_exps = [2,4,9,11];
end

sSpan = 240;
[F_plant, F_yeast, F_empty, plant, yeast, empty] = deal([]);
for trial = 1:ntrials
    if ismember(trial,female_exps)
        for well = 1:4
            x = data(trial).occupancy.temp';
            y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            % plant food well:
            if strfind(data(trial).wellLabels{well},'Plant')
                F_plant = [F_plant; x, y];
            end
            % yeast food well:
            if strfind(data(trial).wellLabels{well},'Yeast')
                F_yeast = [F_yeast; x, y];
            end
            % empty food well:
            if strfind(data(trial).wellLabels{well},'Empty')
                F_empty = [F_empty; x, y];
            end 
        end
    else
        for well = 1:4
            x = data(trial).occupancy.temp';
            y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            % plant food well:
            if strfind(data(trial).wellLabels{well},'Plant')
                plant = [plant; x, y];
            end
            % yeast food well:
            if strfind(data(trial).wellLabels{well},'Yeast')
                yeast = [yeast; x, y];
            end
            % empty food well:
            if strfind(data(trial).wellLabels{well},'Empty')
                empty = [empty; x, y];
            end 
        end
    end
end 

numFitPoints = 1000;
threshHigh = 19.88;
threshLow = 8.02;
plotData = [];
% sort data
for ii = 1:6 %for each well type
    [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); 
    switch ii
        case 1 % plant food
            inputData = plant;
            kolor = pullFoodColor('Plant');
        case 2 % yeast food
            inputData = yeast;
            kolor = pullFoodColor('Yeast');
        case 3 % empty
            inputData = empty;
            kolor = pullFoodColor('Empty');
        case 4 % FM plant
            inputData = F_plant;
            kolor = Color('PaleGreen');
        case 5 % FM yeast
            inputData = F_yeast;
            kolor = Color('LemonChiffon');
        case 6 % FM empty
            inputData = F_empty;
            kolor = Color('lightgrey');
    end
    % cut off the high and low ends of data (to clean):
    loc = inputData(:,1)>threshHigh | inputData(:,1)<threshLow;
    inputData(loc,:) = [];
    % sort all the data by temperature:
    [~,idx] = sort(inputData(:,1));
    sortedData = inputData(idx,:);
    smoothData(:,1) = smooth(sortedData(:,1),sSpan);
    smoothData(:,2) = smooth(sortedData(:,2),sSpan);
    % find the line of best fit:
    coefficients = polyfit(smoothData(:,1), smoothData(:,2),1);
    xFit = linspace(smoothData(1,1), smoothData(end,1), numFitPoints);
    yFit = polyval(coefficients, xFit);
    yfit = polyval(coefficients, smoothData(:,1));
    yResid = smoothData(:,2)-yfit;
    SSresid = sum(yResid.^2);
    SStotal = (length(smoothData(:,2))-1) * var(smoothData(:,2));
    rsq = 1 - SSresid/SStotal; % get the R-square value
    % coefficient of correlation
    R = corrcoef(smoothData(:,1),smoothData(:,2));
    R = R(2,1);
    % save data for plotting:
    plotData(ii).smoothed = smoothData;
    plotData(ii).bestfit = [smoothData(:,1),yfit];
    plotData(ii).color = kolor;
    plotData(ii).rsq = round(rsq,2);
    plotData(ii).R = round(R,2);
end

% FIGURE: 
LW = 2;
fig = getfig; hold on
for ii = 1:6
    scatter(plotData(ii).smoothed(:,1), plotData(ii).smoothed(:,2), 30, plotData(ii).color, 'filled')
end
for ii = 1:6
%     plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'r', 'linewidth', LW+1)
    plot(plotData(ii).bestfit(:,1),plotData(ii).bestfit(:,2), 'color', 'w', 'linewidth', LW)
end
ylim([10,40])

% Labels:
xlabel('temperature (\circC)')
ylabel('distance from well (mm)')
title(strrep(ExpGroup,'_',' '))
formatFig(fig,true);

% label key:
str = {'Plant mixed', 'Yeast Mixed', 'Empty Mixed',...
       'Plant female', 'Yeast female', 'Empty female'};
for ii = 1:6
    legend_str{ii} = [str{ii}  ' R = ' num2str(plotData(ii).R) '  R^2 = ' num2str(plotData(1).rsq)];
end
l = legend(legend_str); 
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','Position', [0.7081 0.8031 0.1952 0.1606]);

save_figure(fig, [figDir ExpGroup ' female-male temp vs distance'], '-png');
clearvars('-except',initial_vars{:})

%% Plot male and female and mixed arenas for distance vs. temp

if strcmpi(ExpGroup,'cooling_warming_ramp')
    female_exps = [7,8];
    male_exps = [5,6];
elseif strcmpi(ExpGroup,'warming_cooling_ramp')
    female_exps = [2,4];
    male_exps = [1,3];
end

sSpan = 360;
[M_plant, M_yeast, M_empty, F_plant, F_yeast, F_empty, plant, yeast, empty] = deal([]);
for trial = 1:ntrials
    if ismember(trial,female_exps)
        for well = 1:4
            x = data(trial).occupancy.temp';
            y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            % plant food well:
            if strfind(data(trial).wellLabels{well},'Plant')
                F_plant = [F_plant; x, y];
            end
            % yeast food well:
            if strfind(data(trial).wellLabels{well},'Yeast')
                F_yeast = [F_yeast; x, y];
            end
            % empty food well:
            if strfind(data(trial).wellLabels{well},'Empty')
                F_empty = [F_empty; x, y];
            end 
        end
    elseif ismember(trial,male_exps)
        for well = 1:4
            x = data(trial).occupancy.temp';
            y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            % plant food well:
            if strfind(data(trial).wellLabels{well},'Plant')
                M_plant = [M_plant; x, y];
            end
            % yeast food well:
            if strfind(data(trial).wellLabels{well},'Yeast')
                M_yeast = [M_yeast; x, y];
            end
            % empty food well:
            if strfind(data(trial).wellLabels{well},'Empty')
                M_empty = [M_empty; x, y];
            end 
        end
    else
        for well = 1:4
            x = data(trial).occupancy.temp';
            y = (data(trial).dist2wells(well).N(:,1))./pix2mm; %convert the data to mm
            % plant food well:
            if strfind(data(trial).wellLabels{well},'Plant')
                plant = [plant; x, y];
            end
            % yeast food well:
            if strfind(data(trial).wellLabels{well},'Yeast')
                yeast = [yeast; x, y];
            end
            % empty food well:
            if strfind(data(trial).wellLabels{well},'Empty')
                empty = [empty; x, y];
            end 
        end
    end
end 

numFitPoints = 1000;
threshHigh = 19.88;
threshLow = 8.02;
plotData = [];
% sort data
for ii = 1:9 %for each well type
    [smoothData,coefficients,xFit,yFit,sortedData,yfit,] = deal([]); 
    switch ii
        case 1 % plant food
            inputData = plant;
            kolor = Color('MediumSeaGreen');
        case 2 % yeast food
            inputData = yeast;
            kolor = pullFoodColor('Yeast');
        case 3 % empty
            inputData = empty;
            kolor = Color('white');
        case 4 % F plant
            inputData = F_plant;
            kolor = Color('PaleGreen');
        case 5 % F yeast
            inputData = F_yeast;
            kolor = Color('LemonChiffon');
        case 6 % F empty
            inputData = F_empty;
            kolor = Color('lightgrey');
        case 7 % M plant
            inputData = M_plant;
            kolor = Color('DarkGreen');
        case 8 % M yeast
            inputData = M_yeast;
            kolor = Color('DarkOrange');
        case 9 % M empty
            inputData = M_empty;
            kolor = Color('Gray');
            
    end
    % cut off the high and low ends of data (to clean):
    loc = inputData(:,1)>threshHigh | inputData(:,1)<threshLow;
    inputData(loc,:) = [];
    % sort all the data by temperature:
    [~,idx] = sort(inputData(:,1));
    sortedData = inputData(idx,:);
    smoothData(:,1) = smooth(sortedData(:,1),sSpan);
    smoothData(:,2) = smooth(sortedData(:,2),sSpan);
    % find the line of best fit:
    coefficients = polyfit(smoothData(:,1), smoothData(:,2),1);
    xFit = linspace(smoothData(1,1), smoothData(end,1), numFitPoints);
    yFit = polyval(coefficients, xFit);
    yfit = polyval(coefficients, smoothData(:,1));
    yResid = smoothData(:,2)-yfit;
    SSresid = sum(yResid.^2);
    SStotal = (length(smoothData(:,2))-1) * var(smoothData(:,2));
    rsq = 1 - SSresid/SStotal; % get the R-square value
    % coefficient of correlation
    R = corrcoef(smoothData(:,1),smoothData(:,2));
    R = R(2,1);
    % save data for plotting:
    plotData(ii).smoothed = smoothData;
    plotData(ii).bestfit = [smoothData(:,1),yfit];
    plotData(ii).color = kolor;
    plotData(ii).rsq = round(rsq,2);
    plotData(ii).R = round(R,2);
end


% FIGURE: 
nrows = 1;
ncol = 6;
sb(1).idx = 1:3;
sb(2).idx = 4:5;
sb(3).idx = 6;
SZ = 120;
LW = 4; 
fig = getfig; 
% avg data
subplot(nrows, ncol, sb(1).idx)
hold on
for ii = 1:9
    scatter(plotData(ii).smoothed(:,1), plotData(ii).smoothed(:,2), 30, plotData(ii).color, 'filled')
end
ylim([10,40])
xlabel('temperature (\circC)')
ylabel('distance from well (mm)')
% best fit line
subplot(nrows, ncol, sb(2).idx)
hold on
for ii = 1:9
    slope(ii).x = plotData(ii).bestfit([1,end],1);
    slope(ii).y = plotData(ii).bestfit([1,end],2);
    slope(ii).m = (slope(ii).y(2)-slope(ii).y(1))/(slope(ii).x(2)-slope(ii).x(1));
    plot(slope(ii).x, slope(ii).y, 'color', plotData(ii).color, 'linewidth', LW)
end
ylim([10,40])
xlabel('temperature (\circC)')
% ylabel('distance from well (mm)')
title(strrep(ExpGroup,'_',' '))
% slope points
subplot(nrows, ncol, sb(3).idx)
hold on
for ii = 1:9
    slope(ii).m = (slope(ii).y(2)-slope(ii).y(1))/(slope(ii).x(2)-slope(ii).x(1));
    scatter(1, slope(ii).m, SZ, plotData(ii).color, 'filled')
end
xlim([0,2])
hline(0,'w:')
% xlabel('temperature (\circC)')
ylabel('slope of best fit line')
% ylim([-0.8,0.2])

formatFig(fig,true, [nrows,ncol], sb);
set(gca, 'XColor', 'k')

subplot(nrows, ncol, sb(1).idx)
%legend color code:
str = {'Plant', 'Yeast', 'Empty'};
% for ii = 1:9
%     legend_str{ii} = [str{ii}  ': R = ' num2str(plotData(ii).R) '  R^2 = ' num2str(plotData(ii).rsq)];
% end
l = legend(str); 
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','Position', [0.4278 0.8448 0.0662 0.0689]);

save_figure(fig, [figDir ExpGroup ' temp vs distance by sex'], '-png');
clearvars('-except',initial_vars{:})




%% SECTIONS UNDER EDIT:
%% Group occupancy across the trials % still in progress!!!

% for now: average across the experiments (assuming they are TEMP locked)
[A,B,C,D] = deal([]);
for trial = 1:ntrials
    A = [A,data(trial).occupancy.occ(:,1)]; %(:,well) % well occupancy probability for each well
    B = [B, data(trial).occupancy.occ(:,2)];
    C = [C, data(trial).occupancy.occ(:,3)];
    D = [D, data(trial).occupancy.occ(:,4)];
end
% average across the trials TODO update the timing according to temp alignment
for trial = 1:ntrials
    temperature(trial).log = data(trial).tempLog(:,2);
end
% TEMP ALIGNED TIME (ONLY WORKS FOR TEMP HOLD OR SAME GROUP TRIALS
TAT = data(1).time();
temperature = data(1).occupancy.temp;
wellLabels = data(1).wellLabels;

nrows = 4;
ncols = 1;
subs(1).idx = 1;
subs(2).idx = 2:4;

fig = getfig;
subplot(nrows, ncols, subs(1).idx)
    plot(TAT,temperature, 'linewidth', 2, 'color', 'w')
    ylim([5,30])
    ylabel('Temp (\circC)')
subplot(nrows, ncols, subs(2).idx)
    sSpan = 240;
    hold on
    Y = [smooth(mean(A,2),sSpan), smooth(mean(B,2),sSpan),...
         smooth(mean(C,2),sSpan), smooth(mean(D,2),sSpan)];
    h = area(TAT,Y);
    for well = 1:4
        h(well).FaceColor = pullFoodColor(data(1).wellLabels{well});
    end
    set(gca, 'tickdir', 'out')
    l = legend(strrep(wellLabels,'_','-'));
    set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k',...
    'position', [0.7457 0.6520 0.0883 0.0772]);% [0.8780 0.8119 0.0963 0.1126])%

ylabel('Occupancy')
xlabel('Time (min)')
formatFig(fig, true, [nrows, ncols], subs)
subplot(nrows, ncols, subs(1).idx)
set(gca, 'XColor', 'k')
title(strrep(ExpGroup,'_','-'),'color', 'w')

save_figure(fig, [figDir 'Quadrant Occupation'], '-png');

clearvars('-except',initial_vars{:})

%% Temperature and occupancy relationship: [takes a few seconds to run]
wellLabels = data(1).wellLabels;
% organize the occupation by fly count in each temperature
% --> temperature for each fly at each timepoint
tic
[sumA,sumB,sumC,sumD,tempAll] = deal(zeros(30,1));
for trial = 1:ntrials
    [A,B,C,D] = deal([]);
    temperature = round(data(trial).occupancy.temp);
    raw = data(trial).occupancy.count;
    for ii = 1:size(raw,1) %for each data point
        A = [A; temperature(ii)*ones(raw(ii,1),1)];
        B = [B; temperature(ii)*ones(raw(ii,2),1)];
        C = [C; temperature(ii)*ones(raw(ii,3),1)];
        D = [D; temperature(ii)*ones(raw(ii,4),1)];
    end
    for temp = 5:30
        sumA(temp) = sumA(temp) + sum(A==temp);
        sumB(temp) = sumB(temp) + sum(B==temp);
        sumC(temp) = sumC(temp) + sum(C==temp);
        sumD(temp) = sumD(temp) + sum(D==temp);
        tempAll(temp) = tempAll(temp) + sum(temperature==temp);
    end
end
toc

% normalize the number of flies for each temperature but the total time
% spent at a given temperature:
tempFrac = tempAll./(sum(tempAll));
ratio = 1; %each point will represent 20 flies normalized for total time spent at each temperature
arenaA = round(sumA.*tempFrac./ratio);
arenaB = round(sumB.*tempFrac./ratio);
arenaC = round(sumC.*tempFrac./ratio);
arenaD = round(sumD.*tempFrac./ratio);


% increase points for each temperature

[yA,yB,yC,yD] = deal([]);
for temp = 1:length(arenaA)
    yA = [yA; temp*ones(arenaA(temp),1)];
    yB = [yB; temp*ones(arenaB(temp),1)];
    yC = [yC; temp*ones(arenaC(temp),1)];
    yD = [yD; temp*ones(arenaD(temp),1)];
end
xA = 1*ones(length(yA),1);
xB = 2*ones(length(yB),1);
xC = 3*ones(length(yC),1);
xD = 4*ones(length(yD),1);

% Plot the occupation by temperature plots
sz = 30;
fig = getfig;
    hold on
    swarmchart(xA,yA,sz, pullFoodColor(wellLabels{1}), 'filled');
    swarmchart(xB,yB,sz, pullFoodColor(wellLabels{2}), 'filled');
    swarmchart(xC,yC,sz, pullFoodColor(wellLabels{3}), 'filled');
    swarmchart(xD,yD,sz, pullFoodColor(wellLabels{4}), 'filled');

% swarmchart(x,y,sz,c)














%% Collapse Well ROIs across all trials for each temperature to create 'occupancyFrames'
gridSize = [50,50];
nvids = length(data(1).temp);

% loop through all trials
for trial = 1:ntrials
  temp = data(trial).temp; 
  for ii = 1:4 % well#
    % crop each of the frames from the video into well rois
    frames = data(trial).frame;
    mask = data(trial).quadMask(ii).mask;
    props = regionprops(~mask, 'BoundingBox');
    mask = repmat(mask,[1,1,size(frames,3)]); %mask for all frames
    frames(mask) = 0; %mask out other data outside ROI
    for vid = 1:nvids
        img = imcrop(frames(:,:,vid), props.BoundingBox); %crop to the circle
        maskedImage(:,:,vid) = img;
        roiImage(:,:,vid) = imresize(img, gridSize); %bin across pixels
    end
    well(ii).img = roiImage;
    well(ii).maskedImg = maskedImage;
    well(ii).temp = temp;
  end
  data(trial).well = well;
end

% align across videos roughly (TODO: refine this later)
Wells = [];
for roi = 1:4
    Wells(roi).name = data(1).trialData(1).params.(['well_' num2str(roi)]);
    Wells(roi).img = zeros(size(data(1).well(roi).img));
    for trial = 1:ntrials
        Wells(roi).img = Wells(roi).img + data(trial).well(roi).img;
        Wells(roi).occ(:,trial) = data(trial).y(:,roi); % determine this...
    end
    % find density range for the rois:
    M(roi) = max(Wells(roi).img,[],'all');
end

% Normalize the well density figure:
gridS = 10;
n = gridS+2;
well = [];
for roi = 1:4
    well(roi).img = Wells(roi).img / max(M);
    dummy = nan(n,n,nvids);
    img = imresize(well(roi).img, [gridS,gridS]);
    dummy(2:gridS+1,2:gridS+1,:) = img;
    well(roi).bin = dummy;
end
% grouped image:
occupancyFrames = [well(1).bin, well(2).bin;...
                   well(4).bin, well(3).bin];
               
occupancyNames = {Wells(1).name, Wells(2).name;...
                   Wells(4).name, Wells(3).name};              

% quick vid preview of any structure over time
fig = figure; set(fig, 'color', 'k');
disp(occupancyNames)
for vid = 1:nvids
    img = occupancyFrames(:,:,vid);
    imagesc(img)
    pause(0.1)
end
     






%% IDEAS for later:




























