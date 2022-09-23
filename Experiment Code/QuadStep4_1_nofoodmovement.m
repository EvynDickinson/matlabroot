

%% LOAD:multiple trials that are grouped into a structure
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
    figDir = [baseFolder 'Data structures/' ExpGroup '/'];
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
        dirc = [baseFolder, trialDate, '/Arena ' trialArena '/analysis/' trialExpID trialArena  ' timecourse data.mat'];
        
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
            catch
                try data(trial).(var2load{ii}) = temp.data.(var2load{ii});
                catch
                    if strcmp(var2load{ii},'dist2wells')
                        try data(trial).(var2load{ii}) = temp.data.dist2well;
                        catch; data(trial).(var2load{ii}) = [];
                        end
                    else
                        data(trial).(var2load{ii}) = [];
                    end
                end
            end
        end
        disp([expID{trial} arenas{trial}])
    end
end

%Pull and reorganize data within the structure:

% pix2mm = 12.8; %conversion from pixels to mm for these videos
initial_vars = {'ExpGroup','baseFolder', 'T', 'data', 'figDir', 'filePath',...
                'initial_vars', 'folder', 'ntrials', 'pix2mm'};
clearvars('-except',initial_vars{:})
if questdlg('Save loaded data?')
    save([figDir ExpGroup ' raw'],'-v7.3')
end
fprintf('Data loaded\n')

%% ANALYSIS: General data organization (forward and backwards compatability
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

%% FIGURE: visual check of temperature alignment across the experiments:

fig = figure; set(fig, 'position', [107 595 1150 235]); hold on
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

save_figure(fig, [figDir 'temperature alignment'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')

%% FIGURE: movement across time during the experiment

r = 3;
c = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;

fig = figure; set(fig,'pos',[246 253 1001 634])

subplot(r,c,sb(1).idx)
x = data(1).occupancy.time;
y = data(1).occupancy.temp;
plot(x,y,'linewidth',1.5,'color','w')
ylabel('Temp (\circ)')

subplot(r,c,sb(2).idx)
hold on
[avgSpeed,avgTime] = deal([]);
for trial = 1:ntrials
    x = data(trial).occupancy.time;
    y = data(trial).speed.avg;
    plot(x,y,'linewidth',1.5)

    % pull avg. speed
    avgTime = autoCat(avgTime,x,false);
    avgSpeed = autoCat(avgSpeed,y,false);
end
% Plot AVG speed
plot(mean(avgTime,2,'omitnan'),mean(avgSpeed,2,'omitnan'),'color','w','linewidth',2)


%labels
xlabel('time (min)')
ylabel('speed (mm/s)')

formatFig(fig,true,[3,1],sb);

save_figure(fig, [figDir 'unsmoothed speed across trials'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')


%% Avg movement or speed vs temp

inputVar =  questdlg('Which data type to compare?','','speed','movement','clustering','speed');
switch inputVar
    case 'movement'
        ylab = 'movement (~mm/s)';
        L_loc = 'southeast';
    case 'clustering'
        ylab = 'Inter-fly-distance (mm)';
        L_loc = 'northwest';
    case 'speed'
        ylab = 'speed (mm/s)';
        L_loc = 'northwest';
    case ''
        return
end

% allocate empty data structure & set params
food = struct;
food.N = [];
for trial = 1:ntrials
    % Screen out data pre/post data 
    tempPoints = getTempTurnPoints(T.TempProtocol{trial});
    roi = sort([tempPoints.DownROI,tempPoints.UpROI,tempPoints.HoldROI]);
    switch inputVar
        case 'movement'
            x = data(trial).occupancy.temp(1:end-1);
            y = data(trial).occupancy.movement;
        case 'clustering'
            y = data(trial).occupancy.IFD';
            x = data(trial).occupancy.temp;
        case 'speed'
            y = data(trial).speed.avg;
            x = data(trial).occupancy.temp;
    end 
    food.N = [food.N; x(roi),y(roi)];
end 

% Temperature range and bin size formatting
[threshHigh, threshLow] = getTempThresholds(T.TempProtocol);
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
x = t_roi(1:end-1);
y = food.avg(1:end-1);
y_err = food.err(1:end-1);
fill_data = error_fill(x, y, y_err);
h = fill(fill_data.X, fill_data.Y, kolor, 'EdgeColor','none');
set(h, 'facealpha', 0.2)
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



%% Could look at distance as a series of boot-strapped selections?

%%




























