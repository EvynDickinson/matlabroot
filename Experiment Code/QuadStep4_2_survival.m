% figure;
% exp = 1;
% trial = 1;
% y = ones(size(sROImask(exp).all(trial).m));
% y = y.*[1:size(y,2)];
% y(~sROImask(exp).all(trial).m) = nan;
% 
% plot(y)

%    GroupDataGUI_v2

%% Load data

clear; clc; warning off
blkbnd = true;
fig_type = '-png';

% Select data structure from folder names out of GroupDataGUI:
% Need both the Data Folder path (pulling data from here) and the
% structures folder (saving data to here)
baseFolder = getDataPath(5,0);
pathNames = getPathNames;
ExpGroup = selectFolder([baseFolder pathNames.grouped_trials],false,'Select data structure');
ExpGroup = ExpGroup{:};
% data structures folder:
figDir = [baseFolder pathNames.grouped_trials ExpGroup '/'];
    
rawDataFolder = [baseFolder, pathNames.single_trial];

% Does the data grouping already exist?
raw_file = [figDir ExpGroup ' raw.mat'];
if exist(raw_file,'file') == 2 %file exists
    oldList =  load(raw_file, 'T');
    newList = load([figDir 'fileList.mat'],'T');
end

% load the directory list:
load([figDir 'fileList.mat'])

% extract information
ntrials = size(T,1);
dates = T.Date;
arenas = T.Arena;
expID = T.ExperimentID;

% load data
n = cell(1,ntrials);
data = struct('dist2wells', n, 'wellLabels', n, 'occupancy', n,'nflies',n,'data', n,'plate', n, 'pix2mm', n);  
var2load = fieldnames(data);

for trial = 1:ntrials
    trial_ID = [dates{trial} '_' expID{trial} '_' arenas{trial}];
    filePath = [rawDataFolder trial_ID '/'];
    % filePath = [baseFolder, dates{trial}, '/Arena ' arenas{trial} '/analysis/'];
    if ~isfolder(filePath)
        warndlg(['File not found ' trial_ID])
    end
    
    temp = load([filePath expID{trial} arenas{trial} ' timecourse data v2.mat'], var2load{:});
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

    % load speed data:
    temp = load([filePath expID{trial} ' speed data v2.mat']);
    data(trial).speed = temp.speed;
    disp([expID{trial} arenas{trial}])
end

initial_vars = {'ExpGroup','baseFolder', 'T', 'data', 'figDir', 'filePath',...
        'initial_vars', 'folder', 'ntrials', 'FPS'};
clearvars('-except',initial_vars{:})

save([figDir ExpGroup ' raw.mat'],'-v7.3')

fprintf('Data loaded\n')
disp('next section')

%% 

sSpan = 180; % 3fps = 1min bin

speed = [];
a = [];
for trial = 1:ntrials 
    speed = autoCat(speed, data(trial).speed.avg, false);
    a = [a, length(data(trial).occupancy.time)];
end 
[b,c] = min(a);
avgspeed = mean(speed(1:b,:),2,'omitnan');

fig = getfig;
yy = smooth(avgspeed,sSpan,'moving');
plot(data(c).occupancy.time, (yy))


%%


% Smooth raw speed data for clearer visualization
binsize = 3;
for trial = 1:ntrials    
    numtracks = size(data(trial).speed.raw,2);
    numframes = length(data(trial).speed.raw);
    a = [];
    for track = 1:numtracks
        a(:,track) = smooth(data(trial).speed.raw(:,track),binsize,'moving');
    end
    data(trial).speed.smoothed_raw = a;
end

min_speed = 0.1; %min speed, up for debate, matched to courtship min
for trial = 1 %:ntrials
    for track = 1 %:numtracks
        still = data(trial).speed.smoothed_raw(:,track) < min_speed;
        loc = find(still);
    end
end

%have location in raw data for one trial and one track where speed is < 0.1
%need place to store locs 
%what percentage of tracks are still over time
%compare that to percentage of flies still over time























%% Load MANUAL data
clear; clc;
% baseFolder = getDataPath(2,2);

blkbnd = true;
fig_type = '-png';

% Load excel file
[excelfile, Excel, xlFile] = load_SurvivalQuadCounts;

% Find the experiments that have been counted
loc = cellfun(@isnan,excelfile(2:end,Excel.counted));
loc = ~loc;
rownums = find(loc)+1; 
eligible_files = excelfile([false;loc],[Excel.date, Excel.protocol, Excel.arena, Excel.counted]);
% FileNames = format_eligible_files(eligible_files);
temptype = unique(eligible_files(:,2));

% Select files 
fileIdx = listdlg('ListString', temptype,'ListSize',[350,450],'promptstring', 'Select temp protocol to process');
if isempty(fileIdx)
    disp('No trials selected')
    return
end

% Path data structures folder
baseFolder = getDataPath(5,2);
pathNames = getPathNames;
% ExpGroup = selectFolder([baseFolder pathNames.grouped_trials],false,'Select data structure');
% ExpGroup = ExpGroup{:};

% If survival curves folder isn't made yet, make it
ExpGroup = 'Berlin survival comparisons no food';
rootDir = createFolder([baseFolder pathNames.group_comparision ExpGroup '/']);
figDir = createFolder([rootDir 'Figures/']);

nexps = length(fileIdx);


%% Pull the protocol name and trials under the selected temp protocol
data = struct;

for exp = 1:nexps
tempprotocol = temptype{fileIdx(exp)};
temptrials = strcmp(tempprotocol, excelfile(:,5)); % binary
data(exp).tempprotocol = tempprotocol;

% Which columns in excel file to plot based on each temp protocol
switch tempprotocol
    case 'survival_curve_40C'
        data(exp).columns = [11,13:52];
        data(exp).color =  Color('red');
        data(exp).temp = 40;
    case 'survival_curve_35C'
        data(exp).columns = [11,13:52];
        data(exp).color =  Color('tomato');
        data(exp).temp = 35;
    case 'survival_curve_10C'
        data(exp).columns = [11:14,16,18,20,22,24,27,29,32,34,37,39,42,44,47,49,52,54,57,59];
        data(exp).color = Color('dodgerblue');
        data(exp).temp = 10;
    case 'survival_curve_5C'
        data(exp).columns = [11:14,16,18,20,22,24,27,29,32,34,37,39,42,44,47,49,52,54,57,59];
        data(exp).color = Color('blue');
        data(exp).temp = 5;
end

data(exp).raw = cell2mat([excelfile(temptrials,data(exp).columns)]);
data(exp).avg = mean(data(exp).raw,'omitnan');
data(exp).timestmp = cell2mat([excelfile(1,data(exp).columns)]);
data(exp).time = data(exp).timestmp*10;
data(exp).ntrials = size(data(exp).raw,1);
data(exp).perc = (data(exp).raw ./ cell2mat([excelfile(temptrials,10)]))*100;
data(exp).avgperc = mean(data(exp).perc,'omitnan');
end


%% Plot number of incapacitated/dead flies overtime
fig = getfig;
for exp = 1:nexps
    hold on
    y_err = (std(data(exp).raw,0,1,'omitnan'))./sqrt(data(exp).ntrials); % calculate SEM
    plot_error_fills(1, data(exp).time, data(exp).avg, y_err, data(exp).color, fig_type, 0.35); % plot error
    plot(data(exp).time,data(exp).avg,'linewidth',3,'color',data(exp).color)
end
formatFig(fig, blkbnd);
    ylim([0 20])
    set(gca,"YTick", 0:3:20, 'FontSize', 20)
    xlabel('Time (min)')
    ylabel('Number of incapaciated flies')
save_figure(fig,[figDir, '/' ExpGroup ' death overtime'], fig_type);

%% Percentage flies incapacitated/dead overtime

fig = getfig;
for exp = 1:nexps
    hold on
    y_err = (std(data(exp).perc,0,1,'omitnan'))./sqrt(data(exp).ntrials); % calculate SEM
    plot_error_fills(1, data(exp).time/60, data(exp).avgperc, y_err, data(exp).color, fig_type, 0.35); % plot error
    plot((data(exp).time)/60,data(exp).avgperc,'linewidth',3,'color',data(exp).color)
end

formatFig(fig, blkbnd);
    ylim([0 100])
    set(gca,"YTick", 0:25:100, 'XTick', 0:6:48, 'FontSize', 15)
    xlabel('Time (hours)')
    ylabel('Percentage of incapaciated flies')

save_figure(fig,[figDir, '/' ExpGroup ' death overtime'], fig_type);

%% Plot number of alive flies overtime
% numflies = cell2mat([excelfile(temptrials,10)]);
% aliveflies = numflies - cell2mat([excelfile(temptrials,11:51)]);
% % perc_alive = (aliveflies/numflies)*100;
% 
% fig = getfig;
%     plot(time,aliveflies,'linewidth',2)
% formatFig(fig, blkbnd);
%     ylim([0 20])
%     % xlim([0 500])
%     xlabel('Time (min)')
%     ylabel('Number of alive flies')


