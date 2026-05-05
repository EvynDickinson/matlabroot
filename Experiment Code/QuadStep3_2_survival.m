
%    GroupDataGUI_v2
clearvars('-except',initial_vars{:})

%% LOAD: Individual survival data structure

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

% Pull temperature protocol info
tPoints = getTempTurnPoints(T.TempProtocol{1});

sSpan = 180; % 3fps = 1min bin
foreColor = formattingColors(blkbnd); % get background colors
initial_vars = {'ExpGroup','baseFolder', 'T', 'data', 'tPoints','figDir', 'filePath',...
        'initial_vars', 'folder', 'ntrials', 'FPS','sSpan','blkbnd','fig_type',...
        'manual','manual_excelfile', 'manual_Excel', 'manual_xlFile',...
        'avgpercstill','mexps','foreColor','mtrials','maxspeed','movingAD','neverDies'};
clearvars('-except',initial_vars{:})

save([figDir ExpGroup ' raw.mat'],'-v7.3')

fprintf('Data loaded\n')
disp('next section')

%% ANALYSIS: Smooth speed data
clearvars('-except',initial_vars{:})

% Smooth raw speed data for clearer visualization
binsize = 3;
for trial = 1:ntrials    
    numtracks = size(data(trial).speed.raw,2);
    a = [];
    for track = 1:numtracks
        a(:,track) = smooth(data(trial).speed.raw(:,track),binsize,'moving');
    end
    data(trial).speed.smoothed_raw = a;
end

%% CHECK: Max speeds for raw and smoothed speed?

maxspeed = [];
for trial = 1:ntrials
    b = max(data(trial).speed.raw);
    c = max(b);
    maxspeed = [maxspeed, c];   
end
maxspeed = maxspeed';

maxsmoothed = [];
for trial = 1:ntrials
    b = max(data(trial).speed.smoothed_raw);
    c = max(b);
    maxsmoothed = [maxsmoothed, c];    
end
maxsmoothed = maxsmoothed';

%% MANUAL LOAD: Manual excel data for ground truthing
% baseFolder = getDataPath(2,2);
% clearvars('-except',initial_vars{:})


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


% % FOR MULTIPLE TEMP COMPARISON
% % Path data structures folder
% baseFolder = getDataPath(5,2);
% pathNames = getPathNames;
% % ExpGroup = selectFolder([baseFolder pathNames.grouped_trials],false,'Select data structure');
% % ExpGroup = ExpGroup{:};
% 
% % If survival curves folder isn't made yet, make it
% ExpGroup = 'Berlin survival comparisons no food';
% rootDir = createFolder([baseFolder pathNames.group_comparision ExpGroup '/']);
% figDir = createFolder([rootDir 'Figures/']);
% 
% 

%% MANUAL ANALYSIS: Pull the protocol name and trials under the selected temp protocol for ground truthing

manual = struct;
mexps = length(fileIdx);

for exp = 1:mexps
tempprotocol = temptype{fileIdx(exp)};
temptrials = strcmp(tempprotocol, excelfile(:,5)); % binary
manual(exp).tempprotocol = tempprotocol;

% Which columns in excel file to plot based on each temp protocol
switch tempprotocol
    case 'survival_curve_40C'
        manual(exp).columns = [11,13:32];
        manual(exp).color =  Color('red');
        manual(exp).temp = 40;
    case 'survival_curve_35C'
        manual(exp).columns = [11,13:32];
        manual(exp).color =  Color('tomato');
        manual(exp).temp = 35;
    case 'survival_curve_10C'
        manual(exp).columns = [11:14,16,18,20,22,24,27,29,32,34,37,39,42,44,47,49,52,54,57,59];
        manual(exp).color = Color('dodgerblue');
        manual(exp).temp = 10;
    case 'survival_curve_5C'
        manual(exp).columns = [11:14,16,18,20,22,24,27,29,32,34,37,39,42,44,47,49,52,54,57,59];
        manual(exp).color = Color('blue');
        manual(exp).temp = 5;
end

manual(exp).raw = cell2mat([excelfile(temptrials,manual(exp).columns)]);
manual(exp).avg = mean(manual(exp).raw,'omitnan');
manual(exp).timestmp = cell2mat([excelfile(1,manual(exp).columns)]);
manual(exp).time = manual(exp).timestmp*10; % time in minutes
manual(exp).ntrials = size(manual(exp).raw,1);
manual(exp).perc = (manual(exp).raw ./ cell2mat([excelfile(temptrials,10)]))*100;
manual(exp).avgperc = mean(manual(exp).perc,'omitnan');
manual(exp).flycount = cell2mat([excelfile(temptrials,10)]);
end

manual_excelfile = excelfile;
manual_Excel = Excel;
manual_xlFile = xlFile;
clearvars('excelfile', 'Excel', 'xlFile')

clearvars('-except',initial_vars{:})

%% MANUAL ANALYSIS: Find when flies died in manual count

maxflies = manual.flycount;
match = (manual.raw==maxflies);
if ~any(match(:,end) == 1) % if no flies are dead by the end of the manual exp
    neverDies = true;
else
    c = [];
    dead = [];
    for trial = 1:size(match,1)
        [~,c(trial)] = find(match(trial,:),1);
        dead = manual.timestmp(c);
    end
    dead = dead'; % video where all flies are counted dead
    manual.TOD = dead*10; % time of death
    manual.TOD_min = min(manual.TOD);
    manual.TOD_max = max(manual.TOD);    
    manual.TOD_mode = mode(manual.TOD);
end

clearvars('-except',initial_vars{:})

%% ANALYSIS: Calculate max average speed recorded in auto after flies die in manual

% Max speed in auto data after TOD in manual data
if neverDies
    disp('No total TOD to calculate (there were still flies alive at end of exp)')
else
    maxspeed = NaN(3,ntrials);
    for trial = 1:ntrials
        % Quickest death
        time = data(trial).data.T.time;
        frames = time >= manual.TOD_min;
        % Save max speed after quickest TOD for each trial into a matrix
        maxspeed(1,trial) = max(data(trial).speed.avg(frames));
       
        % Most common death
        frames = time >= manual.TOD_mode;
        % Save max speed after most common TOD for each trial into a matrix
        maxspeed(2,trial) = max(data(trial).speed.avg(frames));
        
        % Slowest death
        frames = time >= manual.TOD_max;
        % Save max speed after slowest TOD for each trial into a matrix
        maxspeed(3,trial) = max(data(trial).speed.avg(frames));
    end
    
    % for ii = 1:3
    %     mean(maxspeed(ii,:))
    % end
    
    movingAD = NaN(3,ntrials);
    for trial = 1:ntrials
    % Quickest death
    frames = time >= manual.TOD_min;
    % Save avg speed after quickest TOD for each trial into a matrix
    movingAD(1,trial) = mean(data(trial).speed.avg(frames),'omitnan');
    
    % Most commong death
    frames = time >= manual.TOD_mode;
    % Save avg speed after quickest TOD for each trial into a matrix
    movingAD(2,trial) = mean(data(trial).speed.avg(frames),'omitnan');
    
    % Slowest death
    frames = time >= manual.TOD_max;
    % Save avg speed after quickest TOD for each trial into a matrix
    movingAD(3,trial) = mean(data(trial).speed.avg(frames),'omitnan');
    end
    
    % for ii = 1:3
    %     mean(movingAD(ii,:))
    % end
end

%% ANALYSIS: Calculate proportion of still flies
% Find proportion of still flies throughout experiment
clearvars('-except',initial_vars{:})

min_speed = 1; %min speed, based on speed distribution

for trial = 1:ntrials
    numtracks = size(data(trial).speed.raw,2);
    % Create matrix for when fly speed is < min allocated 
    still = [];
    loc = [];
    for track = 1:numtracks
        % Create matrix of logicals for whether speed is > or < min speed
        b = data(trial).speed.raw(:,track) < min_speed;
        still = [still, b];
        % Find locations where speed is < min speed
        c = find(b);
        loc = autoCat(loc, c, false); 
    end
    % Calcuate proportion of still tracks from total number of tracks for each frame
    incap = sum(still,2);
    d = ~isnan(data(trial).speed.raw);
    data(trial).speed.sizetracks = sum(d,2);
    data(trial).speed.still = incap./data(trial).speed.sizetracks;     
end

% % Find places in Still where values are infinite
% infinity_loc = [];
% for trial = 1:ntrials
%     infinity_loc = autoCat(infinity_loc,find(isinf(data(trial).speed.still)),false);
%     loc2 = length(infinity_loc);
%     for i = 1:loc2
%         data(trial).speed.still(infinity_loc(i)) = 0;
%     end
% end

clearvars('-except',initial_vars{:})

%% FIGURE: Histogram of fly speed 

% FIGURE
fig = getfig;
hold on
for trial = 1:ntrials
    histogram(data(trial).speed.smoothed_raw)%,"BinWidth",0.01)
end
xlim([0 2])
for trial = 1:ntrials
    dummy = max(data(trial).speed.smoothed_raw);
end

% Format figure
formatFig(fig,blkbnd);
% set(gca, 'yscale', 'log') % TODO: figure the scaling out on this
xlabel('Speed (mm/s)')
ylabel('Frequency')
title(trial)

% Save figure
save_figure(fig,[figDir, 'Speed histogram'], fig_type);

% Histogram of smoothed raw speed after each TOD (only if all flies die by end)
if ~neverDies   
    deathList = [manual.TOD_min,manual.TOD_max,manual.TOD_mode];
    for ii = 1%:3
        % Pull frames after each TOD
        frames = time >= deathList(ii);
        % FIGURE
        fig = getfig;
        hold on
        for trial = 1:ntrials
            histogram(data(trial).speed.smoothed_raw(frames,:),"BinWidth",0.01)
        end
        xlim([0 2])
        % for trial = 1:ntrials
        %     dummy = max(data(trial).speed.smoothed_raw);
        % end
    
        % Format figure
        formatFig(fig,blkbnd);
        xlabel('Speed (mm/s)')
        ylabel('Frequency')
        title(deathList(ii),'Color',foreColor)
        axis square
        
        % Save figure
        save_figure(fig,[figDir, 'Speed histogram after quickest TOD'], fig_type);
    end
end

%% FIGURE: Smoothed avg speed over time

showTemp = false;

if showTemp
    r = 3;
    c = 1;
    sb(1).idx = 1; % temp protocol
    sb(2).idx = 2:3; % speed
end

% FIGURE
fig = getfig;
if showTemp
    % Plot temp protocol
    subplot(r,c,sb(1).idx)
        plot(data(1).occupancy.time, data(1).data.T.temperature,'LineWidth',1,'Color',foreColor)
% Plot average speed for each trial
subplot(r,c,sb(2).idx)
end
    hold on
    for trial = 1:ntrials        
        y = data(trial).speed.avg;
        % Smooth data
        yy = smooth(y,sSpan,'moving');
        plot(data(trial).occupancy.time, yy,'LineWidth',2)        
    end
    % Plot vertical line at each manual TOD's
    deathList = [manual.TOD_min,manual.TOD_max,manual.TOD_mode];
    for ii = 1:3
        xline(deathList(ii),'Color','r','LineWidth',1.5)
    end
    xlabel('time (min)')
    ylabel('Average speed (mm/s)')
    xlim([0 400])

% Format fig
if showTemp
    formatFig(fig,blkbnd,[r,c],sb)
else
    formatFig(fig,blkbnd)
end

if showTemp
    % Remove x axis for top subplot
    subplot(r,c,sb(1).idx)
    set(gca, 'xcolor', 'none')
    xlim([0 400])
end
title(ExpGroup,'Color',foreColor)

% Save figure
save_figure(fig,[figDir, 'Avg speed over time'], fig_type);

%% FIGURE: Smoothed raw speed over time - trial 4

% FIGURE
fig = getfig;
hold on
% Plot raw speed 
for track = 1:15
    y = data(4).speed.raw(:,track);
    % Smooth data
    yy = smooth(y,sSpan,'moving');
    plot(data(trial).occupancy.time, yy,'LineWidth',2)        
end
% Plot vertical line at each manual TOD's
deathList = [manual.TOD_min,manual.TOD_max,manual.TOD_mode];
for ii = 1:3
    xline(deathList(ii),'Color','r','LineWidth',1.5)
end

% Format fig
formatFig(fig,blkbnd)
% xlim([0 400])
xlabel('time (min)')
ylabel('Raw speed (mm/s)')
axis square
title([ExpGroup, ' trial 4'],'Color',foreColor)

% Save figure
save_figure(fig,[figDir, 'Raw speed over time - trial 4 only'], fig_type);

%% FIGURE: Stillness over time
clearvars('-except',initial_vars{:})

r = 4;
c = 3;
sb(1).idx = 1:2; % temp protocol
sb(2).idx = [4,5,7,8,10,11]; % stillness
sb(3).idx = [6,9,12];

fig = getfig;
 % Plot temp protocol
subplot(r,c,sb(1).idx)
    plot(data(1).occupancy.time, data(1).data.T.temperature,'LineWidth',1,'Color',foreColor)
    ylabel('Temperature (\circC)')
% Plot stillness for each trial
subplot(r,c,sb(2).idx)
    hold on
    for trial = 1:ntrials     
        notmoving = data(trial).speed.still;
        y = smooth(notmoving,sSpan,'moving');
        yy = y*100;
        plot(data(trial).occupancy.time, yy,'LineWidth',2)    
    end
    xlabel('Time (min)')
    ylabel('Percentage of flies still')
    ylim([0 100])
% Plot scatterplot of percentage of flies moving after recovery
subplot(r,c,sb(3).idx)
    hold on
    recovery = (tPoints.hold(3):tPoints.hold(4))';
    for trial = 1:ntrials
        yy = mean(data(trial).speed.avg(recovery),'omitnan');
        scatter(1,yy,50,'filled', 'XJitter','density','XJitterWidth',1);
    end
    ylim([0 25])
    xticks(1)
    xticklabels([])
    ylabel('Percentage of flies moving during recovery')



% Format figure
formatFig(fig,blkbnd,[r,c],sb);
subplot(r,c,sb(1).idx)
set(gca, 'xcolor', 'none')

% Save figure
save_figure(fig,[figDir, 'Percentage flies still overtime timecourse with recovery'], fig_type);

%% FIGURE: Average percentage of still flies over time (all trials combined)
clearvars('-except',initial_vars{:})

% Create matrix with all stillness values for each trial
incap = [];
a = [];
for trial = 1:ntrials 
    incap = autoCat(incap, data(trial).speed.still, false);
    a = [a, length(data(trial).occupancy.time)]; % frame total for each exp
end 
[b,e] = min(a); %[minimun, where column minimum value is located within a]
% Calculate average percentage still 
avgstill = mean(incap(1:b,:),2,'omitnan');
avgpercstill = avgstill*100;

r = 4;
c = 3;
sb(1).idx = 1:2; % temp protocol
sb(2).idx = [4,5,7,8,10,11]; % stillness
sb(3).idx = [6,9,12];

% FIGURE
fig = getfig;
 % Plot temp protocol
subplot(r,c,sb(1).idx)
    plot(data(e).occupancy.time, data(e).data.T.temperature,'LineWidth',1,'Color',foreColor)
    ylabel('Temperature (\circC)')
% Plot stillness for each trial
subplot(r,c,sb(2).idx)
hold on
    y = smooth(avgpercstill,sSpan,'moving');
    % y = smooth(y,sSpan,'moving'); % for extra smoothing
    plot(data(e).occupancy.time, y,'linewidth',1)
    ylabel('Avgerage percentage of flies still')
    xlabel('Time (min)')
    ylim([0 100])
% Plot scatterplot of percentage of flies moving after recovery
subplot(r,c,sb(3).idx)
    hold on
    recovery = (tPoints.hold(3):tPoints.hold(4))';
    for trial = 1:ntrials
        yy = mean(data(trial).speed.avg(recovery),'omitnan');
        scatter(1,yy,50,'filled', 'XJitter','density','XJitterWidth',1);
    end
    ylim([0 25])
    xticks(1)
    xticklabels([])
    ylabel('Percentage of flies moving during recovery')
    
% Format figure
formatFig(fig,blkbnd,[r,c],sb);
subplot(r,c,sb(1).idx)
set(gca, 'xcolor', 'none')

% Save figure
save_figure(fig,[figDir, 'Avg percentage still flies overtime across trials'], fig_type);

%% FIGURE: Compare tracked vs manual data for ground truthing
clearvars('-except',initial_vars{:})

% Calculate average percentage still 
incap = [];
a = [];
for trial = 1:ntrials 
    incap = autoCat(incap, data(trial).speed.still, false);
    a = [a, length(data(trial).occupancy.time)];
end 
[b,c] = min(a); %[minimun, what column minimum value is located within a]
avgstill = mean(incap(1:b,:),2,'omitnan');
avgpercstill = avgstill*100;

% FIGURE
fig = getfig;
hold on
y = smooth(avgpercstill,sSpan,'moving');
z = (manual.avgperc)';
plot(data(c).occupancy.time, y,'linewidth',2)
plot(manual.time, z,'LineWidth',2)
% xline(manual.time,'Color','r','LineWidth',1)

% Format figure
formatFig(fig,blkbnd);
xlabel('Time (min)')
ylabel('Percentage of still flies')
ylim([0 100])


% Save figure
save_figure(fig,[figDir, 'Tracked vs manual comparison of average stillness'], fig_type);

%% MANUAL FIGURE: Plot number of incapacitated/dead flies overtime
fig = getfig;
for exp = 1:nexps
    hold on
    y_err = (std(manual(exp).raw,0,1,'omitnan'))./sqrt(manual(exp).ntrials); % calculate SEM
    plot_error_fills(1, manual(exp).time, manual(exp).avg, y_err, manual(exp).color, fig_type, 0.35); % plot error
    plot(manual(exp).time,manual(exp).avg,'linewidth',3,'color',manual(exp).color)
end
formatFig(fig, blkbnd);
    ylim([0 20])
    set(gca,"YTick", 0:3:20, 'FontSize', 20)
    xlabel('Time (min)')
    ylabel('Number of incapaciated flies')
save_figure(fig,[figDir, '/' ExpGroup ' death overtime'], fig_type);

%% MANUAL FIGURE: Percentage flies incapacitated/dead overtime

fig = getfig;
for exp = 1:nexps
    hold on
    y_err = (std(manual(exp).perc,0,1,'omitnan'))./sqrt(manual(exp).ntrials); % calculate SEM
    plot_error_fills(1, manual(exp).time/60, manual(exp).avgperc, y_err, manual(exp).color, fig_type, 0.35); % plot error
    plot((manual(exp).time)/60,manual(exp).avgperc,'linewidth',3,'color',manual(exp).color)
end

formatFig(fig, blkbnd);
    ylim([0 100])
    set(gca,"YTick", 0:25:100, 'XTick', 0:6:48, 'FontSize', 15)
    xlabel('Time (hours)')
    ylabel('Percentage of incapaciated flies')
    xlim([0 12])

save_figure(fig,[figDir, '/' ExpGroup ' % death overtime'], fig_type);

%% MANUAL FIGURE: Plot number of alive flies overtime
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


