
figure;
exp = 1;
trial = 1;
y = ones(size(sROImask(exp).all(trial).m));
y = y.*[1:size(y,2)];
y(~sROImask(exp).all(trial).m) = nan;

plot(y)

%% Load data
clear; clc;
baseFolder = getDataPath(2,0);

blkbnd = true;

% Load excel file
[excelfile, Excel, xlFile] = load_SurvivalQuadCounts;

% Find the experiments that have been counted
loc = cellfun(@isnan,excelfile(2:end,Excel.counted));
loc = ~loc;
rownums = find(loc)+1; 
eligible_files = excelfile([false;loc],[Excel.date, Excel.expID, Excel.arena, Excel.counted]);
FileNames = format_eligible_files(eligible_files);
temptype = unique(excelfile(2:end,5));

% Select files 
fileIdx = listdlg('ListString', temptype,'ListSize',[350,450],'promptstring', 'Select temp protocol to process');
if isempty(fileIdx)
    disp('No trials selected')
    return
end

% % Establish date, trial, and base directories
% dateDir = eligible_files(fileIdx(1),1);
% trialDir = eligible_files(fileIdx(1),2); 
% baseDir = [baseFolder, dateDir{:} '\', trialDir{:} '\'];

% Pull the protocol name and trials under the selected temp protocol
tempprotocol = temptype{fileIdx};
temptrials = strcmp(tempprotocol, excelfile(:,5)); % binary

% Which columns in excel file to plot based on each temp protocol
switch tempprotocol
    case 'survival_curve_40C'
        columns = 11:51;
        k =  'tomato';
    case 'survival_curve_5C'
        columns = [11:13,15,17,19,23,28,33,38,43,48,53,58];
        k = 'dodgerblue';
end

data = [];
data.raw = cell2mat([excelfile(temptrials,columns)]);
data.avg = mean(data.raw,'omitnan');
timestmp = cell2mat([excelfile(1,columns)]);
time = timestmp*10;


% Plot number of dead flies overtime
fig = getfig;
    plot(time,data.avg,'linewidth',3,'color',Color(k))
formatFig(fig, blkbnd);
    ylim([0 15])
    set(gca,"YTick", 0:3:20, 'FontSize', 20)
    xlabel('Time (min)')
    ylabel('Number of dead flies')

% Plot number of alive flies overtime
numflies = cell2mat([excelfile(temptrials,10)]);
aliveflies = numflies - cell2mat([excelfile(temptrials,11:51)]);
% perc_alive = (aliveflies/numflies)*100;

fig = getfig;
    plot(time,aliveflies,'linewidth',2)
formatFig(fig, blkbnd);
    ylim([0 20])
    % xlim([0 500])
    xlabel('Time (min)')
    ylabel('Number of alive flies')


