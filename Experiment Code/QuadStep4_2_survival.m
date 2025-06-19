
figure;
exp = 1;
trial = 1;
y = ones(size(sROImask(exp).all(trial).m));
y = y.*[1:size(y,2)];
y(~sROImask(exp).all(trial).m) = nan;

plot(y)

%% Load data
clear; clc;
% baseFolder = getDataPath(2,2);

blkbnd = false;
fig_type = '-pdf';

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
    case 'survival_curve_5C'
        data(exp).columns = [11:14,16,18,20,22,24,27,29,32,34,37,39,42,44,47,49,52,54,57,59];
        data(exp).color = Color('dodgerblue');
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


%% Plot number of dead flies overtime
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

%% Percentage flies dead overtime

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


