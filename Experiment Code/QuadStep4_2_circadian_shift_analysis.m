
% Circadian timing comparison across time 

initial_vars{end+1} = 'expstarttime';
initial_vars{end+1} ='clist';
initial_vars{end+1} = 'binedges';
initial_vars{end+1} = 'expstarttime_hr';
initial_vars{end+1} = 'groupidx';
initial_vars{end+1} = 'ngroups';
initial_vars{end+1} = 'foreColor';
initial_vars{end+1} = 'dataString';

%% ANALYSIS: make unique trial ID and pull experiment start time
clearvars('-except',initial_vars{:})
exp = 1; 
% Establish variables for experiment start times
expstarttime = nan([num.trial(exp),1]);

% load excel summary information
[excelfile, Excel, xlFile] = load_QuadBowlExperiments;

T = data(exp).T;

% Make a unique ID for each trial and find experiment start time
for trial = 1:num.trial(exp)
    % Create unique trial ID
    unqID = [T.Date{trial} '_' T.ExperimentID{trial} '_' T.Arena{trial}];
    % Load QuadBowl excel file
    explist = excelfile(:,Excel.trialID);
    % Find row in excel file for each trial
    rownumber = find(strcmp(explist, unqID));

    % extract the zeitgeber time for the experiment: 
    expstarttime(trial) = excelfile{rownumber,Excel.zeitgebertime};
end

% Figure formatting variables
clist = {'Red','Lime','LemonChiffon','Gold','Tomato','DodgerBlue','Teal','Purple'};
% clist = {'Red','Lime','LemonChiffon','Gold','Tomato','DodgerBlue','Teal','MediumPurple'};

% split trials into time bins

% Establish edges for time bins
binedges = 0:3:24;
expstarttime_hr  = floor(expstarttime);

% Create group index, used to split trials into time bins
groupidx = discretize(expstarttime_hr,binedges);

% Make group names for figure legend and labels
dataString = [];
dummy = unique(groupidx);
for idx = 1:length(dummy)
    i = dummy(idx);
    % bincenter = mean(binedges(i:i+1));
    % binwidth = binedges(i+1)-bincenter;
    dataString{idx} = [num2str(binedges(i)) '-' num2str(binedges(i+1))];
end

% Caluclate number of time bin groups
ngroups = length(binedges)-1;

foreColor = formattingColors(blkbgd);

%% FIGURE: plot experienced temperature by chronological time not event aligned
clearvars('-except',initial_vars{:})
exp = 1;
temp_offset = 2; % vertical offset for temperature

roi = 1:66700; % cut off the final release of temp

fig = getfig('',1); hold on
for trial = 1:num.trial(exp)
    pd =  data(exp).data(trial).occupancy;
    % offset time and temp
    x = pd.time(roi) + (60*expstarttime(trial)); % convert zeitbeger time offset to minutes
    y = pd.temp(roi) + (temp_offset*trial);

    % plot
    kolor = Color(clist{groupidx(trial)});
    plot(x,y,'linewidth', 1.5,'Color',kolor)
end
formatFig(fig,blkbgd);
xlabel('Zeitgeber time (min since light on)')
set(gca, 'ycolor','none')

save_figure(fig,[figDir 'Experiment start time timecourse'],fig_type);

%% FIGURE: histogram of experiment start times

fig = getfig('', 1,[460,680]);
    % plot histogram of start times
    h = histogram(expstarttime,0:24);
    % Format figure
    formatFig(fig,blkbgd);
    h.FaceColor = Color('purple');
    h.FaceAlpha = 1;
    h.EdgeColor = foreColor;
    h.LineWidth = 1.5;
    % Create axes labels
    xlabel('Zeitgeber time (hours since light on)')
    ylabel('Count')

save_figure(fig,[figDir 'Experiment start time histogram'],fig_type);

%% FIGURE: total sleep vs start time

fig = getfig('',1,[632 680]);
    x = expstarttime;
    y = sleep.avg_quant;
    scatter(x,y,100,foreColor, "filled",'MarkerFaceAlpha',0.75)
formatFig(fig,blkbgd);
xlabel('Zeitgeber time (hours since light on)')
ylabel('avg sleep per trial')
save_figure(fig,[figDir 'avg sleep per start time'],fig_type);

%% FIGURE: total sleep per timebinned group

LW = 2; % errorbar width
FA = 0.75; % scatter face alpha
SZ = 80; % scatter marker size
BW = 0.8; % bar width

[y,y_sem] = deal(nan([1,ngroups]));

fig = getfig('',1,[632 680]);
hold on
    % extract the avg sleep timing for each group
    for i = 1:ngroups
        y_raw = sleep.avg_quant(groupidx==i);
        y(i) = mean(y_raw,'omitnan');
        y_std(i) = std(y_raw,'omitnan');
        % plot data
        bar(i,y(i),'FaceColor',Color('purple'),'EdgeColor',foreColor,'BarWidth',BW)
        errorbar(i,y(i),y_std(i),'Color',foreColor,'LineWidth',LW)
        scatter(i*ones(size(y_raw)),y_raw, SZ,foreColor,"filled",'MarkerFaceAlpha',FA,'XJitter','density','XJitterWidth',0.6)
    end
formatFig(fig,blkbgd);
set(gca, 'xtick', unique(groupidx),'XTickLabel', dataString)
xlabel('Zeitgeber time (hours since light on)')
ylabel('avg sleep per trial')

save_figure(fig,[figDir 'avg sleep per bined ZGtime'],fig_type);


