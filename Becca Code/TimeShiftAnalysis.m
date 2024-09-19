
%% Load data
clear; clc; close all
path = getDataPath(3, 0, 'Select location of data structure');
folderSelected = selectFolder(path,'Single', 'Select data structure');
figdirectory = [path,folderSelected{:} '/'];
load([figdirectory, folderSelected{:} ' ' 'post 3.1 data.mat'])

initial_vars{end + 1} = 'figdirectory';

%"S:\Evyn\Data structures\Berlin F LRR caviar time shift\Berlin F LRR caviar time shift post 3.1 data.mat"
%% Pull start time data

clearvars('-except',initial_vars{:})

% load excel file:
[excelfile, Excel, XL] = load_QuadBowlExperiments;

initial_vars{end + 1} = 'expstarttime';
initial_vars{end + 1} = 'expstarttime_hr';
initial_vars{end + 1} = 'groupix';
initial_vars{end + 1} = 'binedges';

%% Make unique trial ID
expstarttime = [];
expstarttime_hr = [];
ntrials = size(T,1);

for trial = 1:ntrials
    unqID = [T.Date{trial} '_' T.ExperimentID{trial} '_' T.Arena{trial}];
    explist = excelfile(:,Excel.trialID);
    rownumber = find(strcmp(explist, unqID));
    t = excelfile{rownumber,Excel.starttime};
    dn = excelfile{rownumber,Excel.daynight};

    % create contingency plan for non-duration values
    if ischar(t)
       startdur = duration(t);
    else
        newTime = datetime(t,"ConvertFrom", "excel");
        [h,m,s] = hms(newTime); % Get the hour, minute, and second values
        startdur = duration(h,m,s); % Get the output as a duration() array
    end

    if strcmpi(dn,'N')
        startdur = startdur - duration(7,0,0);
    end

    expstarttime{trial} = startdur;
    expstarttime_hr(trial) = hours(expstarttime{trial});
end


%% FIGURE: Histogram of experiment start times

clearvars('-except',initial_vars{:})

x = expstarttime_hr;
binedges = [4:4:20];

fig = figure('Name','Experiment start time');
h = histogram(expstarttime_hr,binedges);
h.FaceColor = Color('Plum');
h.FaceAlpha = 1;
h.EdgeColor = 'w';
h.LineWidth = 1.5;
xlabel('Time of day')
ylabel('Count')
formatFig(fig,true)

%Save figure
save_figure(fig, [figdirectory 'Experiment start times'], '-png');

%% Split trials into time bins

clearvars('-except',initial_vars{:})

% Split trials into before and after noon

groupix = discretize(expstarttime_hr,binedges);

%% FIGURE: Distance to food between groups with temperature log

clearvars('-except',initial_vars{:}) %clear all variables except initial ones
plot_err = true; %plot error region
[foreColor,backColor] = formattingColors(true); % get background colors

ntrials = size(T,1); %number of trials to include
ngroups = length(binedges)-1;

% set up figure aligments
r = 4; %rows
c = 1; %columns
sb(1).idx = [1]; %temp timecourse
sb(2).idx = [2,3,4]; %distance from food timecourse %TODO: normalize this to something more intuitive?

LW = 1.5; %linewidth
sSpan = 360; %smoothing function

dataString = [];
for i = 1:ngroups
    bincenter = mean(binedges(i:i+1));
    binwidth = binedges(i+1)-bincenter;
    dataString{i} = [num2str(bincenter) '+-' num2str(binwidth)];
end

g = []; %empty variable to add to later and not override things
for i = 1:ngroups
    g(i).name = dataString{i}; 
    [g(i).raw,g(i).temp,g(i).time] = deal([]); %empty variable
end

for trial = 1:ntrials
    loc = data(trial).data.foodwell; %column number that has the food
    y = data(trial).data.dist2well(:,loc); %y variable = distance to food well
    idx = groupix(trial); %group identity for each trial (before or after noon)
    g(idx).raw = autoCat(g(idx).raw,y,false); %organize data into one big data matrix
    x = data(trial).occupancy.time;
    g(idx).time = autoCat(g(idx).time,x,false);
    z = data(trial).occupancy.temp;
    g(idx).temp = autoCat(g(idx).temp,z,false);
end

for i = 1:ngroups %loop through each group
    g(i).distavg = mean(g(i).raw,2,'omitnan'); %mean distance for each group across time
    g(i).std = std(g(i).raw,0,2,'omitnan'); %error for each group 
    g(i).timeavg = mean(g(i).time,2,'omitnan');
    g(i).tempavg = mean(g(i).temp,2,'omitnan');
end

clist = {'DeepPink','Orange','Lime','DodgerBlue'};


 % FIGURE:
fig = getfig('Distance to food',true); %create a figure with customizable size and position on screen
for i = 1:ngroups
subplot(r,c,sb(1).idx)
    hold on
    plot(g(i).timeavg,g(i).tempavg,'color',Color(clist{i}),'LineWidth',LW)
subplot(r,c,sb(2).idx) %dimensions and location of subplot
    hold on %keep what is there to add multiple things
        x = g(i).timeavg;
        y = g(i).distavg;
        y_err =  g(i).std;
       plot(x,y,'color',Color(clist{i}))
       h = plot_error_fills(plot_err, x, y, y_err, Color(clist{i}),'-png', 0.2);
end
    
formatFig(fig,true,[r,c],sb)

subplot(r,c,sb(1).idx)
    set(gca,'xcolor',backColor)
    ylabel('Temperature (\circC)')

 subplot(r,c,sb(2).idx) 
    legend(dataString,'Location','northwest','color',backColor,'box','off','textColor',foreColor)
    xlabel('Time (min)')
    ylabel('Distance to food (mm)')

% Save figure
save_figure(fig, [figdirectory 'Distance to food'], '-png')


%% FIGURE: Distance from food/temperature correlation : TODO -- format this figure

clearvars('-except',initial_vars{:}) %clear all variables except initial ones
plot_err = true; %plot error region
[foreColor,backColor] = formattingColors(true); % get background colors

plotdata = nan([ntrials,2]); %create empty variable of a specific size

for trial = 1:ntrials
    z = data(trial).occupancy.temp;
    loc = data(trial).data.foodwell; %column number that has the food
    y = data(trial).data.dist2well(:,loc); %y variable = distance to food well
    tdmat = [z,y]; %matrix of temperature and distance to food well
    tempPoints = getTempTurnPoints(temp_protocol);
    ROI = [tempPoints.DownROI,tempPoints.UpROI];

    tdmat = tdmat(ROI,:);
    nanloc = isnan(tdmat);
    a = any(nanloc,2);
    tdmat(a,:) = [];

    corr1 = corrcoef(tdmat(:,1),tdmat(:,2));
    plotdata(trial,2) = corr1(1,2);
    plotdata(trial,1) = expstarttime_hr(trial);
end

fig = getfig('Temp distance correlation',true);
    scatter(plotdata(:,1),plotdata(:,2),50,'w',"filled")
    xlabel('Time of day (hr)')
    ylabel('Correlation between distance and temp')

formatFig(fig,true)

% Save figure
save_figure(fig, [figdirectory 'Distance to temp correlation across start times'], '-png')


%% FIGURE: Change in distance during cooling across start times

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true); 

plotdata = nan([ntrials,2]);
windowsize = 1; %in minutes

clist = {'DeepPink','Orange','grey','Purple'};

r = 1;
c = 2;

fig = getfig('Cooling ramp distance difference',true);
hold on

for trial = 1:ntrials
    subplot(r,c,1)
    hold on
    loc = data(trial).data.foodwell; %column number that has the food
    y = data(trial).data.dist2well(:,loc); %y variable = distance to food well
    tempPoints = getTempTurnPoints(temp_protocol);

    nframes = tempPoints.fps*60*windowsize;
    r1start = tempPoints.hold(1,2)-nframes:tempPoints.hold(1,2);
    r1end = tempPoints.down(1,2)-floor(nframes/2):tempPoints.down(1,2)+ceil(nframes/2);
    r3start = tempPoints.hold(3,2)-nframes:tempPoints.hold(3,2);
    r3end = tempPoints.down(3,2)-floor(nframes/2):tempPoints.down(3,2)+ceil(nframes/2);
   
    r1startavg = mean(y(r1start)); %average distance within r1start area
    r1endavg = mean(y(r1end));
    r3startavg = mean(y(r3start));
    r3endavg = mean(y(r3end));

    r1diff = r1endavg-r1startavg;
    r3diff = r3endavg-r3startavg;

    scatter(1,r1diff,35,Color(clist{groupix(trial)}),"filled")
    scatter(3,r3diff,35,Color(clist{groupix(trial)}),"filled")
    plot([1,3],[r1diff,r3diff],'color',Color(clist{groupix(trial)}),'LineWidth',2)

    subplot(r,c,2)
    hold on
    slope = r3diff-r1diff;
    scatter(expstarttime_hr(trial),slope,35,Color(clist{groupix(trial)}),"filled")

end

subplot(r,c,1)
    xlim([0.5,3.5])
    set(gca,'XTick',[1,3])
    xlabel('Ramp')
    ylabel('\Delta Distance during cooling')
    title('Cooling')

subplot(r,c,2)
    xlim([5,20])
    set(gca,'XTick',binedges)
    h_line(0,'Gray','--')
    xlabel('Start time')
    ylabel('Difference between ramp 1 and ramp 3')
    
 formatFig(fig,true,[r,c])

 % Save figure
save_figure(fig, [figdirectory 'Distance difference across start times during cooling'], '-png')

%% FIGURE: Change in distance during heating across start times

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true); 

plotdata = nan([ntrials,2]);
windowsize = 1; %in minutes

clist = {'DeepPink','Red','Green','Blue'};

r = 1;
c = 2;

fig = getfig('Heating ramp distance difference',true);
hold on

for trial = 1:ntrials
    subplot(r,c,1)
    hold on
    loc = data(trial).data.foodwell; %column number that has the food
    y = data(trial).data.dist2well(:,loc); %y variable = distance to food well
    tempPoints = getTempTurnPoints(temp_protocol);

    nframes = tempPoints.fps*60*windowsize;
    r1start = tempPoints.down(1,2)-floor(nframes/2):tempPoints.down(1,2)+ceil(nframes/2);
    r1end = tempPoints.up(1,2):tempPoints.up(1,2)+nframes;
    r3start = tempPoints.down(3,2)-floor(nframes/2):tempPoints.down(3,2)+ceil(nframes/2);
    r3end = tempPoints.up(3,2):tempPoints.up(3,2)+nframes;
   
    r1startavg = mean(y(r1start)); %average distance within r1start area
    r1endavg = mean(y(r1end));
    r3startavg = mean(y(r3start));
    r3endavg = mean(y(r3end));

    r1diff = r1endavg-r1startavg;
    r3diff = r3endavg-r3startavg;

    scatter(1,r1diff,35,Color(clist{groupix(trial)}),"filled")
    scatter(3,r3diff,35,Color(clist{groupix(trial)}),"filled")
    plot([1,3],[r1diff,r3diff],'color',Color(clist{groupix(trial)}),'LineWidth',2)

    subplot(r,c,2)
    hold on
    slope = r3diff-r1diff;
    scatter(expstarttime_hr(trial),slope,35,Color(clist{groupix(trial)}),"filled")

end

subplot(r,c,1)
    xlim([0.5,3.5])
    set(gca,'XTick',[1,3])
    xlabel('Ramp')
    ylabel('\Delta Distance during heating')
    title('Heating')

subplot(r,c,2)
    xlim([5,20])
    set(gca,'XTick',binedges)
    h_line(0,'Gray','--')
    xlabel('Start time')
    ylabel('Difference between ramp 1 and ramp 3')
    
 formatFig(fig,true,[r,c])

 % Save figure
save_figure(fig, [figdirectory 'Distance difference across start times during heating'], '-png')

%% FIGURE: Intitial distance from food before cooling across start times

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true); 

xtickloc = [];
for i = 1:length(binedges)-1
    stats(i).data = [];
end

windowsize = 1; %in minutes
r = 1;
c = 2;

fig = getfig;

clist = {'DeepPink','Orange','grey','Purple'};

for trial = 1:ntrials
    subplot(r,c,1)
    hold on
    loc = data(trial).data.foodwell; %column number that has the food
    y = data(trial).data.dist2well(:,loc); %y variable = distance to food well
    tempPoints = getTempTurnPoints(temp_protocol);

    nframes = tempPoints.fps*60*windowsize;
    r1start = tempPoints.hold(1,2)-nframes:tempPoints.hold(1,2);
    r1end = tempPoints.down(1,2)-floor(nframes/2):tempPoints.down(1,2)+ceil(nframes/2);
    r3start = tempPoints.hold(3,2)-nframes:tempPoints.hold(3,2);
    r3end = tempPoints.down(3,2)-floor(nframes/2):tempPoints.down(3,2)+ceil(nframes/2);
   
    r1startavg = mean(y(r1start)); %average distance within r1start area
    r3startavg = mean(y(r3start));

    stats(groupix(trial)).data(end+1,:) = [r1startavg,r3startavg];
    
    xramp1 = (groupix(trial)*2)-1;
    xramp3 = groupix(trial)*2;

    scatter(xramp1,r1startavg,35,Color(clist{groupix(trial)}),"filled")
    scatter(xramp3,r3startavg,35,Color(clist{groupix(trial)}),"filled")
    plot([xramp1,xramp3],[r1startavg,r3startavg],'color',Color(clist{groupix(trial)}),'LineWidth',2)

    subplot(r,c,2)
    hold on

    xramp1 = groupix(trial);
    xramp3 = groupix(trial) + max(groupix) + 1;

    buff = 0.3;
    nincr = 50;
    poss = linspace(-buff,buff,nincr);
    idx = randi(nincr,1); %random number between 1 and nincrement
    b = poss(idx);

    scatter(xramp1+b,r1startavg,35,Color(clist{groupix(trial)}),"filled")
    scatter(xramp3+b,r3startavg,35,Color(clist{groupix(trial)}),"filled")

    xtickloc = [xtickloc;xramp1;xramp3];
end

buff = 0.5;
xmax = max(unique(groupix))*2 + buff;
xmin = min(unique(groupix))*2 - 1 - buff;



subplot(r,c,1)
    set(gca,'XTick',1:8,'XTickLabel',{'1','3','1','3','1','3'}) %TODO: make tick range dynamic
    xlim([xmin,xmax])
    xlabel('Ramp')
    ylabel('Distance from food (mm)')

dataString = [];
dummy = unique(groupix);
for idx = 1:length(dummy)
    i = dummy(idx);
    % bincenter = mean(binedges(i:i+1));
    % binwidth = binedges(i+1)-bincenter;
    dataString{idx} = [num2str(binedges(i)) '-' num2str(binedges(i+1))];
end

xticklabels = repmat(dataString,[1,2]);

subplot(r,c,2)
    set(gca,'XTick',unique(xtickloc),'XTickLabel',xticklabels,'XTickLabelRotation',30) %TODO: make tick range dynamic
    xlabel('Start time (hr)')
    ylabel('Distance from food (mm)')

formatFig(fig,true,[r,c]);

 % Save figure
save_figure(fig, [figdirectory 'Starting distance from food before cooling'], '-png')

% Stats
h = [];
p = [];
for i = 1:length(binedges)-1 %TODO make into variable 
    z = stats(i).data;
    if isempty(z)
        h(i) = false;
        p(i) = nan;
    else
        [h(i),p(i)] = ttest(z(:,1),z(:,2));
        if h(i)
            disp(['Group ' num2str(i) ' p = ' num2str(p(i))])
        end
    end
end

[r1,r3] = deal([]);

for i = 1:length(binedges)-1
    z = stats(i).data;
    if ~isempty(z)
        r1 = autoCat(r1,z(:,1),false);
        r3 = autoCat(r3,z(:,2),false);
    end
end

[p1,tbl] = anova1(r1);
[p2,tbl2] = anova1(r3);


%% TODO Distance Analyses

% 1) Is there a difference in the initial distance flies are from the food for the 1 & 3rd ramp across time? 
% -- learn about running statistica tests in Matlab

% 2) How does the return to food change during warming between the 1 & 3rd ramps?
% -- what are the implications if they return to different distances for different start times? 

% 3) s there a correlation between start time and distance traveled during cooling? 

% 4) Is there a correlation between start time and the distance to food in the 1&3rd ramps?

%% 
















