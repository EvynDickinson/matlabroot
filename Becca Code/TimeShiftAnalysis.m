
%% Load data

path = getDataPath(3, 0, 'Select data structure');
folderSelected = selectFolder(path);
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


%% Histogram of experiment start times

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



%% Compare distance to food between groups

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
       h = plot_error_fills(plot_err, x, y, y_err, Color(clist{i}),'-png', 0.4);
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


%% Trial by trial scatterplot

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
    scatter(plotdata(:,1),plotdata(:,2))


%% other figure

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true); 

plotdata = nan([ntrials,2]);
windowsize = 1; %in minutes

clist = {'DeepPink','Orange','grey','Purple'};

r = 1;
c = 2;

fig = getfig('Ramp distance difference',true);
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

subplot(r,c,2)
    xlim([5,20])
    set(gca,'XTick',binedges)
    h_line(0,'Gray','--')
    xlabel('Start time')
    ylabel('Difference between ramp 1 and ramp 3')
    
 formatFig(fig,true,[r,c])




















