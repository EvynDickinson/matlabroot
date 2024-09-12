
%% Pull start time data

clearvars('-except',initial_vars{:})

% load excel file:
[excelfile, Excel, XL] = load_QuadBowlExperiments;

initial_vars{end + 1} = 'expstarttime';
initial_vars{end + 1} = 'expstarttime_hr';
initial_vars{end + 1} = 'groupix';

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

fig = figure;
histogram(expstarttime_hr)

%% Split trials into before noon and after noon

clearvars('-except',initial_vars{:})

% Split trials into before and after noon

binedges = [4,12,20];
groupix = discretize(expstarttime_hr,binedges);


%% Compare distance to food between groups

clearvars('-except',initial_vars{:}) %clear all variables except initial ones
plot_err = true; %plot error region
[foreColor,backColor] = formattingColors(true); % get background colors

ntrials = size(T,1); %number of trials to include

% set up figure aligments
r = 4; %rows
c = 1; %columns
sb(1).idx = [1]; %temp timecourse
sb(2).idx = [2,3,4]; %distance from food timecourse %TODO: normalize this to something more intuitive?

LW = 1.5; %linewidth
sSpan = 360; %smoothing function
dataString = {'before noon', 'after noon'}; %legend

g = []; %empty variable to add to later and not override things
g(1).name = 'before noon'; %group 1 = before noon
g(2).name = 'after noon'; %group 2 = after noon
g(1).raw = []; %empty variable
g(2).raw = []; %empty variable
g(1).time = [];
g(2).time = [];
g(1).temp = [];
g(2).temp = [];

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

for i = 1:2 %loop through each group
    g(i).distavg = mean(g(i).raw,2,'omitnan'); %mean distance for each group across time
    g(i).std = std(g(i).raw,0,2,'omitnan'); %error for each group 
    g(i).timeavg = mean(g(i).time,2,'omitnan');
     g(i).tempavg = mean(g(i).temp,2,'omitnan');
end

    clist = {'blue','lime'};
% FIGURE:
fig = getfig('Distance to food',true); %create a figure with customizable size and position on screen
subplot(r,c,sb(2).idx) %dimensions and location of subplot
    hold on %keep what is there to add multiple things
    for i = 1:2 
        x = g(i).timeavg;
        y = g(i).distavg;
        y_err =  g(i).std;
      plot(x,y,'color',Color(clist{i}))
       h = plot_error_fills(plot_err, x, y, y_err, Color(clist{i}),'-png', 0.4);
      
    end

subplot(r,c,sb(1).idx)
    hold on
    for i = 1:2
    plot(g(i).timeavg,g(i).tempavg,'color',Color(clist{i}),'LineWidth',LW)
    end
    
    ylabel('Temperature (\circC)')


formatFig(fig,true,[r,c],sb)

subplot(r,c,sb(2).idx) 
 legend(dataString,'Location','northwest','color',backColor,'box','off','textColor',foreColor)
    xlabel('Time (min)')
    ylabel('Distance to food (mm)')

subplot(r,c,sb(1).idx)
    set(gca,'xcolor',backColor)














