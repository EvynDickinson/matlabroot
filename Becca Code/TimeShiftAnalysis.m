
%% LOAD: data
clear; clc; close all
path = getDataPath(3, 0, 'Select location of data structure');
folderSelected = selectFolder(path,'Single', 'Select data structure');
figdirectory = [path,folderSelected{:} '/'];
load([figdirectory, folderSelected{:} ' ' 'post 3.1 data.mat'])

initial_vars{end + 1} = 'figdirectory';

disp('next section')

%"S:\Evyn\Data structures\Berlin F LRR caviar time shift\Berlin F LRR caviar time shift post 3.1 data.mat"
%% LOAD: start time data

clearvars('-except',initial_vars{:})

% load excel file:
[excelfile, Excel, XL] = load_QuadBowlExperiments;

initial_vars{end + 1} = 'expstarttime';
initial_vars{end + 1} = 'expstarttime_hr';
initial_vars{end + 1} = 'groupix';
initial_vars{end + 1} = 'binedges';
initial_vars{end + 1} = 'ngroups';
initial_vars{end + 1} = 'clist';
initial_vars{end + 1} = 'bkgrd_color';

disp('next section')

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

binedges = [4:2:20];
clist = {'Red','Lime','LemonChiffon','Gold','Tomato','DodgerBlue','Teal','Purple'};

disp('next section')

%% FIGURE: Histogram of start times

clearvars('-except',initial_vars{:})

bkgrd_color = false;

x = expstarttime_hr;

fig = figure('Name','Experiment start time');
h = histogram(expstarttime_hr,binedges);
h.FaceColor = Color('Plum');
h.FaceAlpha = 1;
h.EdgeColor = 'w';
h.LineWidth = 1.5;
xlabel('Time of day')
ylabel('Count')
formatFig(fig,bkgrd_color)

%Save figure
save_figure(fig, [figdirectory 'Experiment start times histogram'], '-png');


%% FIGURE: Polar plot of start times

clearvars('-except',initial_vars{:})

% Convert the 24 hours into a circular radian value
theta = (24-(expstarttime_hr/24)) * 2 * pi; 

% show the experiment start times as a polar plot
fig = getfig('Polar Histogram', 1);
    h = polarhistogram(theta);
    set(h,'FaceColor',Color('teal'), 'EdgeColor', 'k','FaceAlpha',0.8)
    % set labels and formats
    ax = gca;
    set(ax, 'ThetaZeroLocation', 'top')
    set(ax, 'ThetaTick',0:45:359,'ThetaTickLabel', {'0','21','18','15','12','9','6','3'})

% OVERLAY NIGHT SHADING
ax_cart = axes('Position', ax.Position, 'Color', 'none');  % Create an overlay Cartesian axes

% Step 5: Shade the left half (90° to 270°)
hold(ax_cart, 'on');

% Define the angles and radius for shading the left half
shading_theta = linspace(pi/2, 3*pi/2, 100);  % Angles from 90° to 270°
shading_r = max(ax.RLim) * ones(size(shading_theta));  % Use maximum radius of polar plot

% Convert polar coordinates to Cartesian for the shading
[x_shade, y_shade] = pol2cart(shading_theta, shading_r);

% Add the center point to complete the patch
x_shade = [0 x_shade 0];  % Include the origin (0,0)
y_shade = [0 y_shade 0];

% Create the filled patch in grey
fill(ax_cart, x_shade, y_shade, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');  % Grey with transparency

% Adjust Cartesian axis limits to match polar plot
axis(ax_cart, 'equal');  % Ensure equal scaling on both axes
r_max = max(ax.RLim);  % Get the maximum radius from the polar plot
ax_cart.XLim = [-r_max r_max];
ax_cart.YLim = [-r_max r_max];

% Remove the Cartesian axis ticks
ax_cart.XTick = [];
ax_cart.YTick = [];
ax_cart.XColor = 'none';
ax_cart.YColor = 'none';

% Bring the polar plot to the front  %TODO: see if there is a way to do
% this without blocking the new shading
uistack(ax, 'top');

%Save figure
save_figure(fig, [figdirectory 'Experiment start times polar plot'], '-png');

% TODO: update the label sizes and add some titles etc, 
% TODO: try making the binsizes for each hour

%% Split trials into time bins

clearvars('-except',initial_vars{:})

% Split trials into before and after noon

groupix = discretize(expstarttime_hr,binedges);
ntrials = size(T,1); %number of trials to include
ngroups = length(binedges)-1;

disp('next section')

%% FIGURE: Distance to food across start times with temperature log

clearvars('-except',initial_vars{:}) %clear all variables except initial ones
plot_err = true; %plot error region
[foreColor,backColor] = formattingColors(true); % get background colors

% set up figure aligments
r = 4; %rows
c = 1; %columns
sb(1).idx = [1]; %temp timecourse
sb(2).idx = [2,3,4]; %distance from food timecourse %TODO: normalize this to something more intuitive?

LW = 1.5; %linewidth
sSpan = 360; %smoothing function

dataString = [];
dummy = unique(groupix);
for idx = 1:length(dummy)
    i = dummy(idx);
    % bincenter = mean(binedges(i:i+1));
    % binwidth = binedges(i+1)-bincenter;
    dataString{idx} = [num2str(binedges(i)) '-' num2str(binedges(i+1))];
end

g = []; %empty variable to add to later and not override things
for i = 1:ngroups
    g(i).name = [num2str(binedges(i)) '-' num2str(binedges(i+1))]; 
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


%% FIGURE: Distance to food/temperature correlation

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


%% FIGURE: Change in distance during cooling | heating across start times

% TODO:
% display means on right figure
% separate by start time for left figure
% average lines for each start time

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true);

inputVar =  questdlg('Which data type to compare?','','Cooling','Heating','Cancel');

switch inputVar
    case 'Cooling'
        ylab = '\Delta Distance during cooling';
    case 'Heating'
        ylab = '\Delta Distance during heating';
    case 'Cancel'
        return
end

plotdata = nan([ntrials,2]);
windowsize = 1; %in minutes

r = 1;
c = 2;

fig = getfig(['Distance difference across start times during ', inputVar],true);
hold on

for trial = 1:ntrials
    subplot(r,c,1)
        hold on
        loc = data(trial).data.foodwell; %column number that has the food
        y = data(trial).data.dist2well(:,loc); %y variable = distance to food well
        tempPoints = getTempTurnPoints(temp_protocol);
    
        nframes = tempPoints.fps*60*windowsize;
        switch inputVar
            case 'Cooling'
                r1start = tempPoints.hold(1,2)-nframes:tempPoints.hold(1,2);
                r1end = tempPoints.down(1,2)-floor(nframes/2):tempPoints.down(1,2)+ceil(nframes/2);
                r3start = tempPoints.hold(3,2)-nframes:tempPoints.hold(3,2);
                r3end = tempPoints.down(3,2)-floor(nframes/2):tempPoints.down(3,2)+ceil(nframes/2);
            case 'Heating'
                r1start = tempPoints.down(1,2)-floor(nframes/2):tempPoints.down(1,2)+ceil(nframes/2);
                r1end = tempPoints.up(1,2):tempPoints.up(1,2)+nframes;
                r3start = tempPoints.down(3,2)-floor(nframes/2):tempPoints.down(3,2)+ceil(nframes/2);
                r3end = tempPoints.up(3,2):tempPoints.up(3,2)+nframes;
        end
       
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
        avgslope(trial) = mean(slope);
        scatter(expstarttime_hr(trial),slope,35,Color(clist{groupix(trial)}),"filled")

end

subplot(r,c,1)
    xlim([0.5,3.5])
    set(gca,'XTick',[1,3])
    xlabel('Ramp')
    ylabel(ylab)
    title(inputVar)

subplot(r,c,2)
    xlim([5,20])
    set(gca,'XTick',binedges)
    h_line(0,'Gray','--')
    xlabel('Start time')
    ylabel('Difference between ramp 1 and ramp 3')
    
 formatFig(fig,true,[r,c])

 % Save figure
save_figure(fig, [figdirectory 'Distance difference across start times during ', inputVar], '-png')


%% ***OLD FIGURE: Change in distance during cooling across start times 

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true); 

plotdata = nan([ntrials,2]);
windowsize = 1; %in minutes

% clist = {'LemonChiffon','Gold','Tomato','Blue'};

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

%% ***OLD FIGURE: Change in distance during heating across start times

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true); 

plotdata = nan([ntrials,2]);
windowsize = 1; %in minutes

% clist = {'LemonChiffon','Gold','Tomato','Blue'};

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

%% FIGURE: Intitial distance from food before cooling across start times : TODO -- stats formatting

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true); 

xtickloc = [];
for i = 1:ngroups
    stats(i).data = [];
end

windowsize = 1; %in minutes
r = 1;
c = 2;

fig = getfig;

% clist = {'LemonChiffon','Gold','Tomato','Blue'};

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
% save_figure(fig, [figdirectory 'Starting distance from food before cooling'], '-png')

% Stats
h = [];
p = [];

% Compare ramp 1 and ramp 3 for each start time
for i = 1:ngroups %TODO make into variable 
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

% Compare each start time for ramp 1 and ramp 3
for i = 1:ngroups
    z = stats(i).data;
    if ~isempty(z)
        r1 = autoCat(r1,z(:,1),false);
        r3 = autoCat(r3,z(:,2),false);
    end
end

[p1,tbl] = anova1(r1,unique(groupix),'off');
[p2,tbl2] = anova1(r3,unique(groupix),'off');

% Add stats to figure

figure(fig)
subplot(r,c,1)
y = rangeLine(fig,2,true);

for i = 1:ngroups
    xramp1 = (i*2)-1;
    xramp3 = i*2;
    % Plot average line
    avgchange = mean(stats(i).data);
    plot([xramp1,xramp3],avgchange,'color',foreColor,'LineWidth',4,'Marker','o')
    if h(i) % h is a logical, h(i) = 1 = true
        scatter((xramp1+xramp3)/2,y,150,foreColor,"filled","*","MarkerEdgeColor",foreColor)
    end
end










%% FIGURE: Intitial distance from food before heating across start times : TODO -- stats formatting

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true); 

xtickloc = [];
for i = 1:ngroups
    stats(i).data = [];
end

windowsize = 1; %in minutes
r = 1;
c = 2;

fig = getfig;

% clist = {'LemonChiffon','Gold','Tomato','Blue'};

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
for i = 1:ngroups %TODO make into variable 
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

for i = 1:ngroups
    z = stats(i).data;
    if ~isempty(z)
        r1 = autoCat(r1,z(:,1),false);
        r3 = autoCat(r3,z(:,2),false);
    end
end

[p1,tbl] = anova1(r1);
[p2,tbl2] = anova1(r3);

%% ***TODO Distance Analyses

% 1) Is there a difference in the initial distance flies are from the food for the 1 & 3rd ramp across time? 
% -- learn about running statistica tests in Matlab

% 2) How does the return to food change during warming between the 1 & 3rd ramps?
% -- what are the implications if they return to different distances for different start times? 

% 3) Is there a correlation between start time and distance traveled during cooling? 

% 4) Is there a correlation between start time and the distance to food in the 1&3rd ramps?



%% LOAD: sleep data

path = getDataPath(1,0);

for trial = 1:ntrials
    sleeppath = [path, T.Date{trial}, '_', T.ExperimentID{trial}, '_', T.Arena{trial}, '\', T.ExperimentID{trial}, ' sleeping data.mat'];
    dummy = load(sleeppath);
    data(trial).sleep = dummy.sleeping;
    disp( T.ExperimentID{trial})
end

%% FIGURE: Percentage of flies sleeping over time

clearvars('-except',initial_vars{:}) %clear all variables except initial ones
plot_err = true; %plot error region
[foreColor,backColor] = formattingColors(true); % get background colors

% set up figure aligments
r = 4; %rows
c = 1; %columns
sb(1).idx = [1]; %temp timecourse
sb(2).idx = [2,3,4]; %distance from food timecourse %TODO: normalize this to something more intuitive?

LW = 1.5; %linewidth
sSpan = 360; %smoothing function

dataString = [];
dummy = unique(groupix);
for idx = 1:length(dummy)
    i = dummy(idx);
    % bincenter = mean(binedges(i:i+1));
    % binwidth = binedges(i+1)-bincenter;
    dataString{idx} = [num2str(binedges(i)) '-' num2str(binedges(i+1))];
end

% g = group
g = []; %empty variable to add to later and not override things
for i = 1:ngroups
    g(i).name = [num2str(binedges(i)) '-' num2str(binedges(i+1))]; 
    [g(i).raw,g(i).temp,g(i).time] = deal([]); %empty variable
end

for trial = 1:ntrials
    y = data(trial).sleep.sleepNum;
    y = (y/T.NumFlies(trial))*100;
    idx = groupix(trial); 
    g(idx).raw = autoCat(g(idx).raw,y,false); %organize data into one big data matrix
    x = data(trial).occupancy.time;
    g(idx).time = autoCat(g(idx).time,x,false);
    z = data(trial).occupancy.temp;
    g(idx).temp = autoCat(g(idx).temp,z,false);
end

for i = 1:ngroups 
    g(i).sleepavg = mean(g(i).raw,2,'omitnan');
    g(i).std = std(g(i).raw,0,2,'omitnan');  
    g(i).timeavg = mean(g(i).time,2,'omitnan');
    g(i).tempavg = mean(g(i).temp,2,'omitnan');
end

clist = {'LemonChiffon','Gold','Tomato','DodgerBlue'};

 % FIGURE:
fig = getfig('Percentage of flies sleeping',true); %create a figure with customizable size and position on screen
for i = 1:ngroups
subplot(r,c,sb(1).idx)
    hold on
    plot(g(i).timeavg,g(i).tempavg,'color',Color(clist{i}),'LineWidth',LW)
subplot(r,c,sb(2).idx) %dimensions and location of subplot
    hold on %keep what is there to add multiple things
        x = g(i).timeavg;
        y = g(i).sleepavg;
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
    ylabel('Percentage of flies sleeping')
    ylim([0,25])

% Save figure
save_figure(fig, [figdirectory 'Percentage of flies sleeping'], '-png')

% TODO
% compare sleep across hold periods
% distance to food while sleeping (histogram)

%% FIGURE: Average sleep per fly by start time

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true); 

plotdata = nan([ntrials,2]); %create empty variable of a specific size

for trial = 1:ntrials
    y = data(trial).sleep.sleepNum;
    y = (y/T.NumFlies(trial))*100;
    sleepavg = mean(y);
    plotdata(trial,2) = sleepavg;
    plotdata(trial,1) = expstarttime_hr(trial);
end

fig = getfig('Percent flies sleeping across start times',true);
    scatter(plotdata(:,1),plotdata(:,2),50,'y',"filled")
    xlabel('Time of day (hr)')
    ylabel('Percentage of flies sleeping')

formatFig(fig,true)

% Save figure
save_figure(fig, [figdirectory 'Distance to temp correlation across start times'], '-png')


%% FIGURE

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(true); 
toffset = [100/60, 212/60, 324/60];
windowsize = 5; %in minutes

fig = getfig;
hold on

for trial = 1:ntrials
    hold on
    loc = data(trial).data.foodwell; %column number that has the food
    y = data(trial).data.dist2well(:,loc); %y variable = distance to food well
    tempPoints = getTempTurnPoints(temp_protocol);
    
    nframes = tempPoints.fps*60*windowsize;
    for r = 1:3
        r1start = tempPoints.hold(r,2)-nframes:tempPoints.hold(r,2);
        r1end = tempPoints.down(r,2)-floor(nframes/2):tempPoints.down(r,2)+ceil(nframes/2);
        r1startavg = mean(y(r1start));
        zeit = expstarttime_hr(trial) - 8 + toffset(r);
        scatter(zeit,r1startavg,30,Color(clist{groupix(trial)}),'filled')
    end

end

formatFig(fig,true);




