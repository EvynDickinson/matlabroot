%% LOAD: data
clear; clc; close all

% Find data structure folder
path = getDataPath(3, 0, 'Select location of data structure');
% Select specific data structure
folderSelected = selectFolder(path,'Single', 'Select data structure');
% Establish path to specific data structure
figdirectory = [path,folderSelected{:} '/'];

% Load data
disp('data loading...')
tic
load([figdirectory, folderSelected{:} ' ' 'post 3.1 data.mat'])
toc

% Set initial variables
initial_vars{end + 1} = 'figdirectory';
% initial_vars{end + 1} = 'expstarttime'; don't think this is needed as an initial variable
initial_vars{end + 1} = 'expstarttime_hr';
initial_vars{end + 1} = 'groupix';
initial_vars{end + 1} = 'binedges';
initial_vars{end + 1} = 'ngroups';
initial_vars{end + 1} = 'clist';
initial_vars{end + 1} = 'bkgrd_color';

disp('data loaded')

% "S:\Evyn\Data structures\Berlin F LRR caviar time shift\Berlin F LRR caviar time shift post 3.1 data.mat"
%% LOAD: start time data

clearvars('-except',initial_vars{:})

% Load QuadBowl Experiments excel file:
[excelfile, Excel, XL] = load_QuadBowlExperiments;

disp('start time data loaded')

%% ANALYSIS: make unique trial ID and pull experiment start time

% Establish variables for experiment start times
expstarttime = [];
expstarttime_hr = [];
% Calculate number of trials
ntrials = size(T,1);

% Make a unique ID for each trial and find experiment start time
for trial = 1:ntrials
    % Create unique trial ID
    unqID = [T.Date{trial} '_' T.ExperimentID{trial} '_' T.Arena{trial}];
    % Load QuadBowl excel file
    explist = excelfile(:,Excel.trialID);
    % Find row in excel file for each trial
    rownumber = find(strcmp(explist, unqID));
    % Identify start time and whether it was a day or night experiment
    t = excelfile{rownumber,Excel.starttime};
    dn = excelfile{rownumber,Excel.daynight};

    % If start time is reported as a character, create start time variable (h:m:s)
    if ischar(t)
       startdur = duration(t);
    % If not a character, get h,m,s values to create start time variable
    else
        % Convert numeric values to datetime values
        newTime = datetime(t,"ConvertFrom", "excel");
        [h,m,s] = hms(newTime); 
        startdur = duration(h,m,s);
    end

    % If night experiment, subtract 7 hours to align with day start times
    if strcmpi(dn,'N')
        startdur = startdur - duration(7,0,0);
    end

    % Calculate experiment start time (hr) for each trial
    expstarttime{trial} = startdur;
    expstarttime_hr(trial) = hours(expstarttime{trial});
end

% Figure formatting variables
clist = {'Red','Lime','LemonChiffon','Gold','Tomato','DodgerBlue','Teal','Purple'};
bkgrd_color = true;

disp('next section')

%% ANALYSIS: split trials into time bins

clearvars('-except',initial_vars{:})

% Establish edges for time bins
binedges = [4:2:20];

% Create group index, used to split trials into time bins
groupix = discretize(expstarttime_hr,binedges);

% Caluclate number of time bin groups
ngroups = length(binedges)-1;

disp('next section')

%% FIGURE: Histogram of start times

clearvars('-except',initial_vars{:})

% Establish x variable
x = expstarttime_hr;

% FIGURE
% Create figure
fig = figure('Name','Experiment start time');
% Plot frequency of experiment start times across all trials
h = histogram(expstarttime_hr,binedges);

% Format figure
formatFig(fig,bkgrd_color);
h.FaceColor = Color('Plum');
h.FaceAlpha = 1;
h.EdgeColor = 'w';
h.LineWidth = 1.5;
% Create axes labels
xlabel('Time of day')
ylabel('Count')

% Save figure
save_figure(fig, [figdirectory 'Experiment start times histogram'], '-png');


%% FIGURE: Polar plot of start times -- TODO: make background color dynamic

clearvars('-except',initial_vars{:})

% Convert the 24 hours into a circular radian value
theta = (24-(expstarttime_hr/24)) * 2 * pi; 

% Show the experiment start times as a polar plot
fig = getfig('Polar Histogram', 1);
    h = polarhistogram(theta);
    set(h,'FaceColor',Color('Plum'), 'EdgeColor', 'k','FaceAlpha',0.8)
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

% Save figure
save_figure(fig, [figdirectory 'Experiment start times polar plot'], '-png');

% TODO: update the label sizes and add some titles etc, 
% TODO: try making the binsizes for each hour


%% FIGURE: Distance to food across start times with temperature log

clearvars('-except',initial_vars{:})
plot_err = true; % plot error region
[foreColor,backColor] = formattingColors(bkgrd_color); % get background colors

% Set up figure aligments
r = 4; % rows
c = 1; % columns
sb(1).idx = [1]; % temp timecourse
sb(2).idx = [2,3,4]; % distance from food timecourse %TODO: normalize this to something more intuitive?

LW = 1.5; % linewidth
sSpan = 180; % smoothing function

% Make group names for figure legend and labels
dataString = [];
dummy = unique(groupix);
for idx = 1:length(dummy)
    i = dummy(idx);
    % bincenter = mean(binedges(i:i+1));
    % binwidth = binedges(i+1)-bincenter;
    dataString{idx} = [num2str(binedges(i)) '-' num2str(binedges(i+1))];
end

% Create empty variable to hold data
g = []; 

% Set up matrices for data
for i = 1:ngroups
    % Give each group a name based on time bins
    g(i).name = [num2str(binedges(i)) '-' num2str(binedges(i+1))];
    % Set up empty matrices for each group
    [g(i).raw,g(i).temp,g(i).time] = deal([]);
end

% Pull and organize data
for trial = 1:ntrials
    loc = data(trial).data.foodwell; % column number that has the food
    % Establish plot variables
    y = data(trial).data.dist2well(:,loc);
    x = data(trial).occupancy.time;
    z = data(trial).occupancy.temp;
    % Identify the group index for each trial (which time bin it falls into)
    idx = groupix(trial);
    % Load data into each matrix within variable g
    g(idx).raw = autoCat(g(idx).raw,y,false);
    g(idx).time = autoCat(g(idx).time,x,false);    
    g(idx).temp = autoCat(g(idx).temp,z,false);
end

% Calculate average distance, time, and temp for each group, and standard deviation for distance
for i = 1:ngroups
    g(i).distavg = mean(g(i).raw,2,'omitnan'); 
    g(i).timeavg = mean(g(i).time,2,'omitnan');
    g(i).tempavg = mean(g(i).temp,2,'omitnan');
    g(i).std = std(g(i).raw,0,2,'omitnan'); 
end

% FIGURE
% Create figure and establish name, screen position for viewing, and size
fig = getfig('distance to food',true);
    for i = 1:ngroups
        % Plot temperature timecourse
        subplot(r,c,sb(1).idx)
            hold on
            plot(g(i).timeavg,g(i).tempavg,'color',Color(clist{i}),'LineWidth',LW)
        % Plot distance to food timecourse
        subplot(r,c,sb(2).idx)
            hold on 
                x = g(i).timeavg;
                y = smooth(g(i).distavg,sSpan, 'moving');
                % y_err = g(i).std;
            plot(x,y,'color',Color(clist{i}),'LineWidth',0.75)
            % Display error area
            % h = plot_error_fills(plot_err, x, y, y_err, Color(clist{i}),'-png', 0.1);
    end

% Format figure (call fig variable first, it's usually just called fig)
formatFig(fig,bkgrd_color,[r,c],sb);

% Formatting for temperature timecourse
subplot(r,c,sb(1).idx)
    % Make x axis same color as background (gca = get current axis)
    set(gca,'xcolor',backColor)
    % Create y axis label
    ylabel('Temperature (\circC)')

% Formatting for distance to food timecourse
 subplot(r,c,sb(2).idx) 
    % Create legend
    legend(dataString,'Location','northwest','color',backColor,'box','off','textColor',foreColor)
    % Create axes labels
    xlabel('Time (min)')
    ylabel('Distance to food (mm)')

% Save figure
save_figure(fig, [figdirectory 'Distance to food'], '-png');
% NOTE: the first naming of a figure with getfig only names the pop up
% figure before saving. The second naming with save_figure is the name that
% gets saved in the folder

%% FIGURE: Distance to food/temperature correlation during ramps

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(bkgrd_color);

% Create matrix to hold plot data
plotdata = nan([ntrials,2]);

% Create figure
fig = getfig('Temp distance correlation',true);
    hold on

for trial = 1:ntrials
    % Establish variables for correlation calculation
    z = data(trial).occupancy.temp;
    loc = data(trial).data.foodwell; 
    y = data(trial).data.dist2well(:,loc); 
    % Create matrix to hold correlation variables
    tdmat = [z,y];
    % Identify region of interest (frames during ramps)
    tempPoints = getTempTurnPoints(temp_protocol);
    ROI = [tempPoints.DownROI,tempPoints.UpROI];

    % Restrict data within matrix of correlation variables to only include data during ramps 
    tdmat = tdmat(ROI,:);
    % Identify NANs and replace them with []
    nanloc = isnan(tdmat);
    a = any(nanloc,2);
    tdmat(a,:) = [];
    
    % Calculate correlation coefficient
    corr1 = corrcoef(tdmat(:,1),tdmat(:,2));
    % Load correlation value and start time into pot data matrix
    plotdata(trial,2) = corr1(1,2);
    plotdata(trial,1) = expstarttime_hr(trial);
    
    % FIGURE
        scatter(plotdata(trial,1),plotdata(trial,2),50,Color(clist{groupix(trial)}),"filled")
        % Create axes labels
        xlabel('Experiment start time (hr)')
        ylabel('Correlation between distance and temp')
end

% Format figure
formatFig(fig,bkgrd_color);

% Save figure
save_figure(fig, [figdirectory 'Distance to temp correlation during ramps across start times'], '-png');


%% FIGURE: Change in distance during cooling | heating across start times

% TODO:
% display means on right figure
% separate by start time for left figure
% average lines for each start time
% show only mean difference line for each time bin to make left graph less messy

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(bkgrd_color);

% Ask type of data you want to look at: cooling or heating ramps
inputVar =  questdlg('Which data type to compare?','','Cooling','Heating','Cancel','Cooling');
% Change y axis labels depending on answer
switch inputVar
    case 'Cooling'
        ylab = '\Delta Distance during cooling';
    case 'Heating'
        ylab = '\Delta Distance during heating';
    case 'Cancel'
        return
end

% Period of time before/after ramp to look at distance from food
windowsize = 1; % minutes

% Set up figure coordinates
r = 1;
c = 2;

% FIGURE
% Create figure
fig = getfig(['Distance difference across start times during ', inputVar],true);
hold on
for trial = 1:ntrials
    % Plot change in distance during ramp 1 and ramp 3 (left)
    subplot(r,c,1)
        hold on
        % Establish y variable
        loc = data(trial).data.foodwell; 
        y = data(trial).data.dist2well(:,loc); 
        
        % Calculate number of frames within the established windowsize
        tempPoints = getTempTurnPoints(temp_protocol);
        nframes = tempPoints.fps*60*windowsize;
        % Identify frames in ramp start and end windows
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
       
        % Calculate average distance from food within each window
        r1startavg = mean(y(r1start));
        r1endavg = mean(y(r1end));
        r3startavg = mean(y(r3start));
        r3endavg = mean(y(r3end));
        % Calculate change in average distance at start and end for each ramp
        r1diff = r1endavg-r1startavg;
        r3diff = r3endavg-r3startavg;
    
        % Plot change in distance where x = ramp number
        scatter(1,r1diff,35,Color(clist{groupix(trial)}),"filled")
        scatter(3,r3diff,35,Color(clist{groupix(trial)}),"filled")
        % Plot line connecting ramps in the same trial
        plot([1,3],[r1diff,r3diff],'color',Color(clist{groupix(trial)}),'LineWidth',2)
    
    % Plot rate of change in distance between ramps (slope of connecting line) (right)
    subplot(r,c,2)
        hold on
        % Calculate change in distance difference between ramps
        slope = r3diff-r1diff;
        % Plot this change across start times
        scatter(expstarttime_hr(trial),slope,35,Color(clist{groupix(trial)}),"filled")

end

% Format figure
formatFig(fig,bkgrd_color,[r,c]);

% Formatting for change in distance at ramps plot
subplot(r,c,1)
    % Set x axis limits
    xlim([0.5,3.5])
    % Set x axis ticks to only 1 and 3, for each ramp
    set(gca,'XTick',[1,3])
    % Create axes labels
    xlabel('Ramp')
    ylabel(ylab)

% Formatting for rate of change in distance plot
subplot(r,c,2)
    % Set x axis limits
    xlim([5,20])
    % Set x axis ticks to established time bin edges
    set(gca,'XTick',binedges)
    % Display a line at y = 0
    h_line(0,'Gray','--')
    % Create axes labels
    xlabel('Start time')
    ylabel('Difference between ramp 1 and ramp 3')
    
% Save figure
save_figure(fig, [figdirectory 'Distance difference across start times during ', inputVar], '-png');


%% FIGURE: Starting distance to food before cooling | heating across start times

% TODO: fix ramp x axis labels
% TODO: Add figure titles to distinguish between left and right graphs

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(bkgrd_color); 

% Ask type of data you want to look at: cooling or heating ramps
inputVar =  questdlg('Which data type to compare?','','Cooling','Heating','Cancel','Cooling');

% Create variable to set up x axis tick locations 
xtickloc = [];

% Create matrix to hold stats data 
for i = 1:ngroups
    stats(i).data = [];
end

% Period of time before/after ramp to look at distance from food
windowsize = 1; % minutes

% Set up figure coordinates
r = 1;
c = 2;

% FIGURE
% Create figure
fig = getfig;
for trial = 1:ntrials
    % Plot distance from food before each ramp within each time bin (left)
    subplot(r,c,1)
        hold on
        % Establish y variable
        loc = data(trial).data.foodwell; 
        y = data(trial).data.dist2well(:,loc); 
        
        % Calculate number of frames within the established windowsize
        tempPoints = getTempTurnPoints(temp_protocol);
        nframes = tempPoints.fps*60*windowsize;
        % Identify frames in ramp start window
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
        
        % Calculate average distance from food within each window
        r1startavg = mean(y(r1start));
        r3startavg = mean(y(r3start));
    
        % Load average distance values for each time bin into stats matrix
        stats(groupix(trial)).data(end+1,:) = [r1startavg,r3startavg];
        
        % Establish x axis tick locations
        xramp1 = (groupix(trial)*2)-1;
        xramp3 = groupix(trial)*2;
    
        % Plot start distance where x = ramp number for each time bin
        scatter(xramp1,r1startavg,35,Color(clist{groupix(trial)}),"filled")
        scatter(xramp3,r3startavg,35,Color(clist{groupix(trial)}),"filled")
        % Plot line connecting ramps in the same trial
        plot([xramp1,xramp3],[r1startavg,r3startavg],'color',Color(clist{groupix(trial)}),'LineWidth',2)
    
    % Plot distance from food for each time bin within each ramp (right)
    subplot(r,c,2)
        hold on
        % Establish dynamic x axis tick locations
        xramp1 = groupix(trial);
        xramp3 = groupix(trial) + max(groupix) + 1;
        % Offset each data point slightly so they don't overlap
        buff = 0.3; %ex: x = 2, buffer = 1.7,2.3
        nincr = 50;
        poss = linspace(-buff,buff,nincr); % possible locations = linearspace(min range, max range, how possible locations within that space)
        idx = randi(nincr,1); % random number between 1 and nincrement
        b = poss(idx); % position of randomly selected number within the possible locations we've determined
        % b to be added to actual x value to create jitter
    
        % Plot start distance where x = time bin for each ramp
        scatter(xramp1+b,r1startavg,35,Color(clist{groupix(trial)}),"filled")
        scatter(xramp3+b,r3startavg,35,Color(clist{groupix(trial)}),"filled")
    
        % Identify x tick locations of the same time bin in each ramp section
        xtickloc = [xtickloc;xramp1;xramp3];
end

% Dynamically establish x axis limits
buff = 0.5;
xmax = max(unique(groupix))*2 + buff;
xmin = min(unique(groupix))*2 - 1 - buff;

% Format figure
formatFig(fig,bkgrd_color,[r,c]);

% Formatting for distance for each ramp within time bins
subplot(r,c,1)
    % Set x axis tick labels
    set(gca,'XTick',1:(length(unique(groupix))*3),'XTickLabel',repmat({'1','3'},[1,length(unique(groupix))])) %TODO: make tick range dynamic
    % Set x limits
    xlim([xmin,xmax])
    % Create axes labels
    xlabel('Ramp')
    ylabel('Distance from food (mm)')

    % Make time bin names for x axis labels
    dataString = [];
    dummy = unique(groupix);
    for idx = 1:length(dummy)
        i = dummy(idx);
        % bincenter = mean(binedges(i:i+1));
        % binwidth = binedges(i+1)-bincenter;
        dataString{idx} = [num2str(binedges(i)) '-' num2str(binedges(i+1))];
    end
    xticklabels = repmat(dataString,[1,2]);

% Formatting for distance for each time bin within ramps
subplot(r,c,2)
    % Set x axis tick labels
    set(gca,'XTick',unique(xtickloc),'XTickLabel',xticklabels,'XTickLabelRotation',30) %TODO: make tick range dynamic
    % Create axes labels
    xlabel('Start time (hr)')
    ylabel('Distance from food (mm)')



% STATS
% Create varibles for significance (true or false) and p value
h = [];
p = [];

% Compare between ramps for each time bin
for i = 1:ngroups %TODO: make into variable 
    % Pull data from stats matrix for each group into a new variable
    z = stats(i).data;
    % Run a t-test comparing ramps
    if isempty(z)
        h(i) = false;
        p(i) = nan;
    else
        [h(i),p(i)] = ttest(z(:,1),z(:,2));
        % Display comparisons that are significant
        if h(i)
            disp(['Group ' num2str(i) ' p = ' num2str(p(i))])
        end
    end
end

% Create variables to hold data for each ramp
[r1,r3] = deal([]);

% Compare between time bins for each ramp
for i = 1:ngroups
    % Pull data from stats matrix for each group into a new variable
    z = stats(i).data;
    % Load stats data into matrices (false = add data to beginning of matrix, not end)
    if ~isempty(z)
        r1 = autoCat(r1,z(:,1),false);
        r3 = autoCat(r3,z(:,2),false);
    end
end

% Run an ANOVA comparing time bins
[p1,tbl] = anova1(r1,unique(groupix),'off');
[p2,tbl2] = anova1(r3,unique(groupix),'off');

% Add stats to figure
figure(fig)
subplot(r,c,1)
    y = rangeLine(fig,2,true);
    
    % Plot average lines
    for i = 1:ngroups
        % Get x axis locations
        xramp1 = (i*2)-1;
        xramp3 = i*2;
        
        % Calculate average change in distance for each ramp
        avgchange = mean(stats(i).data);

        % Plot average change in distance between each ramp
        plot([xramp1,xramp3],avgchange,'color',foreColor,'LineWidth',4,'Marker','o')
        
        % If change is signficiant, plot asterisk
        if h(i) % h is a logical, h(i) = 1 = true
            scatter((xramp1+xramp3)/2,y,150,foreColor,"filled","*","MarkerEdgeColor",foreColor)
        end
    end

% Save figure
save_figure(fig, [figdirectory 'Starting distance to food before ', inputVar], '-png');


%% FIGURE: Distance to food at start of each ramp across Zeitgeber time

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(bkgrd_color); 

% Time offset from start of experiment to start of each ramp
toffset = [100/60, 212/60, 324/60];
% Period of time before/after ramp to look at distance from food
windowsize = 5; %in minutes

% FIGURE
% Create figure
fig = getfig;
hold on
for trial = 1:ntrials
    hold on
    % Establish y variable
    loc = data(trial).data.foodwell;
    y = data(trial).data.dist2well(:,loc);
    
    % Calculate number of frames within the established windowsize
    tempPoints = getTempTurnPoints(temp_protocol);
    nframes = tempPoints.fps*60*windowsize;
    
    % Plot distance to food at start of each ramp across Zeitgeber time
    for r = 1:3
        % Identify frames in the window at the start of each ramp
        rstart = tempPoints.hold(r,2)-nframes:tempPoints.hold(r,2);
        % Calulate average distance to food at each ramp start
        rstartavg = mean(y(rstart));
        % Calculate Zeitgeber time (in hrs) of each ramp start
        zeit = expstarttime_hr(trial) - 8 + toffset(r); 
        % Create scatter plot
        scatter(zeit,rstartavg,30,Color(clist{groupix(trial)}),'filled')
    end
end

% Format figure
formatFig(fig,bkgrd_color);
% Create axis labels
xlabel('Zeitgeber time')
ylabel('Distance to food (mm)')

% Save figure
save_figure(fig, [figdirectory 'Distance to food across Zeitgeber time'], '-png');

%% FOR LATER -- ANALYSIS: Calculate flies within the outer ring of the region
clearvars('-except',initial_vars{:})

R = data(1).data(1).data.r; % radius of the arena in pixels
innerR = R/sqrt(2); % radius of the inner 50% occupancy space
dist_from_edge = (R - innerR)/pix2mm;

% Find the percent of the flies that are in the outer ring
for exp = 1:num.exp
    counts = []; ring_per = [];
    for trial = 1:num.trial(exp)
        x = data(exp).data(trial).data.x_loc; % x locations for the entire experiment
        y = data(exp).data(trial).data.y_loc; % x locations for the entire experiment
        centre = data(exp).data(trial).data.centre; %distance to center of arena 
        D = sqrt(((x-centre(1)).^2 + (y-centre(2)).^2)); %distance from center of arena
        loc = D<=R & D>=innerR; % find the locations that are between edge and inner R
        ringCount = sum(loc,2);
        counts = autoCat(counts, ringCount, false); %count #flies in the outer ring
        ring_per = autoCat(ring_per,(ringCount./data(exp).T.NumFlies(trial)).*100,false); % convert to percent & combine
    end
    % pool the data
    grouped(exp).ring.all = counts;
    grouped(exp).ring.percent = ring_per; 
    grouped(exp).ring.avg = mean(ring_per,2,'omitnan');
end

% Pull the data together: 
for exp = 1:num.exp
    temps = grouped(exp).position.temp_list; % pre-binned temperatures
    nTemp = length(temps);
    rates = grouped(exp).position.temp_rates; % temperature rates in this experimental group
    cIdx = find(rates<0); %cooling index
    hIdx = find(rates>0); %heating index
    locs = grouped(exp).position.loc;
    [raw_c, raw_h] = deal(nan(nTemp,num.trial(exp))); %empty raw structures to fill in for each exp
    all_ring = [];
    
    % Update the averages for the classic temperature bins 
    for t = 1:nTemp
        % cooling frames for this temp
        c_frames = locs(cIdx,t).frames;
        h_frames = locs(hIdx,t).frames;
        if all(isnan(c_frames)) || all(isnan(h_frames))
            continue
        end
        raw_c(t,:) = mean(grouped(exp).ring.percent(c_frames,:),1,'omitnan');
        raw_h(t,:) = mean(grouped(exp).ring.percent(h_frames,:),1,'omitnan');
    end

    % find the avg and err and save to group structure
    grouped(exp).ring.increasing.raw = raw_h;
    grouped(exp).ring.increasing.avg = mean(raw_h, 2, 'omitnan');
    grouped(exp).ring.increasing.err = std(raw_h, 0, 2, 'omitnan');
    grouped(exp).ring.decreasing.raw = raw_c;
    grouped(exp).ring.decreasing.avg = mean(raw_c, 2, 'omitnan');
    grouped(exp).ring.decreasing.err = std(raw_c, 0, 2, 'omitnan');
    grouped(exp).ring.temps = temps;
end


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
    disp(T.ExperimentID{trial})
end

disp('sleep data loaded')

%% FIGURE: Percentage of flies sleeping over time

clearvars('-except',initial_vars{:}) %clear all variables except initial ones
plot_err = true; %plot error region
[foreColor,backColor] = formattingColors(bkgrd_color); % get background colors

% set up figure aligments
r = 4; %rows
c = 1; %columns
sb(1).idx = [1]; %temp timecourse
sb(2).idx = [2,3,4]; %distance from food timecourse %TODO: normalize this to something more intuitive?

LW = 1.5; %linewidth
sSpan = 180; %smoothing function

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

% clist = {'LemonChiffon','Gold','Tomato','DodgerBlue'};

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
        % y_err =  g(i).std;
       plot(x,y,'color',Color(clist{i}),'LineWidth',0.75)
       % h = plot_error_fills(plot_err, x, y, y_err, Color(clist{i}),'-png', 0.1);
end
    
formatFig(fig,bkgrd_color,[r,c],sb);

subplot(r,c,sb(1).idx)
    set(gca,'xcolor',backColor)
    ylabel('Temperature (\circC)')

 subplot(r,c,sb(2).idx) 
    legend(dataString,'Location','northwest','color',backColor,'box','off','textColor',foreColor)
    xlabel('Time (min)')
    ylabel('Percentage of flies sleeping')
    ylim([0,20])
    set(gca, 'ytick', 0:5:20);

% Save figure
save_figure(fig, [figdirectory 'Percentage of flies sleeping'], '-png');

% TODO
% compare sleep across hold periods
% distance to food while sleeping (histogram)

%% FIGURE: Average sleep per trial by start time

clearvars('-except',initial_vars{:}) 
plot_err = true; 
[foreColor,backColor] = formattingColors(bkgrd_color); 

plotdata = nan([ntrials,2]); %create empty variable of a specific size

% Create figure
fig = getfig('Average sleep per trial by start time',true);
    hold on

for trial = 1:ntrials
    y = data(trial).sleep.sleepNum;
    y = (y/T.NumFlies(trial))*100;
    sleepavg = mean(y);
    plotdata(trial,2) = sleepavg;
    plotdata(trial,1) = expstarttime_hr(trial);

    % FIGURE
    scatter(plotdata(trial,1),plotdata(trial,2),50,Color(clist{groupix(trial)}),"filled")
    % Create axes labels
    xlabel('Experiment start time (hr)')
    ylabel('Percentage of flies sleeping')

end

formatFig(fig,bkgrd_color);

% Save figure
save_figure(fig, [figdirectory 'Average sleep per trial by start time'], '-png');

%% FIX -- FIGURE AND STATS: average experiment speed
clearvars('-except',initial_vars{:})
autoSave = false;
plot_err = true;
[foreColor,backColor] = formattingColors(bkgrd_color);
LW = 1.5;
r = 1;
c = 3;
sb(1).idx = 1:2; %avg speed per temp
sb(2).idx = 3;   %avg speed

dataString = [];
dummy = unique(groupix);

% make group names for figure legend and labels
for idx = 1:length(dummy)
    i = dummy(idx);
    % bincenter = mean(binedges(i:i+1));
    % binwidth = binedges(i+1)-bincenter;
    dataString{idx} = [num2str(binedges(i)) '-' num2str(binedges(i+1))];
end

g = []; %empty variable

matsize = [];
for trial = 1:ntrials
    y = data(trial).speed.avg;
    matsize(trial,:) = size(y);
end

mt = nan([max(matsize,[],'all')],ntrials);

%create matrices with dimensions of num trials and max size of y
for i = 1:ngroups
    g(i).name = [num2str(binedges(i)) '-' num2str(binedges(i+1))]; 
    [g(i).raw,g(i).temp] = deal(mt); %empty variable
    g(i).count = 1;
end

%add data into matrices and add onto group index count
for trial = 1:ntrials
    y = data(trial).speed.avg;
    x = data(trial).occupancy.temp;
    idx = groupix(trial); %group identity for each trial (before or after noon)
    yend = size(y,1);
    g(idx).raw(1:yend,g(idx).count) = y;
    g(idx).temp(1:yend,g(idx).count) = x;
    g(idx).count = g(idx).count + 1;
end

%remove empty columns in each group
for idx = 1:ngroups
    g(idx).raw(:,g(idx).count:end) = [];
    g(idx).temp(:,g(idx).count:end) = [];
end

%plot data and error
for i = 1:ngroups 
    g(i).speedavg = mean(g(i).raw,2,'omitnan'); 
    g(i).std = std(g(i).raw,0,2,'omitnan'); %error
    g(i).tempavg = mean(g(i).temp,2,'omitnan');
end

fig = getfig('Average speed',true);
    for i = 1:ngroups
        if g(i).count==1
            continue
        end
        subplot(r,c,sb(1).idx)
            hold on
                x = g(i).tempavg;
                y = g(i).speedavg;
                y_err =  g(i).std;
            plot(x,y,'color',Color(clist{i}))
            plot(x, y-y_err,'color',Color(clist{i}))
            plot(x, y+y_err,'color',Color(clist{i}))
            % h = plot_error_fills(plot_err, x, y, y_err, Color(clist{i}),'-pdf', 0.2);
        % subplot(r,c,sb(2).idx)
        %     hold on 
        %         x = groupix(trial);
        %         y = g(i).speedavg;
        %         % y = mean(y);
        %         buff = 0.3;
        %         nincr = 50;
        %         poss = linspace(-buff,buff,nincr);
        %         idx = randi(nincr,1); %random number between 1 and nincrement
        %         b = poss(idx);
        %     scatter(x,y,35,Color(clist{groupix(trial)}),"filled")      
    end
    
formatFig(fig,bkgrd_color,[r,c],sb);

subplot(r,c,sb(1).idx)
    legend(dataString,'Location','northwest','color',backColor,'box','off','textColor',foreColor)
    xlabel('Temperature (\circC)')
    ylabel('Speed (mm/s)')
    xlim([17,25])

 subplot(r,c,sb(2).idx) 
    set(gca,'xcolor',backColor)
    ylabel('Average speed (mm/s')

% TODO: change x axis to be time and add temperature over time as a separate graph. 
% Problem is we're plotting multiple y values for each x temperature value



%% FIGURE: Speed across start times with temperature log

clearvars('-except',initial_vars{:})
plot_err = true; 
[foreColor,backColor] = formattingColors(bkgrd_color); 

% Set up figure aligments
r = 5; 
c = 3; 
sb(1).idx = [1,2]; % temp timecourse
sb(2).idx = [4,5,7,8,10,11,13,14]; % speed timecourse
sb(3).idx = [3,6,9,12,15];

LW = 1.5; % linewidth
sSpan = 180; % smoothing function

% Make group names for figure legend and labels
dataString = [];
dummy = unique(groupix);
for idx = 1:length(dummy)
    i = dummy(idx);
    % bincenter = mean(binedges(i:i+1));
    % binwidth = binedges(i+1)-bincenter;
    dataString{idx} = [num2str(binedges(i)) '-' num2str(binedges(i+1))];
end

% Create empty variable to hold data
g = []; 

% Set up matrices for data
for i = 1:ngroups
    % Give each group a name based on time bins
    g(i).name = [num2str(binedges(i)) '-' num2str(binedges(i+1))];
    % Set up empty matrices for each group
    [g(i).raw,g(i).temp,g(i).time] = deal([]);
end

% Pull and organize data
for trial = 1:ntrials
    % Establish plot variables
    y = data(trial).speed.avg;
    x = data(trial).occupancy.time;
    z = data(trial).occupancy.temp;
    % Identify the group index for each trial (which time bin it falls into)
    idx = groupix(trial);
    % Load data into each matrix within variable g
    g(idx).raw = autoCat(g(idx).raw,y,false);
    g(idx).time = autoCat(g(idx).time,x,false);    
    g(idx).temp = autoCat(g(idx).temp,z,false);
end

% Calculate average distance, time, and temp for each group, and standard deviation for distance
for i = 1:ngroups
    g(i).speedavg = mean(g(i).raw,2,'omitnan'); 
    g(i).timeavg = mean(g(i).time,2,'omitnan');
    g(i).tempavg = mean(g(i).temp,2,'omitnan');
    g(i).std = std(g(i).raw,0,2,'omitnan'); 
end

% FIGURE
% Create figure and establish name, screen position for viewing, and size
fig = getfig('average speed over time',true);
    for i = 1:ngroups
        % Plot temperature timecourse
        subplot(r,c,sb(1).idx)
            hold on
            plot(g(i).timeavg,g(i).tempavg,'color',Color(clist{i}),'LineWidth',LW)
        % Plot distance to food timecourse
        subplot(r,c,sb(2).idx)
            hold on 
                x = g(i).timeavg;
                % y = g(i).speedavg;
                y = smooth(g(i).speedavg,sSpan, 'moving');
                % y_err = g(i).std;
            plot(x,y,'color',Color(clist{i}),'LineWidth',0.75)
            % Display error area
            % h = plot_error_fills(plot_err, x, y, y_err, Color(clist{i}),'-png', 0.2);
        % subplot(r,c,sb(3).idx)
    end

% Format figure (call fig variable first, it's usually just called fig)
formatFig(fig,bkgrd_color,[r,c],sb);

% Formatting for temperature timecourse
subplot(r,c,sb(1).idx)
    % Make x axis same color as background (gca = get current axis)
    set(gca,'xcolor',backColor)
    % Create y axis label
    ylabel('Temperature (\circC)')

% Formatting for distance to food timecourse
 subplot(r,c,sb(2).idx) 
    % Create legend
    legend(dataString,'Location','northwest','color',backColor,'box','off','textColor',foreColor)
    % Create axes labels
    xlabel('Time (min)')
    ylabel('Speed (mm/s)')

% Save figure
save_figure(fig, [figdirectory 'Speed over time'], '-png');

        
        
        
        
 