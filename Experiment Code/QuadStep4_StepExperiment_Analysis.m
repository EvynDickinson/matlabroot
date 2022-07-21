% This script looks at data files that had no food and used the step ramp
% temperature protocols

%% Initial ANALYSIS: General data organization (forward and backwards compatability
if ~exist('pix2mm','var')
    pix2mm = 12.8;
end
% fill out the data structure
for trial = 1:ntrials
    if isempty(data(trial).occupancy)
        data(trial).occupancy = data(trial).data.occupancy;
    end
    if isempty(data(trial).nflies)
        data(trial).nflies = data(trial).data.nflies;
    end
    if isempty(data(trial).wellLabels)
        data(trial).wellLabels = data(trial).data.wellLabels;
    end
    if isempty(data(trial).dist2wells)
        try
            data(trial).dist2wells = data(trial).data.occupancy.dist2wells;
        catch
            for well = 1:4
                data(trial).dist2wells(:,well) = data(trial).occupancy.dist2wells(well).N(:,1)./pix2mm;
            end
        end
    elseif isfield(data(trial).dist2wells, 'N')
        temp = data(trial).dist2wells;
        data(trial).dist2wells = [];
        for well = 1:4
            data(trial).dist2wells(:,well) = temp(well).N(:,1)./pix2mm;
        end
    end
end

%% ANALYSIS: Average speed for binned temperatures
sSpan = 180;
binSpace = 1;
binbuffer = 0.5;

[tempAll, plotData] = deal([]); 
% Average speed for each temperature
for trial = 1:ntrials
    % pull data for use:
    raw = data(trial).speed.avg;
    temp = (data(trial).occupancy.temp);
    % ROI of useable data:
    tPoints = getTempTurnPoints(T.TempProtocol{trial});
    ROI = [tPoints.DownROI';tPoints.UpROI';tPoints.HoldROI']; 
    tempdata = temp(ROI);
    speeddata = raw(ROI);

    % determine temperature bins
    tLow = tPoints.threshLow-binbuffer;
    tHigh = tPoints.threshHigh+binbuffer;    
    t_roi = tLow:binSpace:tHigh; 
%     t_roi(1) = tLow; t_roi(end) = tHigh;
    % sort data by temperature:
    binID = discretize(tempdata,t_roi);
    for tt = 1:length(t_roi)-1
        loc = binID==tt;
        plotData(trial,tt) = mean(speeddata(loc),1,'omitnan');
        tempAll(trial,tt) = mean(t_roi(tt:tt+1));
    end
end

speed.avg = plotData;
speed.temp = tempAll;

initial_vars{end+1} = 'speed';
clearvars('-except',initial_vars{:})

%% FIGURE: plot avg speed for each temperature

fig = figure;
hold on
for trial = 1:ntrials
    x = speed.temp(trial,:);
    y = speed.avg(trial,:);
    x(isnan(y)) = [];
    y(isnan(y)) = [];
    
    plot(x,y,'Marker','o','LineWidth',1,'color','w')
end
xlabel('Temperature (\circC)')
ylabel('Avg Speed (mm/s)')
formatFig(fig,true);

save_figure(fig, [figDir ExpGroup ' avg speed vs temperature'], '-png');

clearvars('-except',initial_vars{:})

%% FIGURE: Compare speed for each step to it's control over time

sSpan = 30;
% PULL RAW DATA:
plotData = struct;
for trial = 1:ntrials
    tPoints = getTempTurnPoints(T.TempProtocol{trial});
    ROI = tPoints.hold;
    % pull the speed for each region of interest
    for tt = 1:size(tPoints.hold,1)/2 % for each of the steps in the temp protocol
        % get appropriate time roi
        con = (tt*2)-1; %control
        tes = tt*2; %test temp
        % pull speed data from select time roi
        plotData(tt).control(:,trial) = smooth(data(trial).speed.avg(ROI(con,1):ROI(con,2)),sSpan,'moving');
        plotData(tt).test(:,trial) = smooth(data(trial).speed.avg(ROI(tes,1):ROI(tes,2)),sSpan,'moving');
        % pull temp data for select roi
        plotData(tt).cTemp(trial) = mean(data(trial).occupancy.temp(ROI(con,1):ROI(con,2)),'omitnan');
        plotData(tt).tTemp(trial) = mean(data(trial).occupancy.temp(ROI(tes,1):ROI(tes,2)),'omitnan');
    end
end
    
% PROCESS THE DATA INTO AVERAGES AND TRIALS
nsteps = 4;

fig = figure;
for tt = 1:nsteps %todo  -- change this if more than 4 test temps are used
    subplot(2,2,tt)
    hold on
    plot(plotData(tt).control,'color','w')
    plot(plotData(tt).test,'color',Color('royalblue'))
    str = [num2str(round(mean(plotData(tt).cTemp))) '\circC vs ' num2str(round(mean(plotData(tt).tTemp))) '\circC'];
    title(str)
    ylabel('mm/s')
    xlabel('time')
end

formatFig(fig,true,[2,2]);
save_figure(fig, [figDir ExpGroup ' speed control vs temp drop timecourse'], '-png');

% FIGURE: scatter of avg during each temp period
fig = figure; set(fig, 'pos', [478 339 520 489]); hold on
buff = 0.20;
LW = 1.5;
for tt = 1:nsteps
    for trial = 1:ntrials
        x = [tt-buff,tt+buff];
        y = [mean(plotData(tt).control(:,trial)), mean(plotData(tt).test(:,trial))]; 
        plot(x,y,'color', 'w', 'Marker','o','LineWidth',LW) 
        scatter(tt+buff,y(2),30,Color('royalblue'),'filled')
        scatter(tt-buff,y(1),30,Color('red'),'filled')
    end
    % plot average bars
    cMean = mean(mean(plotData(tt).control));
    tMean = mean(mean(plotData(tt).test));
    plot([tt-(buff*2),tt+(buff*2)],[cMean,cMean],'Color', Color('red'),'LineWidth',LW*2)
    plot([tt-(buff*2),tt+(buff*2)],[tMean,tMean],'Color', Color('royalblue'),'LineWidth',LW*2)
end
ylabel('Speed (mm/s)')
xlabel('Temp Step')
ax = gca;
% set(ax,'XTick',1:nsteps,'XTickLabel',{'23>16','23>14','23>12','23>10'})
set(ax,'XTick',1:nsteps,'XTickLabel',{'23>10','23>12','23>14','23>16'})

formatFig(fig,true);
save_figure(fig, [figDir ExpGroup ' control vs temp drop speed avg'], '-png');

clearvars('-except',initial_vars{:})

%% FIGURE: speed timecourse across all trials
sSpan = 180;
row = 4;
col = 1;
sb(1).idx = 1;
sb(2).idx = 2:4;


fig = figure; set(fig,'pos', [2198 252 622 620])
subplot(row,col,sb(1).idx)
plot(data(1).occupancy.time,data(1).occupancy.temp,'color','w','LineWidth',2)
ylabel('\circC')
axis tight

subplot(row,col,sb(2).idx)
hold on
for trial = 1:ntrials
    x = data(trial).occupancy.time;
    y = smooth(data(trial).speed.avg,sSpan,'moving');
    plot(x,y,'LineWidth',0.5,'Color','w')
end
xlabel('time (min)')
ylabel('speed (mm/s)')
axis tight

formatFig(fig, true,[row,col],sb);
save_figure(fig, [figDir ExpGroup ' speed timecourse'], '-png');

%% FIGURES RELATED TO WORK

% -------------------------
% FIGURE: Work and temp timecourse
sSpan = 180;
fig = figure; set(fig,'pos',[557 383 752 516]); hold on
for trial = 1:ntrials
    y = smooth(data(trial).data.T.tempWork,sSpan,'moving');
    z = smooth(data(trial).occupancy.temp,sSpan, 'moving');
    x = data(trial).occupancy.time;
    yyaxis left
    plot(x,y,'color',Color('white'),'linewidth',1,'LineStyle','-','Marker','none')
    yyaxis right
    plot(x,z,'color',Color('royalblue'),'linewidth',1,'LineStyle','-','Marker','none')
end
formatFig(fig,true);
%labels
xlabel('Time (min)')
axis tight
yyaxis left
ylabel('Temp controller work (%)')
set(gca,'YColor','w')
yyaxis right
ylabel('Temp (\circC)')
set(gca,'YColor',Color('royalblue'))

save_figure(fig, [figDir ExpGroup ' work and temp vs speed timecourse'], '-png');
% -------------------------


% -------------------------
% FIGURE: work histogram
x = [];
for trial = 1:ntrials
    x = [x; data(trial).data.T.tempWork]; % work
end

fig = figure; set(fig,'Position',[2169 502 434 457])
h = histogram(x);
h.EdgeColor = "none";
h.FaceColor = Color('royalblue');
h.FaceAlpha = 1;
xlabel('Temp Work (%)')
ylabel('Count')
formatFig(fig,true);
save_figure(fig, [figDir ExpGroup ' work histogram'], '-png');

% -------------------------


% -------------------------
% FIGURE: work and temp correlation
[x,y] = deal([]);
for trial = 1:ntrials
    x = [x; data(trial).data.T.tempWork]; % work
    y = [y; data(trial).occupancy.temp];  % temp
end
z = corrcoef(x,y);
plotData = [x,y];
plotData = sortrows(plotData,1,"ascend");

fig = figure; hold on
    scatter(x,y,15,Color('grey'),'filled')
    plot(plotData(:,1),smooth(plotData(:,2),180,'moving'),'color','w','linewidth',2)
    %labels
    xlabel('Work (%)')
    ylabel('Temp (\circC)')
    title(['Correlation Coefficient: ' num2str(z(1,2))],'color', 'w')
formatFig(fig, true);
axis tight
save_figure(fig, [figDir ExpGroup ' temp - work correlation'], '-png');


% -------------------------
% FIGURE: Scatter plot of speed for all work points
fig = figure; set(fig,'pos',[680 530 521 448]); hold on
for trial = 1:ntrials
    x = data(trial).data.T.tempWork;
    y = data(trial).speed.avg;

    scatter(x,y,15,Color('grey'))
end
xlabel('Temp controller work (%)')
ylabel('Speed (mm/s)')

formatFig(fig,true);
save_figure(fig, [figDir ExpGroup ' work vs speed scatter'], '-png');
% -------------------------


% -------------------------
% FIGURE: avg and err speed for each work bin
[x,y,z,dataIn] = deal([]);
for trial = 1:ntrials
    dataIn = data(trial).data.T.tempWork;
    x = [x; dataIn]; % work
    z = [z; nan; diff(dataIn)]; % change in work
    y = [y; data(trial).speed.avg]; % speed
end
plotData = [];
w_roi = linspace(-100,100,50);
binID = discretize(x,w_roi);
for tt = 1:length(w_roi)-1
    loc = binID==tt;
    plotData(tt,1) = w_roi(tt);
    plotData(tt,2) = mean(y(loc),1,'omitnan');
    plotData(tt,3) = std(y(loc),0,'omitnan');
end
fig = figure; hold on
scatter(plotData(:,1),plotData(:,2),50,'w','filled')
for tt = 1:length(w_roi)-1
    plot([plotData(tt,1),plotData(tt,1)],...
         [plotData(tt,2)-plotData(tt,3),plotData(tt,2)+plotData(tt,3)],'color', 'w','LineWidth',1.5,'marker', '_')
end
xlabel('Temp controller work (%)')
ylabel('Speed (mm/s)')
formatFig(fig,true);
save_figure(fig, [figDir ExpGroup ' work vs speed avg'], '-png');
% -------------------------


% -------------------------
% FIGURE: avg speed for change in work (binned)
[x,y,dataIn] = deal([]);
for trial = 1:ntrials
    dataIn = data(trial).data.T.tempWork;
    x = [x; nan; diff(dataIn)]; % change in work
    y = [y; data(trial).speed.avg]; % speed
end
plotData = [];
[binID,binedges] = discretize(x,50);

for tt = 1:length(binedges)-1
    loc = binID==tt;
    plotData(tt,1) = binedges(tt);
    plotData(tt,2) = mean(y(loc),1,'omitnan');
    plotData(tt,3) = std(y(loc),0,'omitnan');
end
fig = figure; hold on
scatter(plotData(:,1),plotData(:,2),50,'w','filled')
for tt = 1:length(binedges)-1
    plot([plotData(tt,1),plotData(tt,1)],...
         [plotData(tt,2)-plotData(tt,3),plotData(tt,2)+plotData(tt,3)],'color', 'w','LineWidth',1.5,'marker', '_')
end
xlabel('Change in work (%)')
ylabel('Speed (mm/s)')
formatFig(fig,true);
axis tight
save_figure(fig, [figDir ExpGroup ' change in work vs speed avg'], '-png');
% -------------------------


clearvars('-except',initial_vars{:})

%% ANALYSIS & FIGURES: Look at Decay and Recovery from cold temp steps

clearvars('-except',initial_vars{:})

confidenceInterval = 0.25; % 10-percent
timeStep = 3; % in minutes
stepSize = 30; % test every ten seconds
timeFrames = timeStep*3*60; % convert minutes to number of frames
sSpan = 1; 

% sort and grab data
plotData = struct;
for trial = 1:ntrials
    tPoints = getTempTurnPoints(T.TempProtocol{trial});
    % TEST
    % control temp is the last 5 minute period in the neutral temp hold
    control = struct;
    nsteps = size(tPoints.hold,1)/2;
    for ii = 1:nsteps
        idx = (ii*2)-1;
        control(ii).roi = tPoints.hold(idx,2)-timeFrames:tPoints.hold(idx,2); % frames for control period
        control(ii).speed = smooth(data(trial).speed.avg(control(ii).roi),sSpan,'moving'); % avg speed for this trial in control period
        control(ii).avgSpeed = mean(control(ii).speed); % avg speed for the control period
    end
    % CONTROL
    test = struct;
    stepLength = min(tPoints.hold(:,2) - tPoints.hold(:,1)); % determine the test period  
    cropROIs = [(1:stepSize:stepLength-timeFrames)',(timeFrames:stepSize:stepLength)'];
    round(timeFrames/2):stepSize:stepLength-round(timeFrames/2); % center time points for each window
    for ii = 1:nsteps
        idx = (ii*2);
        roi = tPoints.hold(idx,1) : tPoints.hold(idx,2);
        test(ii).temp = num2str(round(mean(data(trial).occupancy.temp(roi))));
        x = data(trial).speed.avg(roi);
        for tt = 1:length(cropROIs)
            test(ii).speed(:,tt) = smooth(x(cropROIs(tt,1):cropROIs(tt,2)),sSpan,'moving');
            test(ii).roi(:,tt) = [cropROIs(tt,1); cropROIs(tt,2)] + roi(1) - 1;
            % run a set of t-tests to look for significant differences in speed from control to test:
            [h,p] = ttest2(control(ii).speed,test(ii).speed(:,tt),'Alpha',confidenceInterval);
            test(ii).h(tt) = h;
            test(ii).p(tt) = p;
        end
        test(ii).avgSpeed = mean(test(ii).speed,1,'omitnan');
    end
    plotData(trial).test = test;
    plotData(trial).control = control;
    plotData(trial).tPoints = tPoints;
end


% ----------------------------------------
% FIGURE: p-values and significance across samples
colorList = {'royalblue','orangered', 'yellow', 'plum'};
row = 2; col = 1;
fig = figure;  set(fig, 'pos', [2061 553 751 420])
legendStr = [];
for trial = 1:ntrials
    % p-values
    subplot(row,col,1); hold on
        for ii = 1:nsteps
            plot(plotData(trial).test(ii).p,'color', Color(colorList{ii}))
            if trial == 1
                legendStr{ii} = [plotData(trial).test(ii).temp '\circC'];
            end
        end
    % p-values
    subplot(row,col,2); hold on
        for ii = 1:nsteps
            plot(plotData(trial).test(ii).h,'color', Color(colorList{ii}))
        end
end
% Formatting: 
subplot(row,col,1);
    xlabel('time slot')
    ylabel('p-value')
    l = legend(legendStr);
subplot(row,col,2);
    xlabel('time slot')
    ylabel('significance')
    ylim([-0.1,1.1])
formatFig(fig,true,[row,col]);
set(l,'box','off','TextColor','w')

save_figure(fig, [figDir ExpGroup ' p-value and significance control vs test bins'], '-png');

% ----------------------------------------
% FIGURE: p-values and significance across samples
% colorList = {'royalblue','orangered', 'yellow', 'plum'};
[row, col] = subplot_numbers(nsteps,2);
fig = figure; set(fig, 'Position',[2049 553 776 420])
for ii = 1:nsteps
    for trial = 1:ntrials
        subplot(row,col,ii); hold on
        plot(plotData(trial).test(ii).p)
    end
    title([plotData(trial).test(ii).temp '\circC'])
    h_line(0.05,'red','--')
    ylabel('p-value')
    xlabel('time bin')
    set(gca,'xtick', [])
end
% Formatting: 
formatFig(fig,true,[row,col]);

save_figure(fig, [figDir ExpGroup ' p-value control vs test'], '-png');
% ----------------------------------------


% ----------------------------------------
sSpan = 1;
% FIGURE: significance across trials
colorList = {'royalblue','orangered', 'yellow', 'teal'};
% [row, col] = subplot_numbers(nsteps,2);
legendStr = [];
fig = figure; set(fig, 'pos', [2129 355 545 589])
for ii = 1:nsteps
%     subplot(row,col,ii); 
    hold on
    x = [];
    for trial = 1:ntrials
        x(:,trial) = plotData(trial).test(ii).h;
    end
    plot(smooth(mean(x,2),sSpan,'moving'),'color',Color(colorList{ii}),'linewidth',2)
    legendStr{ii} = [plotData(trial).test(ii).temp '\circC'];
    ylim([0.60,1.00])
end
% Formatting: 
formatFig(fig,true);
legend(legendStr,'box', 'off', 'textcolor', 'w','Location','southeast')
xlabel('time bin')
ylabel('avg significant difference from control')

save_figure(fig, [figDir ExpGroup ' difference from control significance'], '-png');
% ----------------------------------------


%% FIGURE: Plot distribution of 5 min bins compared to control (random time bin selected)

% Fig demo: 
LW = 1.5;
sSpan = 6;
binedges = 0:0.5:20;
ndemos = 3;
row = nsteps;
col = 3;
ttList = randi([1,size(test(ii).speed,2)],1,ndemos);

sb = [];
ind = 1;
fig = figure; set(fig,'pos', [1956 343 964 687])
for ii = 1:nsteps
% SPEED COURSE:
    sb((ii*2)-1).idx = [ind,ind+1];
    subplot(row,col,sb((ii*2)-1).idx)
    hold on
    for tt = ttList  % linspace(1,size(test(ii).speed,2),ndemos)
        plot(smooth(test(ii).speed(:,tt),sSpan,'moving'),'linewidth',LW)
    end
    plot(smooth(control(ii).speed,sSpan,'moving'),'color','w','linewidth',LW)
    axis tight
    ylabel('speed')
    title([test(ii).temp '\circC'])
    if ii == nsteps
        xlabel(['time (total ' num2str(timeStep) ' min)'])
    end
% HISTOGRAM
    ind = ind+1;
    sb(ii*2).idx = ind+1;
    subplot(row,col,sb(ii*2).idx)
    hold on
    legStr = {};
    for tt = ttList
        histogram(test(ii).speed(:,tt),binedges)
        % legend:
        if test(ii).h(tt)
            % TODO: update the string to have significance 
            legStr{end+1} = [num2str(tt) ' diff'];
        else
            legStr{end+1} = [num2str(tt)];
        end
    end
    % control histogram
    histogram(control(ii).speed,binedges,'FaceColor','w','FaceAlpha',0.6)
    xlim([0,8])
    ylabel('count')
    xlabel('speed')
    legend(legStr,'box', 'off', 'textcolor', 'w')
% Update subplot index
    ind = ind+2;
end
formatFig(fig,true,[row,col],sb);
for ii = 1:nsteps
    subplot(row, col, sb((ii*2)-1).idx)
    set(gca,'XTick',[])
    subplot(row, col, sb(ii*2).idx)
    set(gca,'XTick',[])
end

save_figure(fig, [figDir ExpGroup ' sample binned roi test vs control'], '-png');


% ----------- :FIGURE: ------------
CList = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

row = 3;
col = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;

trial = 16;
fig = figure; set(fig, 'pos', [2064 364 468 550])
subplot(row,col,sb(1).idx)
    plot(data(trial).occupancy.time, data(trial).occupancy.temp,'color','w','linewidth', 2)
    ylabel('(\circC)')
    axis tight
subplot(row,col,sb(2).idx)
    % plot the full timecourse in grey...
    y = data(trial).speed.avg;
    x = data(trial).occupancy.time;
    plot(x,y,'color',Color('grey'),'linewidth', 1)
    % plot the demo sections in color
    hold on
    for ii = 1:nsteps
        idx = 0;
        time = x(plotData(trial).control(ii).roi);
        plot(time, plotData(trial).control(ii).speed,'color', 'w', 'LineWidth',LW)
        for tt = ttList
            time = x(plotData(trial).test(ii).roi(1,tt):plotData(trial).test(ii).roi(2,tt));
            idx = idx+1;
            plot(time, plotData(trial).test(ii).speed(:,tt),'LineWidth',LW,'color',CList(idx,:))
        end
    end
    axis tight
    ylabel('speed (mm/s)')
    xlabel('time (min)')
ylim([0,20])
formatFig(fig, true,[row,col],sb);
subplot(row,col,sb(1).idx)
set(gca,'xtick',[],'xColor','k')

save_figure(fig, [figDir ExpGroup ' full timecourse sample roi test vs control'], '-png');



%% 








%% 



    



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    