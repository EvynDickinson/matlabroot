% Load raw data (XX.raw) then run the analysis section of QuadStep3_1
% ANALYSIS: General data organization (forward and backwards compatability

%% FIGURE: Alignment of temperature, speed, and work across full experiment timecourse

fig = figure; set(fig, 'position', [2030 476 824 534]);
sSpan = 60;
for trial = 1:ntrials
    % pull colors
    switch T.TempProtocol{trial}
        case 'frequency_sweeps_18'
            kolor = Color('royalblue');
        case 'frequency_sweeps_12'
            kolor = Color('Gold');
    end
    X = data(trial).occupancy.time;
    Y = data(trial).occupancy.temp;
    speed = data(trial).speed.avg;
    work = data(trial).data.T.tempWork;
    
    % Temperature Plot
    subplot(3,1,1); hold on
    plot(X, Y, 'linewidth', 1,'color',kolor)
    
    
    % Speed Plot
    subplot(3,1,2); hold on
    plot(X, smooth(speed,sSpan,'moving'), 'linewidth', 1,'color',kolor)
    
    % Speed Plot
    subplot(3,1,3); hold on
    plot(X, work, 'linewidth', 1,'color',kolor)
    
end
% LABELS | FORMATTING:
formatFig(fig, true,[3,1]);

% temp
subplot(3,1,1)
ylabel('temp (\circC)')
set(gca,'XColor','k')
% speed
subplot(3,1,2)
ylabel('speed (mm/s)')
set(gca,'XColor','k')
% work
subplot(3,1,3)
ylabel('work')
xlabel('time (min)')


save_figure(fig, [figDir 'temperature work speed alignment'], '-png');

clearvars('-except',initial_vars{:})
fprintf('Next\n')


%% FIGURE: 18C speed and work organized by frequency: 
clearvars('-except',initial_vars{:})

% transition points between the different frequencies in the time series
frequencyBorders = [6452; 68621; 99761; 115359]; 

% sort data into frequency groups
nfreq = length(frequencyBorders)-1;
plotData = struct;
for ii = 1:nfreq
   for trial = 1:ntrials %only the first experiment type
       ROI = frequencyBorders(ii):frequencyBorders(ii+1);
       plotData(ii).speed(:,trial) = data(trial).speed.avg(ROI);
       plotData(ii).work(:,trial) = data(trial).data.T.tempWork(ROI);
       plotData(ii).temp(:,trial) = data(trial).occupancy.temp(ROI);
       plotData(ii).time(:,trial) = data(trial).occupancy.time(ROI);
   end
end

% --------------------------------------------------------------------------------------------------

switch questdlg('Analyze data from 18C or 12C?','','18C', '12C', 'cancel','18C')
    case '18C'
        tempTag = '18C';
        trialROI = 1:4;
        tempColor = 'purple';
    case '12C'
        tempTag = '12C';
        trialROI = 5:8;
        tempColor = 'royalblue';
    case 'cancel'
        return
end

switch questdlg('Analyze speed or work?','','speed', 'work', 'cancel','speed')
    case 'speed'
        yTag = 'speed';
        yLab = 'speed (mm/s)';
    case 'work'
        yTag = 'work';
        yLab = 'work (%)';
    case 'cancel'
        return
end
   

sSpan = 60;
kolors = {'gold', 'darkorange', 'red'};
row = 3; col = 2; 
sb(1).idx = [1,2];
sb(2).idx = [3,5];
sb(3).idx = [4,6];
% binedges = linspace(0, 10, 11);

fig = figure; set(fig,'pos',[-1050 362 1045 770])
    % timecourse
    subplot(row,col,sb(1).idx); hold on
    for ii = 1:nfreq
        x = mean(plotData(ii).time(:,trialROI),2);
        y = smooth(mean(plotData(ii).(yTag)(:,trialROI),2),sSpan,'moving');
        z = smooth(mean(plotData(ii).temp(:,trialROI),2),sSpan,'moving');
        yyaxis left
        plot(x,z,'color',Color(tempColor),'linewidth',1,'linestyle', '-')
        yyaxis right
        plot(x,y,'color',Color(kolors{ii}),'linewidth',1.5,'linestyle', '-')
    end
    xlabel('time (min)')
    yyaxis right
    ylabel(yLab)
    yyaxis left
    ylabel('temp (\circC)')
    
    % speed histogram
    subplot(row,col,sb(2).idx); hold on
    for ii = 1:nfreq
        histogram(plotData(ii).(yTag)(:,trialROI),'facecolor', Color(kolors{ii}),'edgecolor','none')
    end
    xlabel(yLab) 
    ylabel('count')

    % speed PDF
    subplot(row,col,sb(3).idx); hold on
    for ii = 1:nfreq
        x = reshape(plotData(ii).(yTag)(:,trialROI),[1,numel(plotData(ii).(yTag)(:,trialROI))]);
        h = cdfplot(x);
        set(h,'color', Color(kolors{ii}),'linewidth', 1.5)
    end
    set(gca, 'xscale', 'linear')
    axis tight
    xlabel(yLab)

    % Format
    formatFig(fig, true,[row,col],sb);
    subplot(row,col,sb(1).idx)
    yyaxis right
    set(gca,'ycolor','white')
    yyaxis left
    set(gca,'ycolor',Color(tempColor))
    
save_figure(fig, [figDir yTag ' distribution temp ' tempTag], '-png');

%% FIGURES: how do temp and work relate across the three frequencies?
clearvars('-except',initial_vars{:})

% transition points between the different frequencies in the time series
frequencyBorders = [6452; 68621; 99761; 115359]; 

% sort data into frequency groups
nfreq = length(frequencyBorders)-1;
plotData = struct;
for ii = 1:nfreq
   for trial = 1:ntrials %only the first experiment type
       ROI = frequencyBorders(ii):frequencyBorders(ii+1);
       plotData(ii).speed(:,trial) = data(trial).speed.avg(ROI);
       plotData(ii).work(:,trial) = data(trial).data.T.tempWork(ROI);
       plotData(ii).temp(:,trial) = data(trial).occupancy.temp(ROI);
       plotData(ii).time(:,trial) = data(trial).occupancy.time(ROI);
   end
end


% --------------------------------------------------------------------------------------------------
% How much time is spent at each temp across the three frequencies?

sSpan = 360;
row = 2; col = 2; 
kolors = {'gold', 'darkorange', 'red'};

fig = figure; set(fig,'pos',[-1050 362 1045 770])
    for idx = 1:2
        if idx==1 
            trialROI = 1:4;
            tempTag = 'temp mean: 18\circC';
        else
            trialROI = 5:8;
            tempTag = 'temp mean: 12\circC';
        end
        % TEMP HISTOGRAM
        subplot(row,col,idx); 
        hold on
        for ii = 1:nfreq
            % prep data
            x = plotData(ii).temp(:,trialROI);
            x = x(:);
            pd = fitdist(x,'normal');
            y = pdf(pd,x);
            % plot data
            yyaxis left
            histogram(plotData(ii).temp(:,trialROI),'facecolor', Color(kolors{ii}),'edgecolor','none','FaceAlpha',0.5)
            yyaxis right
            plot(x,y,'color',Color(kolors{ii}),'linewidth',2,'linestyle', '-')
        end
        title(tempTag)
        xlabel('temperature (\circC)')
        % TEMP VS WORK
        subplot(row,col,idx+2); 
        hold on
        for ii = 1:nfreq
            y = plotData(ii).temp(:,trialROI); 
            x = plotData(ii).work(:,trialROI); 
            z = [x(:),y(:)];
            z = sortrows(z,1);
            plot(z(:,1), smooth(z(:,2),sSpan,'moving'),'color',Color(kolors{ii}),'linewidth',2,'linestyle', '-')
        end
        xlabel('work (%)')
        ylabel('temperature (\circC)')
    end

    % Format
    formatFig(fig, true,[row,col]);
    for idx = 1:2
        subplot(row,col,idx)
        yyaxis right
        set(gca,'ycolor','white')
        ylabel('CDF')
        yyaxis left
        set(gca,'ycolor','white')
        ylabel('count')
    end
    
save_figure(fig, [figDir 'temp vs work'], '-png');










