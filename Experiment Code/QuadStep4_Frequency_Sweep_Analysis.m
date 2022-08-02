

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


%% Data sorting by frequency region: 

% transition points between the different frequencies in the time series
frequencyBorders = [6452; 68621; 99761; 115359]; 

% sort data into frequency groups
plotData = struct;
for ii = 1:length(frequencyBorders)-1
   for trial = 1:ntrials
       ROI = frequencyBorders(ii):frequencyBorders(ii+1);
       plotData(ii).speed(:,trial) = data(trial).speed.avg(ROI);
       plotData(ii).work(:,trial) = data(trial).data.T.tempWork(ROI);
       plotData(ii).temp(:,trial) = data(trial).occupancy.temp(ROI);
   end
end

% FIGURE: histogram of speed, work, temp for each frequency grouping





























