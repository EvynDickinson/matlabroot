

initial_var{end+1} = 'FV';

%% Extract the periods for food visits
clearvars('-except',initial_var{:})

FV = struct; % food visits (this will later be added to the DATA structure)
fps =  fly(M).fps;

% how many frames can be skipped before a sustained fly on food period is ended
frameDropAllowance = ceil((1/3) * fps);  % 1/3 of a second

% find instances of the flies on the food
% test for a single fly first: 
for trial = 1:num.trials
    for sex = 1:2
    
        % find periods of sustained time on the food aka, when there are long
        % periods of adjactent frames with flies on the food 
        frames = find(data.FlyOnFood(:,sex,trial));
        onFood = diff(frames)<=frameDropAllowance; % allow for a small gap in frames (dropped etc) [logical about 'frames']
        
        idx = find(onFood==0); % locations where a running list of flies on food frames exceed the max skip allowance
        frame_loc_stop = [frames(idx); frames(end)];
        frame_loc_start = [frames(1); frames(idx+1)];
        
        onfoodROI = [frame_loc_start, frame_loc_stop]; % frame indexes of when the fly started on the food and left the food
        nOnFood = size(onfoodROI,1); % how many times is the fly on the food total in the experiment
        onFoodDuration = (diff(onfoodROI,1,2))/fps; % duration of time (s) fly spent on the food

        % save data into the FV struct
        FV(trial,sex).ROI = onfoodROI;
        FV(trial,sex).nROI = nOnFood;
        FV(trial,sex).duration = onFoodDuration;

        % WORKING HERE **** subdivide the protocol into the four different types of temp regimes:
        % temps approaching 25 from each direction **** 

        % subdivide by heating/cooling temperature regime first (then subdivide later)
        % if start of food visit is during a temp regime it counts 
        % probably could add a buffer to the turn points later, but this will
        % suffice for now
        temp_regimes = {'cooling', 'warming', 'hold'};
        for tt = 1:3 % cooling, warming, holds
            type_str = temp_regimes{tt};
            % find which food visits start in each temp regime
            locs =  ismember(frame_loc_start, find(data.(type_str)));    
            % sort into groups for the temp regime types: 
            FV(trial,sex).(type_str).ROI = onfoodROI(locs,:);
            FV(trial,sex).(type_str).nROI = size(FV(trial,sex).(type_str).ROI,1);
            FV(trial,sex).(type_str).duration = onFoodDuration(locs);
        end



    end
end

%% TODO: group trends across the flies for the full experiment to see how they trend
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors


% Giant histogram across the two different fly sexes
fig = getfig('',1); hold on
% histogram
for sex = 1:2
    plotData = [];
    for trial = 1:num.trials
        plotData = [plotData; FV(trial, sex).duration];
    end
    max(plotData)
    histogram(plotData,'FaceColor',data.color(sex,:),'BinEdges',0:2:120,'FaceAlpha',0.6)
end
% formatting the figure
set(gca, 'YScale','log')
xlabel('food visit duration (s)')
ylabel('instances (#)')
formatFig(fig, blkbgd);

save_figure(fig, [figDir 'food visit duration histogram'],fig_type)

% Cumulative visit duration distribution -- find the 95% number and mean
% for each of the flies in the trial
% does the 95% number or the mean change across temperature regimes?
% scatter plot of the 95% cumulative distribution time for food visits --
% could ask if this changes (and the mean) for each temp regime...
r = 1;
c = 2;
offset = 0.3;
LW = 1;
sz = 50;
fig = getfig('',1);

for sex = 1:2
    kolor = data.color(sex,:);
    plotData = [];
    for trial = 1:num.trials
        y = FV(trial, sex).duration;
        [f,x] = ecdf(y);
        % find closest to 95%
        all_idx = find(f<=0.95);
        loc = all_idx(end); % max bin closest but less than 95%
        visit_dur = x(loc); % duration of visits that encapsulates 95% of visits
        % find mean visit duration
        mean_dur = mean(y,'omitnan');
        plotData = [plotData; visit_dur, mean_dur];
    end
    mean_dur = mean(plotData(:,2),'omitnan');
    mean_CDF = mean(plotData(:,1),'omitnan');
    x = shuffle_data(linspace(sex-offset+0.1, sex+offset-0.1, size(plotData,1)));
    
    % plot 95% cumulative dist. data
    subplot(r, c, 1)
    hold on
    scatter(x, plotData(:,1),sz,kolor,'filled')
    plot([sex-offset, sex+offset], [mean_CDF, mean_CDF],"Color",kolor, 'linewidth', LW)
    % plot avg visit duration data
    subplot(r, c, 2)
    hold on
    scatter(x, plotData(:,2),sz,kolor,'filled')
    plot([sex-offset, sex+offset], [mean_dur, mean_dur],"Color",kolor, 'linewidth', LW)
end
% formatting
formatFig(fig, blkbgd, [r c]);
subplot(r, c, 1)
    title('95% CDF','color', foreColor)
    ylabel('food visit duration (s) capturing 95%')
    set(gca, 'xcolor', 'none')
    xlim([0, 2.5])
subplot(r, c, 2)
    title('Mean','color', foreColor)
    ylabel('avg food visit duration (s)')
    set(gca, 'xcolor', 'none')
    xlim([0, 2.5])

save_figure(fig, [figDir 'food visit duration and CDF'],fig_type)

% check significance?
% quick t-test

%% Duration of food visits for each temperature regime...
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors


% Giant histogram across the two different fly sexes
fig = getfig('',1); hold on
% histogram
for sex = 1:2
    plotData = [];
    for trial = 1:num.trials
        plotData = [plotData; FV(trial, sex).duration];
    end
    max(plotData)
    histogram(plotData,'FaceColor',data.color(sex,:),'BinEdges',0:2:120,'FaceAlpha',0.6)
end
% formatting the figure
set(gca, 'YScale','log')
xlabel('food visit duration (s)')
ylabel('instances (#)')
formatFig(fig, blkbgd);

save_figure(fig, [figDir 'food visit duration histogram'],fig_type)

% Cumulative visit duration distribution -- find the 95% number and mean
% for each of the flies in the trial
% does the 95% number or the mean change across temperature regimes?
% scatter plot of the 95% cumulative distribution time for food visits --
% could ask if this changes (and the mean) for each temp regime...
r = 1;
c = 2;
offset = 0.3;
LW = 1;
sz = 50;
fig = getfig('',1);

for sex = 1:2
    kolor = data.color(sex,:);
    plotData = [];
    for trial = 1:num.trials
        y = FV(trial, sex).duration;
        [f,x] = ecdf(y);
        % find closest to 95%
        all_idx = find(f<=0.95);
        loc = all_idx(end); % max bin closest but less than 95%
        visit_dur = x(loc); % duration of visits that encapsulates 95% of visits
        % find mean visit duration
        mean_dur = mean(y,'omitnan');
        plotData = [plotData; visit_dur, mean_dur];
    end
    mean_dur = mean(plotData(:,2),'omitnan');
    mean_CDF = mean(plotData(:,1),'omitnan');
    x = shuffle_data(linspace(sex-offset+0.1, sex+offset-0.1, size(plotData,1)));
    
    % plot 95% cumulative dist. data
    subplot(r, c, 1)
    hold on
    scatter(x, plotData(:,1),sz,kolor,'filled')
    plot([sex-offset, sex+offset], [mean_CDF, mean_CDF],"Color",kolor, 'linewidth', LW)
    % plot avg visit duration data
    subplot(r, c, 2)
    hold on
    scatter(x, plotData(:,2),sz,kolor,'filled')
    plot([sex-offset, sex+offset], [mean_dur, mean_dur],"Color",kolor, 'linewidth', LW)
end
% formatting
formatFig(fig, blkbgd, [r c]);
subplot(r, c, 1)
    title('95% CDF','color', foreColor)
    ylabel('food visit duration (s) capturing 95%')
    set(gca, 'xcolor', 'none')
    xlim([0, 2.5])
subplot(r, c, 2)
    title('Mean','color', foreColor)
    ylabel('avg food visit duration (s)')
    set(gca, 'xcolor', 'none')
    xlim([0, 2.5])

save_figure(fig, [figDir 'food visit duration and CDF'],fig_type)

% check significance?
% quick t-test





%% TODO: Plot out the trajectories of the flies from before to after they visit the food
% define the time periods that want to be plotted before and after the food visit
clearvars('-except',initial_var{:})

fps =  fly(M).fps;

preD = 15; % pre period to plot in seconds
postD = preD; % post period to plot in seconds;
%convert to frames
preD = preD * fps; % 
postD = postD * fps; %


% figure; histogram(onFoodDuration/fps)


%% TODO: which of these are within certain temperature regimes? 

%% TODO: pull out the trajectory before and after the fly is on food to see
% what region they came from/go to


%% TODO: what is the inter-eating-bout frequency? How does this change or
% not change across temperature?












