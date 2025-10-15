

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
        
    end
end

%% TODO: group trends across the flies to see how they trend
clearvars('-except',initial_var{:})

r = 1;
c = 2;

% Giant histogram across the two different fly sexes
fig = getfig('',1); 

% histogram
subplot(r,c,1); hold on
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

% scatter plot of the 95% cumulative distribution time for food visits --
% could ask if this changes (and the mean) for each temp regime...
subplot(r,c,2); hold on
for sex = 1:2
    plotData = [];
    for trial = 1:num.trials
        plotData = [plotData; FV(trial, sex).duration];
    end
    max(plotData)
    histogram(plotData,'FaceColor',data.color(sex,:),'BinEdges',0:2:120,'FaceAlpha',0.6)
end

formatFig(fig, blkbgd);


% Cumulative visit duration distribution -- find the 95% number and mean
% for each of the flies in the trial

% does the 95% number or the mean change across temperature regimes?

[f,x] = ecdf(plotData);


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

%% TODO: compare the behavior between males and females

%% TODO: what is the inter-eating-bout frequency? How does this change or
% not change across temperature?











