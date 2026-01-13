

% create plots about the order of behaviors / priority of behaviors

%% Time delay to each behavior
clearvars('-except',initial_var{:})
[foreColor, backColor] = formattingColors(blkbgd); % get background colors

% from each transition point (temp regime type) how long does it take for
% each behavior to emerge?  [female flies won't have courtship as a
% possibility]


% find the time ROIs for each transition point ... 

% currently hard coded for the LTS data...
% transition points:
% WT: warming & above 25
% WS: cooling & above 25
% CT: cooling & below 25
% CS: warming & below 25

switch groupName
    case 'Berlin LTS caviar'
        % order: WT, WS, CT, CS, WT, WS
        transitions = [data.warming_idx(1,1);...
                   data.cooling_idx(1,1);...
                   data.cooling_idx(1,1) + round((data.warming_idx(2,1)-data.cooling_idx(1,1))/2);...
                   data.warming_idx(2,1);...
                   data.warming_idx(2,1) + round((data.cooling_idx(2,1)-data.warming_idx(2,1))/2);...
                   data.cooling_idx(2,1);...
                   data.cooling_idx(2,2)];
        roi_safe = [0,1,0,1,0,1]; % safe = true, unsafe = false;
        trans_cat = {'WT', 'WS', 'CT', 'CS', 'WT', 'WS'}; % names for the categories
        nTrans = length(transitions)-1;
end


              
% quick demo figure to show what the selected regions are: 
fig = getfig; hold on
    plot(data.temp,'color', foreColor,'linewidth', 2);
    v_line(transitions)
    formatFig(fig, blkbgd);
    ylabel('temperature (\circC)')
    set(gca,'xcolor', 'none')
save_figure(fig, [figDir 'temp regions'],fig_type);


% -----------------------------------------------------------------------------------
% Make a demo logical structure that includes the behavior heirarchy
% of sleep>courtship>all such that a fly which is sleeping would
% behaviorally take precident over a fly being in the escape ring since it
% cannot be actively escaping while sleeping (aka the non-location based
% behaviors take front by default) 

% -----------------------------------------------------------------------------------
% ** 'data' structure has all the data for each behavior while 'exc_data'
% ** is the behaviorally exclusive structure where there can only be one
% ** 'behavior' at a time 
% -----------------------------------------------------------------------------------

exc_data = struct; % exclusive data
fields = {'wing_ext', 'wing_ext_all', 'court_chase', 'chase_all', 'circling_all',...
    'circling_1sec', 'CI', 'FlyOnFood', 'OutterRing', 'foodQuad', 'sleep'};
% load the new structure with the behavioral state data as is -- no
% substitutions
for i = 1:length(fields)
    exc_data.(fields{i}) = data.(fields{i});
end
% courtship index all is any of the other non-time dependent courtship
% behaviors combined
exc_data.CI_all = replaceNaN(exc_data.chase_all,false) |...
                  replaceNaN(exc_data.circling_all,false) |...
                  replaceNaN(exc_data.wing_ext_all,false);


% set order priority first for courtship for the male fly only (female fly
% doesn't have courtship behavior recorded in this data set) 
fields = {'FlyOnFood', 'OutterRing', 'foodQuad'};
for i = 1:length(fields)
    temp =  squeeze(exc_data.(fields{i})(:,M,:));
    % set the location where there is courtship to false for the other
    % behaviors 
    loc = replaceNaN(exc_data.CI,false); % 
    loc = logical(loc);
    temp(loc) = false;
    temp = double(temp);
    % put the newly edited data back into the data structure
    exc_data.(fields{i})(:,M,:) = temp;
end

% set order priority first for sleep for male courtship lists (since there
% isn't female courtship behaviors)
fields = {'wing_ext', 'wing_ext_all', 'court_chase', 'chase_all',...
              'circling_all','circling_1sec', 'CI','CI_all'};
for i = 1:length(fields)
    temp =  (exc_data.(fields{i}));
    % set the locations of courtship to zero when the male fly is asleep
    loc = squeeze(exc_data.sleep(:,M,:));
    temp(loc) = false;
    temp = double(temp);
    exc_data.(fields{i}) = temp;
end

% set order priority first for sleep for M and F data types
fields = {'FlyOnFood', 'OutterRing', 'foodQuad'};
for i = 1:length(fields)
    for sex = 1:2
        temp =  squeeze(exc_data.(fields{i})(:,sex,:));
        % set all other behaviors to zero when the fly is sleeping
        loc = squeeze(exc_data.sleep(:,sex,:));
        temp(loc) = false;
        exc_data.(fields{i})(:,sex,:) = temp;
    end
end
% -----------------------------------------------------------------------------------

% when does each behavior happen within each region? 
% 'behavior_onset' structure shows the first instance of a behavior within
% each of the different temperature regimes

d_fields = { 'FlyOnFood', 'OutterRing', 'foodQuad', 'sleep'}; % fields for all the flies
s_fields = {'wing_ext', 'wing_ext_all', 'court_chase', 'chase_all', 'circling_all',...
    'circling_1sec', 'CI','CI_all'}; % only fields for the male fly

% Pull out the event onset times for each of the behaviors for each of the
% temperature regimes
behavior_onset = struct; % initialize as empty structure

% double fields
for i = 1:length(d_fields) % for each different behavior
    for sex = 1:2 % for each of the flies
        field_name = d_fields{i}; % this is the behavior that is being analyzed this time
        raw = squeeze(exc_data.(field_name)(:,sex,:)); % extract the data for the choice behavior
        % set an empty matrix that will be filled with first behavior timings
        behavior_onset(sex).(field_name) = nan([num.trials,nTrans]); 
        for roi = 1:nTrans % loop all temperature regions
            locs = raw(transitions(roi):transitions(roi+1),:); % all locations within the ROI that the behavior is seen
            for f = 1:num.trials % find the first event time in this regime for each fly (if there are any)
                idx = find(locs(:,f));
                if ~isempty(idx)
                    behavior_onset(sex).(field_name)(f,roi) = idx(1);
                end
            end
        end
    end
end

% single (Male only) fields (aka courtship)
sex = M;
for i = 1:length(s_fields)
    field_name = s_fields{i};
    raw = (exc_data.(field_name));
    behavior_onset(sex).(field_name) = nan([num.trials,nTrans]);
    for roi = 1:nTrans % loop all temperature regions
        locs = raw(transitions(roi):transitions(roi+1),:); % all locations within the ROI that the behavior is seen
        for f = 1:num.trials % find the first event time in this regime for each fly (if there are any)
            idx = find(locs(:,f));
            if ~isempty(idx)
                behavior_onset(sex).(field_name)(f,roi) = idx(1);
            end
        end
    end
end

% ------------------------------ DATA VISUALIZATION -----------------------------
% scatter plot of the onset timing within each of the regions
switch groupName
    case 'Berlin LTS caviar'
        r = 2; 
        c = 4;
end


%  FLY BEHAVIOR SEQUENCE 
fields = {'foodQuad', 'FlyOnFood', 'OutterRing', 'sleep', 'CI'};
kolors = {'gold', 'red', 'purple', 'dodgerblue', 'green'}; % colors for the diff behaviors

sz = 75; 
y = 1:num.trials;
sexList = {'male', 'female'};
for sex = 1:2
    fig = getfig; 
    for i = 1 : nTrans
        subplot(r,c,i)
        hold on
        % time to onset vs row = fly
        for f = 1:length(fields)
            % skip the courtship plot points for the female fly
            if sex == F && strcmp(fields{f},'CI')
                continue
            end
            x = behavior_onset(sex).(fields{f})(:,i); % length of time to behavior start
            x = (x/30)/60; % convert to minutes
            scatter(x,y,sz, Color(kolors{f}),'filled', 'square')
        end
        title(trans_cat{i},'Color',foreColor)
    end
    % formatting
    formatFig(fig, blkbgd,[r,c]);
    matchAxis(fig, true);
    for i = 1:nTrans
        subplot(r,c,i)
        set(gca, 'ycolor', 'none')
    end
    for i = 1:nTrans
        subplot(r,c,i)
        xlabel('onset (min)')
        % xlim([0, 10])
    end
    save_figure(fig, [figDir sexList{sex} ' behavior_onset per region'],fig_type);
end

%% FIGURE: time of first behavior zoom on specific temp region and time
fields = {'foodQuad', 'FlyOnFood', 'OutterRing', 'sleep', 'CI'};
kolors = {'gold', 'red', 'purple', 'dodgerblue', 'green'}; % colors for the diff behaviors

sz = 75; 
y = 1:num.trials * 2; % include both male and female flies stacked
sexList = {'male', 'female'};

% trans_cat
i = 4; % temperature regime
zoom_roi = [0 10]; %time of zoom 
% zoom_roi = [0 70]; %time of zoom 

fig_folder = createFolder([figDir 'behavior onset zoom in/']);

fig = getfig('Time of first behavior', 1, [715 710]); 
hold on
% time to onset vs row = fly
for f = 1:length(fields)
    y = 1:num.trials * 2; % include both male and female flies stacked
    x = behavior_onset(M).(fields{f})(:,i); % length of time to behavior start
    % skip the courtship plot points for the female fly
    if strcmp(fields{f},'CI')
        y = 1:num.trials; % male only data
    else % this field exists for the female fly, so plot it
        x = [x; behavior_onset(F).(fields{f})(:,i)];
    end
    x = (x/30)/60; % convert to minutes
    scatter(x,y,sz, Color(kolors{f}),'filled', 'square')
end
% title(trans_cat{i},'Color',foreColor)

% formatting
formatFig(fig, blkbgd);
set(gca, 'ycolor', 'none','FontSize', 25)
xlabel('onset (min)')
xlim(zoom_roi)
% test lines between fly rows
% h_line(0.5:(2*num.trials)+0.5, 'grey', ':',0.5)
h_line(num.trials+0.5, 'grey', '--',1)


fig_str = [trans_cat{i} ' time ' num2str(zoom_roi(1)) ' to '...
           num2str(zoom_roi(2)) ' min'];
save_figure(fig, [fig_folder fig_str],fig_type);


%% Figure: Plot out all the time delays for each behavior by type not by fly
% do this as well with both the male and female flies

fields = {'FlyOnFood', 'OutterRing', 'sleep', 'CI'};
c = length(fields);
r = 1;

fig = getfig; 
    for f = 1:length(fields)
        subplot(r,c,f)
        fieldName = fields{f};
        hold on
        
        for sex = 1:2
            switch sex
                case 1 % male
                    offset = -0.2;                 
                    kolor = Color('dodgerblue');
                case 2 % female
                    offset = 0.2;
                    if strcmp(fieldName, 'CI')
                        continue
                    end
                    kolor = Color('pink');
            end
            x = repmat(1+offset:nTrans+offset, [num.trials,1]);
            y = behavior_onset(sex).(fieldName);
            y = y./(30*60); % convert from frames to minutes
            
            scatter(x, y, 25, kolor, 'filled','XJitter', 'density',...
                'XJitterWidth',0.2,'MarkerFaceAlpha',0.7)
            scatter(1+offset:nTrans+offset, mean(y,1,'omitnan'),100, kolor, "filled")
            errorbar(1+offset:nTrans+offset, mean(y,1,'omitnan'),std(y,0,1,'omitnan')./(sqrt(num.trials)),'color', ...
                kolor,'linewidth', 1.5,'linestyle', 'none')
        end
        set(gca, 'xtick', 1:nTrans, 'xticklabel', trans_cat)
        xlim([0,nTrans+1])
        title(fieldName)
    end

formatFig(fig, blkbgd,[r,c]);
matchAxis(fig,true);
for i = 2:c
    subplot(r,c,i)
    set(gca, 'ycolor', 'none')
end
subplot(r,c,1)
ylabel('time to first event (min)')

save_figure(fig, [figDir 'behavior_onset per sex scatter'],fig_type);


%% Figure : time delays to each behavior, combined across sex and multiple reps of the same temp regime
% currently only set up for this temp protocol
if ~strcmp(groupName,'Berlin LTS caviar')
    return
end

fields = {'FlyOnFood', 'OutterRing', 'sleep', 'CI'};
c = length(fields);
r = 1;

temp_regimes = unique(trans_cat);
kolor = Color('dodgerblue');
offset = 0.1;
FA = 0.4; % scatter plot face alpha level

fig = getfig; 
for f = 1:length(fields)
    subplot(r,c,f)
    fieldName = fields{f};
    hold on       
    for i = 1:length(temp_regimes)
        loc = find(strcmp(trans_cat,temp_regimes{i}));
        % pull the data for all regions with the same type
        plotData = [];
        for sex = 1:2
            if strcmp(fieldName,'CI') & sex==2 %skip courtship for female flies
                continue
            end
            for t = 1:length(loc) % for each of the different temp regimes
                plotData = [plotData; behavior_onset(sex).(fieldName)(:,loc(t))];
            end
        end
        % plot the data for this temp regime: 
        x = i*ones(size(plotData));
        y = plotData/(30*60); % convert to minutes
        scatter(x-offset, y, 35, foreColor, 'filled','XJitter', 'density',...
            'XJitterWidth',0.2,'MarkerFaceAlpha',FA)
        scatter(i+offset, mean(y,1,'omitnan'),100, kolor, "filled")
        errorbar(i+offset, mean(y,1,'omitnan'),std(y,0,1,'omitnan')./(sqrt(num.trials)),'color', ...
            kolor,'linewidth', 1.5,'linestyle', 'none')
    end
    % Labels
    set(gca, 'xtick', 1:length(temp_regimes),'XTickLabel',temp_regimes)
    xlim([0.5,length(temp_regimes)+0.5])
    title(fieldName)
end
formatFig(fig, blkbgd, [r,c]);
for f = 2:length(fields)
    subplot(r,c,f)
    set(gca, 'ycolor', 'none')
end
matchAxis(fig, true);
subplot(r,c,1)
ylabel('Time to first behavior (min)')

save_figure(fig, [figDir 'behavior_onset scatter'],fig_type);

%% FIGURE : time delays to each behavior, combined across sex for each temperature regime
% currently only set up for this temp protocol

if ~strcmp(groupName,'Berlin LTS caviar')
    return
end

fields = {'FlyOnFood', 'OutterRing', 'sleep', 'CI'};
field_names = {'Food', 'Escape', 'Sleep', 'Courtship'};
c = length(fields);
r = 1;

temp_regimes = unique(trans_cat);
kolor = Color('dodgerblue');
offset = 0.1;
FA = 0.4; % scatter plot face alpha level

fig = getfig('time to first behavior', 1, [ 991 900]); 
for i = 1:length(temp_regimes) % subplot for each of the temp regime types
    subplot(r,c,i)

    for f = 1:length(fields)
        fieldName = fields{f};
        hold on       
    
        loc = find(strcmp(trans_cat,temp_regimes{i}));
        % pull the data for all regions with the same type
        plotData = [];
        for sex = 1:2
            if strcmp(fieldName,'CI') & sex==2 %skip courtship for female flies
                continue
            end
            for t = 1:length(loc) % for each of the different temp regimes
                plotData = [plotData; behavior_onset(sex).(fieldName)(:,loc(t))];
            end
        end
        % plot the data for this temp regime: 
        x = f*ones(size(plotData));
        y = plotData/(30*60); % convert to minutes
        scatter(x-offset, y, 35, foreColor, 'filled','XJitter', 'density',...
            'XJitterWidth',0.2,'MarkerFaceAlpha',FA)
        scatter(f+offset, mean(y,1,'omitnan'),100, kolor, "filled")
        errorbar(f+offset, mean(y,1,'omitnan'),std(y,0,1,'omitnan')./(sqrt(num.trials)),'color', ...
            kolor,'linewidth', 1.5,'linestyle', 'none')
    end
    % Labels
    set(gca, 'xtick', 1:length(fields),'XTickLabel',field_names)
    xlim([0.5,length(fields)+0.5])
    title(temp_regimes{i})
end
formatFig(fig, blkbgd, [r,c]);
for f = 2:length(temp_regimes)
    subplot(r,c,f)
    set(gca, 'ycolor', 'none')
end
matchAxis(fig, true);
subplot(r,c,1)
ylabel('Time to first behavior (min)')

save_figure(fig, [figDir 'behavior_onset scatter by temp regime'],fig_type);



%% TODO 12/30 Do more flies eat before sleep in each temp regime?
% Figure: percent of flies that eat before sleeping vs percent of flies
% that sleep before eating within each temperature regime
% Do more flies eat before sleeping than not? 





% follow up:  what is the time delay between sleeping and eating or eating
% then sleeping?

% how to do this: for each regime look at the trials that have both eating
% and sleeping within the regime and then see which comes first -- add to
% the tally for that behavior type

% for each temp regime -- what fraction of flies go to food before they
% sleep and in what region do they sleep and what is the time delay between
% food visit and sleep
% another way: what percent of food visits are followed by sleep across the
% different temp regimes?

% for each regime: fraction of flies that go to food before they sleep:
r = 2;
c = 4;
sz = 50;

fig = getfig;
for t = 1:nTrans
    subplot(r,c,t)
    food_time = [behavior_onset(M).FlyOnFood(:,t); behavior_onset(F).FlyOnFood(:,t)];
    sleep_time = [behavior_onset(M).sleep(:,t); behavior_onset(F).sleep(:,t)]; 
    onsetDiff = sleep_time-food_time;
    nflies = length(food_time); % number of flies in the set
    % is there sleep in this regime for any of the flies?
    conditions = {'S>F', 'F>S', 'S', 'F', '--'}; 
    perc = []; % fraction of flies that are in each category (perc short for percent) 
    perc(1) = sum(onsetDiff<0)/nflies; % no food before sleep
    perc(2) = sum(onsetDiff>0)/nflies; % food before sleep
    perc(3) = sum(isnan(food_time) & ~isnan(sleep_time))/nflies; % flies that didn't go to food but slept
    perc(4) = sum(isnan(sleep_time) & ~isnan(food_time))/nflies; % flies that didn't sleep but went to food
    perc(5) = sum(isnan(food_time) & isnan(sleep_time))/nflies; % flies that didn't sleep or go to food
    
    % scatter(1:5, perc*100,sz,foreColor,'filled')
    switch roi_safe(t)
        case true
            kolor = Color('mediumspringgreen');
        case false
            kolor = Color('deeppink');
    end
    bar(1:5, perc*100,'faceColor', kolor)
    set(gca, 'xtick', 1:length(perc), 'XTickLabel', conditions)
    title(trans_cat{t}, 'color', foreColor)
    ylabel('% flies')
end
% formatting
formatFig(fig, blkbgd, [r c]);
matchAxis(fig, true);
   
save_figure(fig, [figDir 'food vs sleep order for all temp regions'],fig_type);
        
% same idea but collapse across the different temperature regimes and safe
% zones to just 'thermal threat' and 'safe' regions
r = 1;
c = 2;
fig = getfig('threat and safety summary',1,[1072 669]);
for type = 1:2 % safe and threat
    subplot(r,c,type)
    % extract the locations of the safe vs threat temp regions
    switch type
        case 1 % safe
            kolor = Color('mediumspringgreen');
            t_idx = find(roi_safe==true);
            title_str = 'safe';
        case 2 % threatening
            kolor = Color('deeppink');
            t_idx = find(roi_safe==false);
            title_str = 'threat';
    end
    [food_time, sleep_time] = deal([]); %initialize empty structures
    for t = 1:length(t_idx)
        food_time = [food_time; behavior_onset(M).FlyOnFood(:,t_idx(t));...
                     behavior_onset(F).FlyOnFood(:,t_idx(t))];
        sleep_time = [sleep_time; behavior_onset(M).sleep(:,t_idx(t));...
                      behavior_onset(F).sleep(:,t_idx(t))]; 
    end
    onsetDiff = sleep_time-food_time;
    nflies = length(food_time); % number of flies in the set
    % is there sleep in this regime for any of the flies?
    conditions = {'S>F', 'F>S', 'S', 'F', '--'}; 
    perc = []; % fraction of flies that are in each category 
    perc(1) = sum(onsetDiff<0)/nflies; % no food before sleep
    perc(2) = sum(onsetDiff>0)/nflies; % food before sleep
    perc(3) = sum(isnan(food_time) & ~isnan(sleep_time))/nflies; % flies that didn't go to food but slept
    perc(4) = sum(isnan(sleep_time) & ~isnan(food_time))/nflies; % flies that didn't sleep but went to food
    perc(5) = sum(isnan(food_time) & isnan(sleep_time))/nflies; % flies that didn't sleep or go to food

    bar(1:length(perc), perc*100,'faceColor', kolor,'FaceAlpha',1)
    set(gca, 'xtick', 1:length(perc), 'XTickLabel', conditions)
    title(title_str, 'color', foreColor)
    ylabel('% flies')
end
% formatting
formatFig(fig, blkbgd, [r c]);
matchAxis(fig, true);

save_figure(fig, [figDir 'food vs sleep order collapsed safe and threat regions'],fig_type);

%% Zoom in on specific behavior period to show the order of behaviors:


% [TODO? find the time start and stops of each behavior in sequence...this
% minimizes the number of symbols that need to be plotted ]

% create an 'image' that is the size of each frame by the number of trials (flies)
% beh_timing = struct;
%  % for f = 1:length(fields)
%         %     % find the start and stop of each conscutive behavior
%         %     raw = exc_data.(fields{f})(:,trial);
%         %     start_loc = find([false; (diff(raw)==1)]); % behavior onset
%         %     stop_loc = find([false; (diff(raw)==-1)]); % behavior offset
%         %     % save timing indexes to a new data structure
%         %     beh_timing(trial,sex).(fields{f}) = [start_loc,stop_loc];
%         % end


%% TODO: Plot out the sequence of behaviors color coded and then sort them
% by the the most of one behavior, then the next ext. so that it shows as a
% continuum 

% ------------------------------ DATA VISUALIZATION -----------------------------

% Make the figure as an imagesc where each 'pixel' is a discrete
% point in time and it color coded by the behavior that is being
% displayed... would still want to chunk these and order them by the
% durations of each behavior eventually TODO

% scatter plot of the onset timing within each of the regions
switch groupName
    case 'Berlin LTS caviar'
        r = 2; 
        c = 4;
end

% BEHAVIOR SEQUENCE 
fields = {'OutterRing','foodQuad', 'FlyOnFood', 'CI', 'sleep'};
% set up a color code for each behavior
kolors = {'purple','yellow','red', 'lime', 'dodgerblue'}; % colors for the diff behaviors
cmap = backColor;
for i = 1:length(kolors)
    cmap = [cmap; Color(kolors{i})];
end

% fill each pixel with the behavioral state of the fly
real_img = []; % this will be filled with the combined fly data later
for sex = 2:-1:1 % reverse direction so female is on top
    raw_image = nan(size(exc_data.CI)); 
    for trial = 1:num.trials
        for f = 1:length(fields) % each behavior gets coded with a number
            if sex==F && strcmp(fields{f}, 'CI') % skip CI for female fly
                continue
            end
            raw = exc_data.(fields{f});
            if size(raw,2) == 2 % double data
                loc = exc_data.(fields{f})(:,sex,trial); % frames with this behavior
            else
                loc = exc_data.(fields{f})(:,trial);
            end
            loc = logical(replaceNaN(loc, false));
            raw_image(loc,trial) = f;
        end
    end
    raw_image = replaceNaN(raw_image, 0);
    real_img = [real_img, raw_image]; % add the data from this fly sex
end

sub_dir = createFolder([figDir, 'behavior states/']);
for i = 1:nTrans
    fig = getfig;

    imagesc(real_img')
    title(trans_cat{i},'color', foreColor,'FontSize',12)
    % label the plot
    set(gca, 'xtick', [],'ytick', [])
    xlabel('time','FontSize',12)
    ylabel('fly','FontSize',12)
    C = colorbar;
    colormap(cmap) % set the image colors by behavior
    % set the color bar labels to the appropriate behaviors
    num_increments = (size(cmap,1)*2)+1; % how many increments to divide colorbar
    inc = linspace(C.Limits(1),C.Limits(2),num_increments);
    C.Ticks = inc(2:2:end-1);
    C.TickLabels = ['inner arena', fields];
    C.Color = foreColor;

    % zoom into different transition sections using xlims
    xlim([transitions(i), transitions(i+1)])
    
    set(fig, 'color', backColor)
    
    fig_str = [sub_dir 'all behavior states ' num2str(i) ' ' trans_cat{i}];
    saveas(fig, [fig_str '.png'])
    close(fig)

    % save_figure(fig, fig_str,fig_type);
end









% how much of each behavior is there time wise for each temp regime?
totals = nan([length(fields),nTrans,num.trials,sex]);
for sex = 1:2
  for f = 1:length(fields) %for each behavior
    raw = exc_data.(fields{f}); 
    for t = 1:nTrans % for each temp regime
        if sex==2 & strcmp(fields{f},'CI')
            continue
        end
        if size(raw,2)==2 % for each sex
            totals(f,t,:,sex) = squeeze([sum(raw(transitions(t):transitions(t+1),sex,:))]);
        else % courtship related behaviors
            totals(f,t,:,sex) = squeeze([sum(raw(transitions(t):transitions(t+1),:))]);
        end
    end
  end
end

% figure showing relative amounts of each behavior or location for each
% temp regime
for sex = 1:2
    fig = getfig;
    
    for t = 1:nTrans
        subplot(r,c,t); hold on
        % plot the raw time sums for each behavior in a given temp regime
        raw = squeeze(totals(:,t,:,sex))./(30*60); % add conversion to minutes
        x = repmat((1:length(fields))',[1,num.trials]);
        scatter(x,raw,35,foreColor,"filled",'MarkerFaceAlpha', 0.6,...
            'XJitter', 'density','XJitterWidth', 0.3)
        % plot the avg across flies
        x_avg = 1:length(fields);
        y_avg = mean(raw,2,'omitnan');
        scatter(x_avg, y_avg, 75, cmap(2:end,:), 'filled')
        title(trans_cat{t}) 
    end
    formatFig(fig, blkbgd, [r,c]);
    for t = 1:nTrans
         subplot(r,c,t); hold on
        % set(gca, 'xcolor', 'none')
        set(gca, 'xtick', 1:length(fields),'xticklabel', fields)
        ylabel('duration of behavior (min)')
    end
    matchAxis(fig, true);
        
    save_figure(fig, [figDir sexList{sex} ' behavior duration per region'],fig_type);

end


%% TODO: sex based differences in the onset of different behavior times? 


%% TODO: where is the majority of courtship taking place? 
% Plot the locations for all instances of courtship to see if there is a
% spatial bias

% Location regions of interest: outer ring, food quadrant, inner quadrants

allPoints = false(size(data.OutterRing)); % initialize empty structure 

% remove outer ring from the food quadrant points: 
ring = data.OutterRing; 
foodQuad = data.foodQuad;  
foodQuad(data.OutterRing) = false; % TODO: update this to full logical matrix

data.OutterRing




%% TODO: How often do flies go from each behavior to the next? Make this a chord chart? 
% pull out the chunks of each behavior (like start and stop points for each
% 'time' of a behavior


