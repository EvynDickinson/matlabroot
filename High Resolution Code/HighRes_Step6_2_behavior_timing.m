

% create plots about the order of behaviors / priority of behaviors

%% Time delay to each behavior
clearvars('-except',initial_var{:})
foreColor = formattingColors(blkbgd); % get background colors

% from each transition point (temp regime type) how long does it take for
% each behavior to emerge?  [female flies won't have courtship as a
% possibility]

%
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
fig = figure; hold on
    plot(data.temp,'color', foreColor,'linewidth', 2);
    v_line(transitions)
    formatFig(fig, blkbgd);
    ylabel('temperature (\circC)')
    set(gca,'xcolor', 'none')
save_figure(fig, [figDir 'temp regions'],fig_type);


% -----------------------------------------------------------------------------------
% Make a demo logical structure that includes the behavior heirarchy
% of sleep>courtship>all such that a fly which is sleeing would
% behaviorally take precident over a fly being in the escape ring since it
% cannot be actively escaping while sleeping (aka the non-location based
% behaviors take front by default) 

exc_data = struct; % exclusive data
fields = {'wing_ext', 'wing_ext_all', 'court_chase', 'chase_all', 'circling_all',...
    'circling_1sec', 'CI', 'FlyOnFood', 'OutterRing', 'foodQuad', 'sleep'};
for i = 1:length(fields)
    exc_data.(fields{i}) = data.(fields{i});
end
exc_data.CI_all = replaceNaN(exc_data.chase_all,false) |...
                  replaceNaN(exc_data.circling_all,false) |...
                  replaceNaN(exc_data.wing_ext_all,false);


% set order priority first for courtship 
fields = {'FlyOnFood', 'OutterRing', 'foodQuad'};
for i = 1:length(fields)
    temp =  squeeze(exc_data.(fields{i})(:,M,:));
    loc = replaceNaN(exc_data.CI,false);
    loc = logical(loc);
    temp(loc) = false;
    temp = double(temp);
    exc_data.(fields{i})(:,M,:) = temp;
end

% set order priority first for sleep for male courtship lists
fields = {'wing_ext', 'wing_ext_all', 'court_chase', 'chase_all',...
              'circling_all','circling_1sec', 'CI','CI_all'};
for i = 1:length(fields)
        temp =  (exc_data.(fields{i}));
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
        loc = squeeze(exc_data.sleep(:,sex,:));
        temp(loc) = false;
        exc_data.(fields{i})(:,sex,:) = temp;
    end
end
% -----------------------------------------------------------------------------------

% when does each behavior happen within each region? 

d_fields = { 'FlyOnFood', 'OutterRing', 'foodQuad', 'sleep'}; % only a field for the male flies
s_fields = {'wing_ext', 'wing_ext_all', 'court_chase', 'chase_all', 'circling_all',...
    'circling_1sec', 'CI','CI_all'}; % only a field for the female flies

% Pull out the event onset times for each of the behaviors for each of the
% temperature regimes
behavior_onset = struct; % initialize as empty structure

% double fields
for i = 1:length(d_fields)
    for sex = 1:2 % for each of the flies
        field_name = d_fields{i};
        raw = squeeze(exc_data.(field_name)(:,sex,:));
        behavior_onset(sex).(field_name) = nan([num.trials,nTrans]);
        for roi = 1:nTrans % loop all temperature regions
            locs = raw(transitions(roi):transitions(roi+1),:); % all locations within the ROI that the behavior is seen
            for fly = 1:num.trials % find the first event time in this regime for each fly (if there are any)
                idx = find(locs(:,fly));
                if ~isempty(idx)
                    behavior_onset(sex).(field_name)(fly,roi) = idx(1);
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
        for fly = 1:num.trials % find the first event time in this regime for each fly (if there are any)
            idx = find(locs(:,fly));
            if ~isempty(idx)
                behavior_onset(sex).(field_name)(fly,roi) = idx(1);
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


% MALE FLY BEHAVIOR SEQUENCE 
fields = {'FlyOnFood', 'OutterRing', 'sleep', 'CI'};
kolors = {'gold', 'red', 'dodgerblue', 'green'}; % colors for the diff behaviors

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
    formatFig(fig, blkbgd,[r,c]);
    matchAxis(fig, true);
    for i = 1:nTrans
        subplot(r,c,i)
        set(gca, 'ycolor', 'none')
    end
    
    save_figure(fig, [figDir sexList{sex} ' behavior_onset per region'],fig_type);
end


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

save_figure(fig, [figDir 'behavior_onset per sex scatter'],fig_type);


%% Figure : combine across sex and multiple reps of the same temp regime
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
            'XJitterWidth',0.2,'MarkerFaceAlpha',0.5)
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


%% TODO: Plot out the sequence of behaviors color coded and then sort them
% by the the most of one behavior, then the next ext. so that it shows as a
% continuum 


%% TODO: sex based differences in the onset of different behavior times? 


%% TODO: where is the majority of courtship taking place? 
% Plot the locations for all instances of courtship to see if there is a
% spatial bias


%% TODO: How often do flies go from each behavior to the next? Make this a chord chart? 
% pull out the chunks of each behavior (like start and stop points for each
% 'time' of a behavior


