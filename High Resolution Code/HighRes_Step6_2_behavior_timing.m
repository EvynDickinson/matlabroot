

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
        behavior_onset(sex).(field_name) = nan([num.trials,length(transitions)-1]);
        for roi = 1:length(transitions)-1 % loop all temperature regions
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
    behavior_onset(sex).(field_name) = nan([num.trials,length(transitions)-1]);
    for roi = 1:length(transitions)-1 % loop all temperature regions
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

% subplots for each regime

fields = { 'FlyOnFood', 'OutterRing', 'sleep', 'CI'};

fig = getfig; 
for i = 1 : length(transitions)-1
    subplot(r,c,i)
    hold on
    % working here: set this up as a time vs row=fly
end








% behaviors: 
% sleep, on food, near food, courtship (M only), escape (ring)



data.