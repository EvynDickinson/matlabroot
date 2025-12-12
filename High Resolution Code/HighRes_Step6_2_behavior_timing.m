

% create plots about the order of behaviors / priority of behaviors

%% Time delay to each behavior

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
figure; 
plot(data.temp);
v_line(transitions)

% when does each behavior happen within each region? 

% TODO: make a demo logical structure that includes the behavior heirarchy
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
   
% set order priority
exc_data. % TODO here


behavior_onset = struct; % initialize as empty structure

for sex = 1:2 % for each of the flies
    % sleep 
    raw = squeeze(data.sleep(:,sex,:));
    behavior_onset(sex).sleep = nan([num.trials,length(transitions)]);
    for roi = 1:length(transitions)-1 % loop all temperature regions
        locs = (raw(transitions(roi):transitions(roi+1),:)); % all locations within the ROI that the behavior is seen
        for fly = 1:num.trials
            idx = find(locs(:,fly));
            if ~isempty(idx)
                behavior_onset(sex).sleep(fly,roi) = idx(1);
            end
        end
    end
    % sleep 
    raw = squeeze(data.sleep(:,sex,:));
    behavior_onset(sex).sleep = nan([num.trials,length(transitions)]);
    for roi = 1:length(transitions)-1 % loop all temperature regions
        locs = (raw(transitions(roi):transitions(roi+1),:)); % all locations within the ROI that the behavior is seen
        for fly = 1:num.trials
            idx = find(locs(:,fly));
            if ~isempty(idx)
                behavior_onset(sex).sleep(fly,roi) = idx(1);
            end
        end
    end

end



% scatter plot of the onset timing within each of the regions
fig = getfig; 
hold on


% behaviors: 
% sleep, on food, near food, courtship (M only), escape (ring)



data.