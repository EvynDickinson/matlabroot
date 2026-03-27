

%% LOAD data
clear; clc;
warning off

auto_run = true;

% updates for updating the fly data with no swaps:
[excelfile, Excel, xlFile] = load_HighResExperiments;

% find trials that need to be updated:
done_loc = strcmpi('Y', excelfile(:,Excel.ToUpdate));
% done_loc = (strcmpi('Y', excelfile(:,Excel.swapcorrected)) | strcmpi('NA', excelfile(:,Excel.swapcorrected)));
loc_ready = strcmpi('Y', excelfile(:,Excel.groupready)) & ~done_loc;
excel_loc_list = find(loc_ready);

fprintf('\n %i unprocessed experiments remaining\n', length(excel_loc_list))

% select trial to update: 
trial_options = excelfile(excel_loc_list,Excel.trialID);
path = getDataPath(6,2);
baseFolder = [path,'Trial Data/'];

tic
for trial_idx = 1:length(trial_options)
    toc % time update for each trial
    fprintf('\n\n========== Processing trial %i of %i: %s ==========\n', ...
        trial_idx, length(trial_options), trial_options{trial_idx});

    excel_loc = excel_loc_list(trial_idx); % location in excel file of trials to process
    
    trialDir = trial_options{trial_idx};
    baseDir = [baseFolder, trialDir '/']; % full folder directory for that trial
    figDir = createFolder([baseDir,'Figures/']);
    
    processed_path = [baseDir 'post-5.1.2 data.mat'];
    % 
    % if isfile(processed_path)
    %     % Update Excel
    %     isExcelFileOpen(xlFile);
    %     writecell({'NA'},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.swapcorrected) num2str(excel_loc)]);
    %     fprintf('Trial already processed, skipping...\n');
    %     continue
    % end
    
    % Load basic data
    try
        load(processed_path)
    catch
        fprintf('ERROR: Could not load basic data for trial %s, skipping...\n', trialDir);
        continue
    end

    % Initialize variables
    auto_run = true;
    nvids = parameters.nVids;

%%
initial_var{end+1} = 'max_gap';
max_gap = 1;

%% ANALYSIS: Calculate male wing extension
clearvars('-except',initial_var{:})

% (accounts for suble drops in things like tracking)
Lwing = [];
Rwing = [];
% Determine which positions require which wing to be extended 
L_items = {'L1', 'L4', 'GX1', 'GX4', 'GY2', 'GY3'};
R_items = {'L2', 'L3', 'GX2', 'GX3', 'GY1', 'GY4'};
% Identify if male is in an appropriate position for each wing direction across each item
for i = 1:length(L_items)
    Lwing = [Lwing, position.(L_items{i})];
    Rwing = [Rwing, position.(R_items{i})];
end
% Condense to identify if male is in any of the appropriate positions
Lwing = any(Lwing,2);
Rwing = any(Rwing,2);

% Pull wing angles equal or greater than extension minimum for L and R
wa_cutoff = 50; % minimum wing extension angle for courtship
wing_ext = (Lwing & (m.wing.angle(:,1) >= wa_cutoff)) | (Rwing & (m.wing.angle(:,2) >= wa_cutoff)); % wing must be facing the female fly
wing_ext_filled = imclose(wing_ext, ones(max_gap + 1, 1)); % fill gaps in the behavior with those less than the max gap
% Each value subtracted by the value before it (1 = ext starts, -1 = ext stops, 0 = no state change)
a = diff(wing_ext_filled); 
% Add the first extension value to the list to account for the starting condition
b = [wing_ext_filled(1); a]; 
% Locations in wing_ext where extension period starts/end
ext_start = find(b == 1); 
ext_stop = find(b == -1);
% If wing ext doesn't stop by end, add stop location at end of ext_stop (loc = length of experiment value)
if wing_ext_filled(end)
    ext_stop(end + 1) = length(time);
end
% Calculate the length of each wing ext bout
ext_dur = ext_stop - ext_start;
% Find where wing ext lasts longer than 1sec
dur_loc = find(ext_dur > fps);

% Create new courtship matrix with only true wing ext for bouts longer than 1sec
mt = false(size(time));
for i = 1:length(dur_loc)
    ii = dur_loc(i);
    mt(ext_start(ii):ext_stop(ii)) = true;
end
T.wing_ext = mt;
T.wing_ext_filled = wing_ext_filled;
T.wing_ext_all = wing_ext;

%% ANALYSIS: Chase identification
% < 120 deg area behind female x
% facing female x
% 7mm between M center and F center x
% female speed > 0 x
% if all are 1 then courtship
% diff to see if chasing lasts > 2sec x
clearvars('-except',initial_var{:})

x = 1;
y = 2;

% positions of M head and F head and center
P1 = [m.pos(:,body.head,x),m.pos(:,body.head,y)]; % male head
P2 = [f.pos(:,body.center,x),f.pos(:,body.center,y)]; % female center
P3 = [f.pos(:,body.head,x),f.pos(:,body.head,y)]; % female head

% 1) Calculate body vectors
v1 = P3 - P1;  % Nx2 matrix, vector for female head to male head
v2 = P3 - P2;  % Nx2 matrix, vector for female head to female center

% 2) Calculate the dot product of v1 and v2 for each time step
dotProduct = v1(:,1) .* v2(:,1) + v1(:,2) .* v2(:,2);

% 3) Compute the magnitudes of the vectors
mag_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2); 
mag_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2); 

% 4) Calculate the cosine of the angle
cosTheta = dotProduct ./ (mag_v1 .* mag_v2);

% 5) Compute the angle in radians and convert to degrees
angleRadians = acos(cosTheta);  % angle in radians
angleDegrees = rad2deg(angleRadians);  % convert to degrees
mfpos_angle = angleDegrees;

% Identify when male position angle is less than 60 degrees from female
pos_angle = abs(mfpos_angle) <= 60;

% Identify when male is facing female
facing = [];
m_items = {'L3', 'L4', 'GX3', 'GX4', 'GY3', 'GY4'};
% Identify if male is in an appropriate position for each wing direction across each item
for i = 1:length(m_items)
    facing = [facing, position.(m_items{i})];
end
facing = any(facing,2);

% Identify when male is behind female AND facing her
mbehindf = (facing & pos_angle);

% Identify when male is within 7mm of female
close_dist = T.IFD <= 7; % mm

% Identify when female is moving
fmoving = f.speed >= 0.1; % min speed up for debate

% Identify when male is moving
mmoving = m.speed >= 0.1; % min speed diff than F in order to include true chase bouts

chase = (mbehindf & close_dist & fmoving & mmoving); % requirements for 'chase all'
chase_filled = imclose(chase, ones(max_gap + 1, 1)); % fill gaps in the behavior
% chase_filled = imclose(chase, ones(2, 1)); % fill gaps in the behavior

a = diff(chase_filled); 
% Add the first chase value to the list to account for the starting condition
b = [chase_filled(1); a]; 
% Locations in chase where chasing period starts/end
ch_start = find(b == 1); 
ch_stop = find(b == -1);
% If chasing doesn't stop by end, add stop location at end of ch_stop (loc = length of experiment value)
if chase_filled(end)
    ch_stop(end + 1) = length(time);
end
% Calculate the length of each chasing bout
ch_dur = ch_stop - ch_start;
% Find where chasing lasts longer than 2sec
dur_loc = find(ch_dur > (2*fps));

% Create new courtship matrix with only true chasing bouts longer than 2sec
mt = false(size(time));
if isempty(dur_loc)
    m.chaseroi = [];
else
    for i = 1:length(dur_loc)
        ii = dur_loc(i);
        mt(ch_start(ii):ch_stop(ii)) = true;
        m.chaseroi(i,:) = [ch_start(ii), ch_stop(ii)];
    end
end
T.court_chase = mt; % time restriction 2 seconds
T.chase_all = chase; % NO time limit
T.chase_filled = chase_filled;


%% ANALYSIS: Circling behavior
% M head within 3mm? [head_dist]
% M facing female [position.likely]
% M velocity constant within 1sec [const_var]
% F not moving
% if all are 1 then courtship
% diff to see if circling lasts > 1sec

clearvars('-except',initial_var{:})

% Distance between male head and female
x1 = m.pos(:,body.head,1); % x location for male head
y1 = m.pos(:,body.head,2);
x2 = f.pos(:,body.center,1); % x location for female center
y2 = f.pos(:,body.center,2);
% Calculate male head - female center distance
d = ((sqrt((x1-x2).^2 + (y1-y2).^2))./pix2mm);
head_dist = d <= 3; % mm
y = nan(size(head_dist));
y(head_dist) = 1;

% (for fun) how much time did the male fly spend really close to the female fly?
percent_time_within_3mm = (sum(head_dist)/length(time))*100;
disp(['The male fly spent ' num2str(percent_time_within_3mm) '% of the experiment close to the female'])

% position.all_likely

% female fly is not moving
f_speed_cut = (f.speed<=0.2);

% Determine if semiconstant male velocity
% 0.5sec smooth, then 2 sec std
binwidth = fps; %  1 second bins to look for avg speed
n = length(m.speed) - binwidth; % number of iterations
offset = ceil(binwidth * 0.5);
smoothspeed = smooth(m.speed,floor(0.5*fps),'moving');
smoothspeed = smooth(smoothspeed,floor(0.5*fps),'moving');
smoothspeed = smooth(smoothspeed,floor(0.3*fps),'moving');
smoothspeed = smooth(smoothspeed,floor(0.2*fps),'moving');
sp_var = nan(size(time)); % initialize the speed variability variable
for i = 1:n % loop through each bin
    roi = i:i+binwidth; % frames from i to one second later
    sp_var(i+offset) = std(smoothspeed(roi));
end

testlim = 0.7;
const_var = sp_var>testlim;
ok_var = m.speed;
ok_var(const_var) = nan;% all data points with acceptable variance

% -- View the threshold for consistent speed --
% r = 2;
% c = 1;
% 
% fig = figure('Position',[1 10 1478 845]);
% subplot(r,c,1); hold on
%     plot(time,m.speed,'color', Color('red'))
%     plot(time, ok_var, 'color', foreColor)
% subplot(r,c,2)
%     plot(time,sp_var,'color', foreColor)
% % formatting
% xlimits = [1,2];
% subplot(r,c,1)
%     xlim(xlimits)
%     ylabel('male fly speed (mm/s)')
% subplot(r,c,2)
%     xlim(xlimits)
%     h_line(testlim, 'r','-',1)
%     ylabel('speed variance')
% formatFig(fig,blkbnd,[r,c])

% Full selection criteria: 
V = position.likely & head_dist & const_var & f_speed_cut;
V_filled = imclose(V, ones(max_gap + 1, 1)); % fill micro gaps in the behavior

%  --- Find periods longer than 1 second ---
a = diff(V_filled); % when does the speed switch between stability and instability 
% Add the first chase value to the list to account for the starting condition
b = [V_filled(1); a]; 
% Find when stability periods starts/end
v_start = find(b == 1); 
v_stop = find(b == -1);
% If speed stability doesn't stop by end, add stop location at end of v_stop (loc = length of experiment value)
if V_filled(end)
    v_stop(end + 1) = length(time);
end
% Calculate the length of each speed bout
v_dur = v_stop - v_start;
% Find where speed stability lasts longer than 1sec
dur_loc = find(v_dur > fps);
% save the constant speed regions into a new metric
constant_velocity = false(size(time));
for i = 1:length(dur_loc)
    idx = dur_loc(i);
    constant_velocity(v_start(idx): v_stop(idx)) = true;
end

T.circling_all = V; % when the male fly is circling the female (no time restriction)
T.circling_1sec = constant_velocity; % when circling is longer than 1 second
T.circling_filled = V_filled;

% COURTSHIP INDEX
CI = any([T.court_chase,T.wing_ext,T.circling_1sec],2); % courtship index
T.CI = CI;


%% ANALYSIS: Behavior probability map
clearvars('-except',initial_var{:})

% 1 = on food
% 2 = sleeping
% 3 = edge-occupancy
% 4 = courtship

% Behavior states
b_list = {'food', 'sleeping', 'edge','m_courtship'};
nStates = length(b_list);
ntime = length(time);
time_buff = 10; % number of frames that can be 'skipped' before a state switch

for sex = 1:2 % male and female

    % ---- Initialize behavior state matrix -----
    behavior = nan(ntime, nStates);
    % Food: 
    behavior(T.FlyOnFood(:, sex), 1) = 1;
    % Sleep:
    if sex == 1
        behavior(m.sleep, 2) = 2;
    else
        behavior(f.sleep, 2) = 2;
    end
    % Escape: 
    behavior(data(sex).OutterRing, 3) = 3;
    % Courtship: 
    behavior(T.CI, 4) = 4;

    % Sleep takes priority when overlapping with other states
    sleep_override = ~isnan(behavior(:, 2));
    behavior(sleep_override, [1, 3, 4]) = nan;

    % Assign a single state per time point
    beh = nan(ntime,1);
    for ii = [4, 3, 1, 2] % Sleep is assigned last to take priority
        beh(~isnan(behavior(:, ii))) = ii;
    end

    % ------ Find state transitions -----

    % Remove NaNs and get state changes
    valid_idx = find(~isnan(beh));
    valid_states = beh(valid_idx);
    
    % Find where state changes (accounting for time_buff)
    state_changes = [true; (diff(valid_states) ~= 0) | (diff(valid_idx) > time_buff)];
    
    % Get transition indices
    trans_start = valid_idx(state_changes);
    trans_end = [trans_start(2:end)-1; ntime];
    states = valid_states(state_changes);
    
    % Build state transition table
    n_trans = length(trans_start) - 1;
    ST = table(...
        states(1:end-1), trans_start(1:end-1), trans_end(1:end-1), ...
        states(2:end), trans_start(2:end), trans_end(2:end), ...
        states(1:end-1)*10 + states(2:end), ...
        'VariableNames', {'state1', 'state1_start', 'state1_end', ...
                          'state2', 'state2_start', 'state2_end', 'transition'});

    % ===== COMPUTE TRANSITION MATRIX  =====
    % Use accumarray instead of loops
    transition_rows = floor(ST.transition / 10);
    transition_cols = mod(ST.transition, 10);
    transitions = accumarray([transition_rows, transition_cols], 1, [nStates, nStates]);
    
    % Create transition list for reference
    [i_grid, j_grid] = ndgrid(1:nStates, 1:nStates);
    transition_list = i_grid(:)*10 + j_grid(:);

    % ===== SAVE DATA =====
    data(sex).states.ST = ST;
    data(sex).states.b_list = b_list;
    data(sex).states.nstates = nStates;
    data(sex).states.beh = beh;
    data(sex).states.behavior = behavior;
    data(sex).states.transitions = transitions;
    data(sex).states.transition_list = transition_list;

end

% How long to each behavior state after each temp change and what order?
for sex = 1:2  % for each of the two flies
    nStates = data(sex).states.nstates;
    nRegimes = numel(tRate); % to account for each temperature regime
    ST = data(sex).states.ST;
    beh = data(sex).states.beh;

    freq = struct('name', cell(nRegimes, nStates), ...
                  'n_states', {0}, ...
                  'n_instances', {0}, ...
                  'perc_time_in_state', {0}, ...
                  'time_to_first_instance', {nan}, ...
                  'temp_at_first_instance', {nan});

    for r = 1:nRegimes  % for each temperature regime
        roi = tRate(r).idx;  % get the temperature regime frame indices
        total_frames = diff(roi);
        roi_range = roi(1):roi(2);
        for ii = 1:nStates  % for each type of behavior/state
            freq(r, ii).name = b_list{ii};
            
            % Find states in this regime
            state_in_regime = (ST.state1 == ii) & ...
                              (ST.state1_start >= roi(1)) & ...
                              (ST.state1_start <= roi(2));

            if ~any(state_in_regime)
                continue;  % Skip if behavior never occurred
            end
            
            % Find first occurrence
            first_idx = find(state_in_regime, 1, 'first');
            first_frame = ST.state1_start(first_idx);

            % Number of occurrences
            freq(r, ii).n_states = sum(state_in_regime);
            freq(r, ii).n_instances = sum(beh(roi_range) == ii);
            freq(r, ii).perc_time_in_state = ...
                            (freq(r, ii).n_instances / total_frames) * 100;
            freq(r, ii).time_to_first_instance = ...
                            (first_frame - roi(1)) / parameters.FPS;
            
            % Safe indexing for temperature
            temp_idx = first_frame:min(first_frame+3, ntime);
            freq(r, ii).temp_at_first_instance = mean(T.temperature(temp_idx));

        end
    end
    
    % Store results in data structure
    data(sex).states.freq = freq;
end

%%  Save the 'analyzed' data package:
    clearvars('-except',initial_var{:})
    
    try
        save([baseDir 'post-5.1.3 data.mat'],'-v7.3')
        % Update Excel
        isExcelFileOpen(xlFile);
        writecell({swap_status},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.ToUpdate) num2str(excel_loc)]);
        fprintf('✓ Saved %s data file\n', trialDir)
    catch ME
        fprintf('ERROR saving trial %s: %s\n', trialDir, ME.message);
    end
    
    % Clean up for next iteration
    clearvars('-except', 'trial_idx', 'trial_options', 'excel_loc_list', 'excel_loc', ...
              'xlFile', 'Excel', 'excelfile', 'path', 'baseFolder')

end

