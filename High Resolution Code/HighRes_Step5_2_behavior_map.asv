

%% Behavior probability map

% 1 = on food
% 2 = sleeping
% 3 = edge-occupancy
% 4 = courtship

behavior = nan([length(time),4]);

behavior(T.FlyOnFood(:,M),1) = 1;
behavior(m.sleep,2) = 2;
behavior(data(M).OutterRing,3) = 3;
behavior(T.CI,4) = 4;

figure;
for i = 1:4
    scatter(time, behavior(:,i))
end

sum(data(M).OutterRing)

% dummy data set: 
ntime = length(time);

behavior = nan([ntime,4]);
p = randperm(ntime);
loc = 1:floor(ntime*0.2):ntime;
for i = 1:4
    rloc = loc(i):loc(i+1);
    behavior(p(rloc),i) = i;
end

% 
figure; hold on
for i = 1:4
    scatter(time, behavior(:,i))
end


% what are the transitions? 

% 
% if sleep co-occurs with other behavior markers, sleep takes priority 
sleep_override = sum(~isnan(behavior(:,1:3)),2)>1 & ~isnan(behavior(:,2));
behavior(sleep_override,[1,3,4]) = nan; % override the other states since they can't be coincident with sleep




% TODO: set everything in these locations to just sleeeeeep
beh = nan(size(time));
for i = 4:-1:1
    loc = find(~isnan(behavior(:,i)));
    beh(loc) = i;
end
beh(1:17) = [NaN 1 1 1 NaN 1 NaN NaN NaN 3 3 2 NaN 4 4 4 1]; % easy start test

state_starts = find(~isnan(beh)); % each time point that is an assigned state 

% TODO
% 1) for each state point, find the next state
% 2) if that state is the same one, determine if it is within the same
% continuous period (need some preset allowable 'buffer') and if so, skip
% to the end and look for the next state
% 3) note the 'type' of transition -- need to create a number key for each
% type of transition (could be simply #1 and #2, e.g., from state 1 to
% state 2, that would be transition 12)...

time_buff = 2; % number of frames that can be 'skipped' before the state is considered switched

ST = table('Size', [length(state_starts),7],...
        'VariableNames',{'state1', 'state1_start', 'state1_end', 'state2', 'state2_start','state2_end', 'transition'},...
        'VariableTypes', repmat({'double'},[1,7]));
idx = 1;
for i = 1:length(state_starts)-1 % for each time point with a state
    curr_state = beh(state_starts(i)); % current behavior state
    % find the next state
    for i_dt = 1:length(state_starts(i)+1:state_starts(end))
        t = state_starts(i+i_dt); % t = index for state_type moving forward from current one
        next_state = beh(t);
        frames_since_last_known_state = t-state_starts(i+i_dt-1); % how many nans from last state point

        % if the next state is the same as the current state and the last
        % time gap is within the acceptable range, skip to the next (i) state
        if curr_state==next_state && (frames_since_last_known_state<=time_buff)
            continue % same state, look at next statepoints
        
        % Register this as a state transition if either:
        % -- the next state is the same as the current one and the time gap is outside the allowable range
        % -- the state is different from the current one
        elseif (next_state==curr_state && (t-state_starts(i)<=time_buff)) || ~(next_state==curr_state)
            ST.state1(idx) = beh(state_starts(i));
            ST.state1_start(idx) = state_starts(i);
            
            ST.transition(idx) = str2double([num2str(curr_state) num2str(next_state)]);
            ST.state2(idx) = next_state;
            ST.state2_start(idx) = t;
            idx = idx+1;
        end
        
    end
end

            
i_dt = i_dt+1;
   

    next_state = (state_starts(i)+1:ntime); % next behavior state












































