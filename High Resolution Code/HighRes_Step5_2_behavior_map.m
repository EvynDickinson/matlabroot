% todo
% - update this to include a map for what the female is doing -- have to
% think about how the courtship as a behavior state would be for her...
% most likely that would be best as simply a Male state and having her
% behavior being the other three states...

%% Behavior probability map
clearvars('-except',initial_var{:})

% 1 = on food
% 2 = sleeping
% 3 = edge-occupancy
% 4 = courtship
b_list = {'food', 'sleeping', 'edge','courtship'};
nstates = length(b_list);
ntime = length(time);

% male behavior states
behavior = nan([length(time),4]);
behavior(T.FlyOnFood(:,M),1) = 1;
behavior(m.sleep,2) = 2;
behavior(data(M).OutterRing,3) = 3;
behavior(T.CI,4) = 4;

% quick figure on the total set of observed behaviors
fig = getfig('',1,[655 406]);
bar(1:nstates,sum(~isnan(behavior)))
set(gca, "XTickLabel", b_list)
xlabel('behavior state')
ylabel('full trial count')
formatFig(fig);
save_figure(fig,[figDir 'Behavior state bar graph full ramp M'],fig_type);

% time course of the behavior states
fig = figure; hold on
for i = 1:4
    scatter(time, behavior(:,i),'filled')
end
xlabel('time (min)')
ylabel('behavior state')
ylim([0.5,nstates+0.5])
set(gca, "YTick",1:nstates, 'YTickLabel', b_list,'ydir', 'reverse')
formatFig(fig);
save_figure(fig,[figDir 'Behavior state timecourse M'],fig_type);


%% Build probabilty map
% % ======== build a dummy data set for testing ==========
% behavior = nan([ntime,4]);
% p = randperm(ntime);
% loc = 1:floor(ntime*0.2):ntime;
% for i = 1:4
%     rloc = loc(i):loc(i+1);
%     behavior(p(rloc),i) = i;
% end
% % preload the start of the behavior vector with specific test scenarios:
% a = [NaN 1 1 1 NaN 1 NaN NaN NaN 3 3 2 NaN 4 4 4 1 nan 1 1 nan 1 1 1 nan nan 1 1 1 1 1]; % easy start test
% beh(1:length(a))  = a;

% if sleep co-occurs with other behavior markers, sleep takes priority 
sleep_override = sum(~isnan(behavior(:,1:3)),2)>1 & ~isnan(behavior(:,2));
behavior(sleep_override,[1,3,4]) = nan; % override the other states since they can't be coincident with sleep

% Condense all the states into a single vector with numbers for each state
% at every point in time. nan represent frames with undetermined states
% TODO: refine this to check for more potential over laps etc.
beh = nan(size(time));
for i = [4,3,1,2] % ordered so that sleep is the last to fill and overrides the other states 
    loc = find(~isnan(behavior(:,i)));
    beh(loc) = i;
end

state_starts = find(~isnan(beh)); % each time point (frame #) that is an assigned state 

time_buff = 10; % number of frames that can be 'skipped' before the state is considered switched

% create a table that holds information on each state transition: 
tic 
ST = table('Size', [length(state_starts),7],...
        'VariableNames',{'state1', 'state1_start', 'state1_end', 'state2', 'state2_start','state2_end', 'transition'},...
        'VariableTypes', repmat({'double'},[1,7]));
idx = 1; % initialize the state transition #
for i = 1:length(state_starts)-1 % for each time point with a state (though the ultimate list will be shorter since repeat states count as 1)
    curr_state = beh(state_starts(i)); % current behavior state
    
    % skip this (i) if the current state start is in the preceding on-going state
    if idx-1>0
        if ST.state1_end(idx-1)>=state_starts(i) 
            continue
        end
    end

    % find the next state 
    for i_dt = 1:length(state_starts(i)+1:state_starts(end))
        t = state_starts(i+i_dt); % t = index for state_type moving forward from current one
        next_state = beh(t);
        frames_since_last_known_state = t-state_starts(i+i_dt-1); % how many nans from last state point
        
        % STATE CONTINUATION:
        % if the next state is the same as the current state and the last
        % time gap is within the acceptable range, skip to the next (i)
        % state as this one (current t) is a continuation of the same state
        if curr_state==next_state && (frames_since_last_known_state<=time_buff)
            continue % same state, look at next statepoints
        
        % NEW STATE
        % Register this as a state transition if either:
        % -- the next state is the same as the current one and the time gap is outside the allowable range 
        % -- or the state is different from the current one
        elseif (next_state==curr_state && (frames_since_last_known_state>time_buff)) || ~(next_state==curr_state)
            ST.state1(idx) = beh(state_starts(i)); % current state of the fly
            ST.state1_start(idx) = state_starts(i); % frame at which the current state started
            ST.state1_end(idx) = state_starts(i+i_dt-1);
            
            ST.transition(idx) = str2double([num2str(curr_state) num2str(next_state)]); % state transition number ID
            ST.state2(idx) = next_state; % next state of the fly
            ST.state2_start(idx) = t; % frame at which the next state begins
            if idx>1
                ST.state2_end(idx-1) = ST.state1_end(idx);
            end
            if idx==2
                ST.state2_end(1) = ST.state1_end(2);
            end
            idx = idx+1; 
            break % break i_dt loop since next state has been found
        end
    end
end
toc

% remove extra locations from the table
z = (ST.state1==0 & ST.state2==0);
ST(z,:) = [];

ntrans = size(ST,1); % how many transition exist in the experiment 

%% Histogram of the different state transitions & total count of the behaviors (~time)

bin_edges = min(ST.transition)-1:max(ST.transition)+1;
fig = figure; 
histogram(ST.transition,bin_edges)
xlabel('transition')
ylabel('count (#)')
fig = formatFig(fig);


%% Make the probability map as a heatmap with the transition probability as a shaded square with a number
% then it will also be easy to make a 'difference' map between each of the
% temperature periods (yes yes) 

% need to make a master transition list of all the possible ones so that a
% missing transitions can just be a zero not change the shape entirely.

transition_list = nan([nstates^2,1]);
idx = 1;
for i = 1:nstates
    for j = 1:nstates
        transition_list(idx) = str2double([num2str(i) num2str(j)]); 
        idx = idx + 1;
    end
end

transitions = [];
for i = 1:length(transition_list)
    transitions(i) = sum(ST.transition==transition_list(i));
end
transitions = reshape(transitions,[nstates, nstates])';

% quick transition matrix for the states
fig = getfig('',1,[492 331]);
heatmap(transitions)
xlabel('State 2')
ylabel('State 1')
set(gca, 'XDisplayLabels',b_list,'YDisplayLabels',b_list)
ax = gca;
ax.XDisplayLabels

save_figure(fig,[figDir 'Behavior state transitions heatmap M'],fig_type);


%% What are the transition probabilities within each temperature condition?

% start hold



%% 






































