% todo
% - update this to include a map for what the female is doing -- have to
% think about how the courtship as a behavior state would be for her...
% most likely that would be best as simply a Male state and having her
% behavior being the other three states...

% % ======== dummy data set for testing ==========
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


%% Behavior probability map
clearvars('-except',initial_var{:})

% 1 = on food
% 2 = sleeping
% 3 = edge-occupancy
% 4 = courtship
b_list = {'food', 'sleeping', 'edge','m_courtship'};
nstates = length(b_list);
ntime = length(time);
time_buff = 10; % number of frames that can be 'skipped' before the state is considered switched

for sex = 1:2 % male and female
   
    % Fill in behavior states %TODO HERE: working to align everything for
    % the loop
    behavior = nan([length(time),4]);
    behavior(T.FlyOnFood(:,sex),1) = 1;
    switch sex
        case 1
            behavior(m.sleep,2) = 2;
        case 2
            behavior(f.sleep,2) = 2;
    end
    behavior(data(sex).OutterRing,3) = 3;
    behavior(T.CI,4) = 4;
    
    % Build probabilty map
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
    % remove extra locations from the table
    z = (ST.state1==0 & ST.state2==0);
    ST(z,:) = [];
    toc

    % Make the probability map as a heatmap with the transition probability as a shaded square with a number
    % then it will also be easy to make a 'difference' map between each of the
    % temperature periods (yes yes) 
    
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
    
    % Save some of the data to the data structure for future use:
    data(sex).states.ST = ST;
    data(sex).states.b_list = b_list;
    data(sex).states.nstates = nstates;
    data(sex).states.beh = beh;
    data(sex).states.behavior = behavior;
    data(sex).states.transitions = transitions;
    data(sex).states.transition_list = transition_list;

end

% 

%% FIGURES: number and distribution of behaviors that were observed
% quick figure on the total set of observed behaviors
clearvars('-except',initial_var{:})

sexes = {'male', 'female'};
nstates = data(sex).states.nstates;
b_list = strrep(data(1).states.b_list,'_',' ');
fig = getfig('',1,[ 937 406]);
for sex = 1:2
    subplot(1,2,sex)
    bar(1:nstates, sum(~isnan(data(sex).states.behavior)),'FaceColor',data(sex).color)
    set(gca, "XTickLabel", b_list)
    xlabel('behavior state')
    ylabel('full trial count')
end
formatFig(fig,false,[1,2]);
save_figure(fig,[figDir 'Behavior state bar graph full ramp  M and F'],fig_type);

% time course of the behavior states
fig = getfig('',1,[ 937 406]);
for sex = 1:2
    subplot(1,2,sex);
    hold on
    h = rectangle('Position', [roi(1), ylims(1), diff(roi),diff(ylims)], 'FaceColor', tRate(i).color);

    for i = 1:4
        y = data(sex).states.behavior(:,i);
        scatter(time, y, 35, data(sex).color, 'filled')
    end
    xlabel('time (min)')
    ylabel('behavior state')
    ylim([0.5,nstates+0.5])
    set(gca, "YTick",1:nstates, 'YTickLabel',  b_list, 'ydir', 'reverse')
end
formatFig(fig,false, [1,2]);
save_figure(fig,[figDir 'Behavior state timecourse M and F'],fig_type);

% FIGURE: transition matrix heatmap of the states
fig = getfig('',1);
for sex = 1:2
    subplot(1,2,sex);
    heatmap(data(sex).states.transitions)
    xlabel('State 2')
    ylabel('State 1')
    set(gca, 'XDisplayLabels',b_list,'YDisplayLabels',b_list)
    title(sexes{sex})
end
save_figure(fig,[figDir 'Behavior state transitions heatmap M and F'],fig_type);

% Figures: Histogram of the different state transitions & total count of the behaviors (~time)
fig = getfig('',1,[1064 473]);
for sex = 1:2
    subplot(1,2,sex)
    bin_edges = min(data(sex).states.ST.transition)-1:max(data(sex).states.ST.transition)+1;
    histogram(data(sex).states.ST.transition,bin_edges,'FaceColor',data(sex).color)
    xlabel('transition')
    ylabel('count (#)')
end
fig = formatFig(fig,false,[1,2]);
save_figure(fig,[figDir 'Behavior state transitions histogram M and F'],fig_type);

% Figure:  pie chart of different behaviors across the full experiment
fig = getfig('',1);
for sex = 1:2
    subplot(1,2,sex)
    y = data(sex).states.beh;
    x = sum(isnan(y));
    for i = 1:data(sex).states.nstates
        x = [x, sum(y==i)];
    end
    xplode = true(size(x));
    xplode(1) = false; % do not explode out the unknown states/positions
    pie(x,xplode)
    title(sexes{sex})
end


%% Save data 
clearvars('-except',initial_var{:})

switch questdlg('Save processed data?')
    case 'Yes'
        save([baseDir 'post-5.1 data.mat'],'-v7.3')
        disp('Saved data file')
    case 'No'
        return
    case 'Cancel'
        return
end








%% How long to each behavior state after each temp change and what order?
freq = [];
for sex = 1:2
    for r = 1:size(tRate,2)
        roi = tRate(r).idx; % get the temperature regime frame indeces 
        for i = 1:nstates
            loc = find(data(sex).states.ST.state1==i);
            % -- make a list of all the times for this behavior (so we can do
            % a frequency graph etc) 

            % -- pull time of the first one

            % -- pull the temp for the first one


        end
    end
end
% when (if) does each state arise in each temp regime? At what temp is the
% behavior observed?


%% Pointed Question: when after the onset of heat does the 

% how long after the onset of each temperature event does it take for flies
% to return to food?

for sex = 1:2
    
    data(sex).states.





































