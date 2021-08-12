
function [behavior, figHandle] = DLC_behaviorpredictions(angles,param,ROI)
% [behavior, figHandle] = DLC_behaviorpredictions(Angles,param)



%% Behavior predictions from the FeTi joint angle of a single leg
% then prediction of full fly behavior from all the leg information

% controlROI = 1:param.basler_delay*fps; % all frames before the laser
controlROI = ROI;
num = NUM(param);

% ----- Adjustable Input -----
iJ = 2;             % joint used for step freq calculations
min_freq = 5;       % minimum step frequency for 'walking' %hard lim 5 = 2 steps of the L1 leg
min_diff = 30;      % minimum joint angles oscillation amp for 'walking'
max_std = 15;       % maximum step frequency std for 'walking
max_variance = 0.1; % max variance in joint angle for 'stationary'
max_range = 5;      % max range in joint angle for 'stationary'
% -----------------

% Determine the behavior state of each leg:
for cond = 1:num.conds
   for rep = 1:num.reps
      for leg = 1:6 %number of legs
            raw = angles.Leg(leg).data(cond,rep).all(:,iJ);
            behavior.Leg(leg).control_behavior{cond,rep} = []; %set blank
            % Find the peaks on the leg angle trace:
            max_loc = findPeaks('-x', raw, '-extrema', 1);
            max_val = raw(round(max_loc));
            min_loc = findPeaks('-x', raw, '-extrema', -1);
            min_val = raw(round(min_loc));
            % check that the peaks qualify for the minimum oscillation amplitude:
            oscillation_diff = max(max_val) - min(min_val);
            if oscillation_diff <= min_diff
               [max_loc, max_val, min_loc, min_val] = deal([]);
            end
            % add the peak information to the structure:
            behavior.Leg(leg).step(cond,rep).max_loc = max_loc;
            behavior.Leg(leg).step(cond,rep).max_val = max_val;
            behavior.Leg(leg).step(cond,rep).min_loc = min_loc;
            behavior.Leg(leg).step(cond,rep).min_val = min_val;
            
            % Find control step locations and frequency:
            frames = round(max_loc(max_loc<controlROI(end))); %loc of step
            ROI_time = length(controlROI)/param.Basler_fps;
            freq = length(frames)/ROI_time; % step frequency
            if freq >= min_freq
                step_std = std(diff(frames));
            else
                step_std = nan;
            end
            % add information to the structure:
            behavior.Leg(leg).control(cond,rep).frames = frames; 
            behavior.Leg(leg).control(cond,rep).freq = freq;
            behavior.Leg(leg).control(cond,rep).step_std = step_std;
            
            % Categorize single leg behavior:
            %walking
            if freq >= min_freq && step_std <= max_std
                behavior.Leg(leg).control_behavior{cond,rep} = 'walking';
            end
            %stationary  
            cntl_angles = raw(controlROI);
            a = abs(mean(diff(cntl_angles))) < max_variance; 
            b = (max(cntl_angles)- min(cntl_angles)) < max_range;
            if a && b == true
                behavior.Leg(leg).control_behavior{cond,rep} = 'stationary';
            end
      end
   end    
end



% Determine the fly's behavior based on the leg state:
% set blank behavior group:
[behavior.walking, behavior.stationary, behavior.other,...
        behavior.grooming] = deal(false(num.conds,num.reps)); 
behavior.behavior = cell(num.conds,num.reps);
behavior.class = nan(num.conds,num.reps);
for cond = 1:num.conds
    for rep = 1:num.reps
        MT = 4*ones(6,1); % other = 4;
        for leg = 1:6
            data(leg) = behavior.Leg(leg).control_behavior(cond,rep);
        end
        % pull the state for each leg:
        MT(strcmp(data, 'walking')) = 1; % walking = 1;
        MT(strcmp(data, 'stationary')) = 2; % stationary = 2;
        behavior.leg_behavior(cond).data(:,rep) = MT;

    % Categorize the trial behavior: 
        %walking:
        if sum(MT==1)>3 && sum(MT==2)==0 % 3+ legs walking & 0 stationary
            behavior.walking(cond,rep) = true;
            behavior.behavior{cond,rep} = 'walking';
            behavior.class(cond,rep) = 1;
        %stationary:
        elseif sum(MT==2)==6 % all legs stationary
            behavior.stationary(cond,rep) = true;
            behavior.behavior{cond,rep} = 'stationary';
            behavior.class(cond,rep) = 2;
        %grooming
        elseif sum(MT==2)==3 || sum(MT==2)==4 % 3 or 4 stationary legs
            behavior.grooming(cond,rep) = true;
            behavior.behavior{cond,rep} = 'grooming';
            behavior.class(cond,rep) = 3;
        %other
        else
            behavior.other(cond,rep) = true;
            behavior.behavior{cond,rep} = 'other';
            behavior.class(cond,rep) = 4;
        end
    end
end


%% Figure with leg behavior labels for each condition and rep for a given fly
a = (6:-1:1)';
MT_Y = repmat(a,1,num.reps);
a = 1:num.reps;
MT_X = repmat(a,6,1);

% compare behavior labels across legs:
kolor.walk = Color('darkcyan');
kolor.stat = Color('firebrick');
kolor.groom = Color('purple');
kolor.other = Color('grey');
SZ = 50;
LW = 10;
figHandle = getfig('',1);
set(figHandle, 'color', 'k')
[nrows, ncols] = subplot_numbers(num.conds);
for cond = 1:num.conds
    subplot(nrows,ncols,cond)
    set(gca, 'color', 'k')
    xlim([0,num.reps+1])
    ylim([0,7]) %number of legs
    hold on
    for rep = 1:num.reps
        % plot the trial behavior
        switch behavior.class(cond,rep)
            case 1; C = kolor.walk; 
            case 2; C = kolor.stat; 
            case 3; C = kolor.groom;
            case 4; C = kolor.other;
        end
        plot([rep,rep], [0.5,6.5], 'color', C, 'linewidth', LW)
        
        % plot the leg behavior
        input = behavior.leg_behavior(cond).data; 
        % plot walking points:
        loc = input==1;
        scatter(MT_X(loc), MT_Y(loc), SZ, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', kolor.walk)
        
        % plot stationary points:
        loc = input==2;
        scatter(MT_X(loc), MT_Y(loc), SZ, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', kolor.stat)
        
        % plot other points:
        loc = input==4;
        scatter(MT_X(loc), MT_Y(loc), SZ, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', kolor.other)
    end
    
    condtitle = sprintf('\\color{%s}%s', 'cyan', num2str(cond));
    title_list = ['\color{white}Cond: ' condtitle];
    title(title_list)
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    ylabel('Leg')
    xlabel('Trial')
    set(gca,'XColor', 'w','YColor','w')
    
end



end













