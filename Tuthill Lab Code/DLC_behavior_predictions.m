


%% Behavior predictions from the FeTi joint angle of a single leg
% then prediction of full fly behavior from all the leg information

controlROI = 1:param.basler_delay*fps; % all frames before the laser

% ----- Input -----
iJ = 2;             % joint used for step freq calculations
min_freq = 5;       % minimum step frequency for 'walking' %hard lim 5 = 2 steps of the L1 leg
min_diff = 30;      % minimum joint angles oscillation amp for 'walking'
max_std = 15;       % maximum step frequency std for 'walking
max_variance = 0.1; % max variance in joint angle for 'stationary'
max_range = 5;      % max range in joint angle for 'stationary'
behavior_leg = 1;   % leg used to assign the behavior
% -----------------

% Determine the behavior state of each leg:
for ifly = 1:num.flies
 for cond = 1:num.conds
    for rep = 1:num.reps
        for leg = 1:6 %number of legs
            raw = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
            angles(ifly).Leg(leg).control_behavior{cond,rep} = []; %set blank
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
            angles(ifly).Leg(leg).step(cond,rep).max_loc = max_loc;
            angles(ifly).Leg(leg).step(cond,rep).max_val = max_val;
            angles(ifly).Leg(leg).step(cond,rep).min_loc = min_loc;
            angles(ifly).Leg(leg).step(cond,rep).min_val = min_val;
            
            % Find control step locations and frequency:
            frames = round(max_loc(max_loc<controlROI(end))); %loc of step
            freq = length(frames)/param.basler_delay; % step frequency
            if freq >= min_freq
                step_std = std(diff(frames));
            else
                step_std = nan;
            end
            % add information to the structure:
            angles(ifly).Leg(leg).control(cond,rep).frames = frames; 
            angles(ifly).Leg(leg).control(cond,rep).freq = freq;
            angles(ifly).Leg(leg).control(cond,rep).step_std = step_std;
            
            % Categorize single leg behavior:
            %walking
            if freq >= min_freq && step_std <= max_std
                angles(ifly).Leg(leg).control_behavior{cond,rep} = 'walking';
            end
            %stationary  
            cntl_angles = raw(controlROI);
            a = abs(mean(diff(cntl_angles))) < max_variance; 
            b = (max(cntl_angles)- min(cntl_angles)) < max_range;
            if a && b == true
                angles(ifly).Leg(leg).control_behavior{cond,rep} = 'stationary';
            end
        end
     end    
  end
end


% Determine the fly's behavior based on the leg state:
for ifly = 1:num.flies
    % set blank behavior group:
    [behavior(ifly).walking, behavior(ifly).stationary, behavior(ifly).other,...
            behavior(ifly).grooming] = deal(false(num.conds,num.reps)); 
    behavior(ifly).behavior = cell(num.conds,num.reps);
    behavior(ifly).class = nan(num.conds,num.reps);
    for cond = 1:num.conds
        for rep = 1:num.reps
            MT = 4*ones(6,1); % other = 4;
            for leg = 1:6
                data(leg) = angles(ifly).Leg(leg).control_behavior(cond,rep);
            end
            % pull the state for each leg:
            MT(strcmp(data, 'walking')) = 1; % walking = 1;
            MT(strcmp(data, 'stationary')) = 2; % stationary = 2;
            behavior(ifly).leg_behavior(cond).data(:,rep) = MT;
            
        % Categorize the trial behavior: 
            %walking:
            if sum(MT==1)>3 && sum(MT==2)==0 % 3+ legs walking & 0 stationary
                behavior(ifly).walking(cond,rep) = true;
                behavior(ifly).behavior{cond,rep} = 'walking';
                behavior(ifly).class(cond,rep) = 1;
            %stationary:
            elseif sum(MT==2)==6 % all legs stationary
                behavior(ifly).stationary(cond,rep) = true;
                behavior(ifly).behavior{cond,rep} = 'stationary';
                behavior(ifly).class(cond,rep) = 2;
            %grooming
            elseif sum(MT==2)==3 || sum(MT==2)==4 % 3 or 4 stationary legs
                behavior(ifly).grooming(cond,rep) = true;
                behavior(ifly).behavior{cond,rep} = 'grooming';
                behavior(ifly).class(cond,rep) = 3;
            %other
            else
                behavior(ifly).other(cond,rep) = true;
                behavior(ifly).behavior{cond,rep} = 'other';
                behavior(ifly).class(cond,rep) = 4;
            end
        end
    end
end


clearvars('-except',initial_vars{:})


%% Figure with leg behavior labels for each condition and rep for a given fly
% ---- input -----
ifly = 9;
% ----------------
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
fig = getfig('',1);
set(fig, 'color', 'k')
for cond = 1:num.conds
    subplot(4,7,cond)
    set(gca, 'color', 'k')
    xlim([0,num.reps+1])
    ylim([0,7]) %number of legs
    hold on

    for rep = 1:num.reps
        % plot the trial behavior
        switch behavior(ifly).class(cond,rep)
            case 1; C = kolor.walk; 
            case 2; C = kolor.stat; 
            case 3; C = kolor.groom;
            case 4; C = kolor.other;
        end
        plot([rep,rep], [0.5,6.5], 'color', C, 'linewidth', LW)
        
        % plot the leg behavior
        input = behavior(ifly).leg_behavior(cond).data; 
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
    
    title(['Cond ' num2str(cond)])
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
%     ylabel('Leg')
%     xlabel('Trial')
    set(gca,'XColor', 'none','YColor','none')
end

save_figure(fig, [fig_dir '\' FilePath.structure_name ' joint behavior label example'], '-png');

clearvars('-except',initial_vars{:})


%% compare determined behaviors with those manually labeled:
% Add behavior class fields to the structure:
load(['C:\matlabroot\behavior class\' FilePath.cross ' behavior class']);

for ifly = 1:num.flies
    input = group(ifly).behavior;
    group(ifly).walking = strcmpi(input, 'Walking');
    group(ifly).stationary = strcmpi(input, 'Stationary');
    group(ifly).other = strcmpi(input, 'Other');
    group(ifly).grooming = strcmpi(input, 'Grooming');
end



% compare the manually labeled data and the auto-labeled data:
for ifly = 1:num.flies
   % find locations with matches between manual and auto:
   behavior(ifly).walk_data.match_loc = ...
       (behavior(ifly).walking & group(ifly).walking == true);
   % sum of matches
   behavior(ifly).walk_data.match_tot = sum(sum(behavior(ifly).walk_data.match_loc));
   % total manual 'walking' count
   behavior(ifly).walk_data.match_tot = sum(sum(group(ifly).walking));
   % total auto 'walking' count
   behavior(ifly).walk_data.match_tot = sum(sum(behavior(ifly).walking));
end


% Visualize the labeling information:
fig_title = [FilePath.structure_name ' auto-labeled vs manual-labeled behavior'];
fig = getfig(fig_title,1);
[nrows,ncols] = subplot_numbers(num.flies);
for ifly = 1:num.flies
    subplot(nrows,ncols,ifly)
    % convert the logical data into an array:
    MT = nan(num.conds, num.reps);
    MT(behavior(ifly).walking) = 1;                 % auto labeled = blue
    MT(group(ifly).walking) = 3;                    % manual labeled = yellow
    MT(behavior(ifly).walk_data.match_loc) = 2;     % overlap = green
    %plot the data
    imAlpha=MT;
    imAlpha(isnan(MT)) = 0;
    imagesc(MT,'AlphaData',imAlpha);
    set(gca,'color',1*[1 1 1]); %set background color to white
    clims = [1, 3];
    caxis(clims);
    % Adjust colorbar
%     c = colorbar;
%     c.Label.String = 'Labeling';
    xlabel('rep')
    ylabel('cond')  
    tag = FilePath.locations{ifly,4};
    title(strrep(tag, '_', '-'))
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

save_figure(fig, [fig_dir '\' fig_title]);



% Filter the automatic labeling data to eliminate auto-labeled 'walking' to
% 'other' if visual labeling doesn't show walking

for ifly = 1:num.flies
   man = group(ifly).walking; %manual labeled points
   aut = behavior(ifly).walking; %auto-labeled
   
   loc = (man==false) & (aut==true); % find locations
   behavior(ifly).walking(loc) = false; % false for walking
   behavior(ifly).other(loc) = true; % true for other
   behavior(ifly).class(loc) = 3; % change class to other
    
end

clearvars('-except',initial_vars{:})








































