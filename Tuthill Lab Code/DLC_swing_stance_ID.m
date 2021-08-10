% Use this script to find the virtual sphere of the treadmill to easily
% identify the location of the 'ground' to better classify swing / stance.
% This uses the predictions of walking behavior from the
% DLC_behavior_predictions script to isolate frames that are 'on' the ball
% to use for fitting a sphere


clearvars('-except',initial_vars{:})

% ----- input ----- 
ifly = 6;
LP = 5; % leg point (1-5)  5 = tarsus 
controlROI = 1:param.basler_delay*fps;
% -----------------
var_list = [];
var_list = initial_vars;
% variables that need to be saved during the extraction for loops
new_vars = {'var_list';'cond'; 'ifly'; 'rep'; 'trial';...
            'controlROI'; 'LP'; 'idx'; 'ntotal'; 'walkLoc'};

var_list(end+1:end+length(new_vars)) = new_vars;
% clearvars('-except',var_list{:})

% build a folder for this specific fly: 
temp_folder = [fig_dir '\' FilePath.locations{ifly,4}];
if ~exist(temp_folder, 'dir')
    mkdir(temp_folder)
end

%% extract frames from WALKING trials:
ntotal = sum(sum(behavior(ifly).walking)); % total walking trials
% pull the list of walking: 
walkLoc = [];
for rep = 1:num.reps
   a(:,1) = find(behavior(ifly).walking(:,rep)==true);
   a(:,2) = rep*ones(length(a),1);
   walkLoc = [walkLoc; a];
   clear a
end
%check that the walk length matches the expected length
if ~(size(walkLoc,1)==ntotal)
    warndlg('mismatched number of walking trials')
    return
end

for idx = 1:ntotal %1:ntotal % trial by trial
    cond = walkLoc(idx,1);
    rep = walkLoc(idx,2);
    disp(['begin ' num2str(idx) '/' num2str(ntotal)])
    % ---- save trial information ------
    trial(idx).ID = [ifly,cond,rep];
    % ----------------------------------

    % Filter out non-walkers
    if behavior(ifly).walking(cond,rep)==false
      continue; % skip other behaviors
    end

    % Pull the tarsus position data from the Pose 3d data files:
    min_vals = [];
    for leg = 1:6
        % find the raw data from the pose 3d
        pos = [leg_labels{leg} Alphabet(LP) '_x'];
        labels = pose_3d(ifly).labels{cond,rep};
        label_loc = find(strcmpi(pos, labels));
        raw = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);

        % save the raw data
        data(leg).raw = raw(controlROI,:);
        % find the leg minimum 
        data(leg).min = min(raw);
        % offset the individual leg data
        data(leg).ind = data(leg).raw-data(leg).min;

        % save leg min for group offset
        min_vals = [min_vals; min(raw)];
    end
    % find the group XYZ value offsets to move whole fly to positive range
    offset = min(min_vals);
    for leg = 1:6
        data(leg).input = data(leg).raw+abs(offset);
    end
    % % PLOT THE DATA FOR EACH TARSUS IN SPACE (2 figures)
    % [fig1, fig2] = DLC_plot_tarsus(data);

    [cloud_1, cloud_2, cloud_3]  = DLC_getStancePointsWALK(data, false);
    cloud = [cloud_1; cloud_2; cloud_3]; % triple up on 'good' data points
    
    % ---- save trial information ------
    trial(idx).cloud = cloud;
    trial(idx).data = data;
    trial(idx).offset = offset; % data shift to positive space
    % ----------------------------------

    %  Find the sphere of best fit: 

    % RMSE with minimum optimization function: 
    objective = @(XX) sqrt(mean((pdist2([XX(1),XX(2),XX(3)],cloud)-XX(4)).^2,2));
    x0 = [0,0,0,1.3]; % starting guess [center-x,center-y,center-z,radius]
    [A,b,Aeq,beq] = deal([]);   %empty values
    lb = [-inf,-inf,-inf,1.0];  %lower bounds
    ub = [inf,inf, 0, 2.0];     %upper bounds
    XX = fmincon(objective,x0,A,b,Aeq,beq,lb,ub);
    Center = [XX(1),XX(2),XX(3)];
    Radius = XX(4);
    disp(['Radius: ' num2str(Radius)])

    % disp number of points that hit the threshold:
    R = Radius*0.01;                % 3-percent threshold
    dist = pdist2(Center, cloud);  % find the euclidian distances
    err = (sum(dist>(Radius-R) & dist<(Radius+R),2)/length(cloud))*100; %percent of points 'on' the sphere
    disp(['Points on ball: ' num2str(err) '%'])

%     if err <70
%         warndlg('Less than 70% of points within the stance bounds')
%     end
    clear objective x0 A b Aeq beq lb ub

    % ---- save trial information ------
    trial(idx).XX = XX;
    trial(idx).Center = Center;
    trial(idx).Radius = Radius;
    trial(idx).onepercent_err = err;
    % ----------------------------------

    % % older verison: 
    % [center, fig, err] = DLC_MANUALsphereFit(X, r)

    % Swing - stance prediction and threshold setting:
    % Turn figures back on:
    set(0,'DefaultFigureVisible','on');
    % --- Figure of ball with full step cycles ---- 
    LW = 1.5;

    fig = getfig;
    subplot(1,2,1)
    for leg = 1:6
        input = data(leg).input;
        plot3(input(:,1), input(:,2), input(:,3), 'linewidth', LW); hold on
    end
    axis tight
    box off
    set(gca,'visible','off')
    [x,y,z] = sphere;
    s = surf(Radius*x+Center(1), Radius*y+Center(2), Radius*z+Center(3));
    set(s, 'FaceColor', Color('grey'))
    alpha 0.8
    axis vis3d % sets the aspect ratio for 3d rotation

    R = Radius*1.02; % increase threshold

    % find the distance from the center of the sphere to the tarsus and greater
    % than the radius value equates to swing:
    subplot(1,2,2)
    for leg = 1:6

        input = data(leg).input;
        e_dist = pdist2(Center,input);
        stance_loc = e_dist<=R;
        swing_loc = e_dist>R; 
        scatter3(input(stance_loc,1),input(stance_loc,2),input(stance_loc,3), 20, Color('grey'), 'filled')
        hold on
        scatter3(input(swing_loc,1),input(swing_loc,2),input(swing_loc,3), 100, 'p', 'filled')
    end


    uiwait(fig) % don't continue until the figure is closed

    answer = questdlg('Fit okay?');
    switch answer
        case 'Yes'
            trial(idx).fit = true;
        case 'No'
            trial(idx).fit = false;
    end
    close all %clears hidden figures generated by the DLC_getStancePointsWALK function
    clearvars('-except',var_list{:})
    idx = idx+1;

end % trial
disp('Done')

save([fig_dir '\' FilePath.locations{ifly,4} '\Walking ball fit'], 'trial')


angles(ifly).ballFit = trial;
% find the avg ball location across all trials? Use it for the best fit??
for idx = 1:ntotal
    if trial(idx).fit == false
        continue
    end
    % pull fly trial information:
    cond = trial(idx).ID(2);
    rep = trial(idx).ID(3);
    % save ball location information into the angles data structure
    offset = abs(trial(idx).offset);
    newCenter = trial(idx).Center-offset;
    angles(ifly).ball.Center{cond,rep} = newCenter;
    angles(ifly).ball.Radius(cond,rep) = trial(idx).Radius;
end


clearvars('-except',var_list{:})

%% Plot the swing-stance differentiation with various thresholds:
set(0,'DefaultFigureVisible','on');
thresh = 0.02; % percent
for idx = 1:ntotal
    cond = walkLoc(idx,1);
    rep = walkLoc(idx,2);
    Center = angles(ifly).ball.Center{cond,rep};
    Radius = angles(ifly).ball.Radius(cond,rep);
    R = Radius+Radius*thresh;
    if isempty(Center); continue; end
    % make figure
    fig = getfig;
    for leg = 1:6
        % find the raw data from the pose 3d
        pos = [leg_labels{leg} Alphabet(LP) '_x'];
        labels = pose_3d(ifly).labels{cond,rep};
        label_loc = find(strcmpi(pos, labels));
        input = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);
        e_dist = pdist2(Center,input);
        stance_loc = e_dist<=R;
        swing_loc = e_dist>R; 
        scatter3(input(stance_loc,1),input(stance_loc,2),input(stance_loc,3), 20, Color('grey'), 'filled')
        hold on
        scatter3(input(swing_loc,1),input(swing_loc,2),input(swing_loc,3), 100, 'p', 'filled')
    end
    xlabel('x')
    ylabel('y')
    zlabel('z')
end


%% Swing-Stance Plots: 
% inputs: data, Center, Radius, threshold, laser logical (optional)
% outputs: figure with swing/stance, logical of each leg and frame
% swing = black, stance = white
LP = 5; % tarsus selection
thresh = 0.02; % error threshold for 'ground'

for tt = 1:ntotal
    cond = walkLoc(tt,1);
    rep = walkLoc(tt,2);
    Center = angles(ifly).ball.Center{cond,rep};
    Radius = angles(ifly).ball.Radius(cond,rep);
    if isempty(Center); continue; end
    
    laserROI = [param.basler_delay*fps,ceil((param.conds_matrix(cond).opto+param.basler_delay)*fps)];

    % pull the tarsus position data: 
    for leg = 1:6
        % find the raw data from the pose 3d
        pos = [leg_labels{leg} Alphabet(LP) '_x'];
        labels = pose_3d(ifly).labels{cond,rep};
        label_loc = find(strcmpi(pos, labels));
        raw = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);

        % save the raw data
        data(leg).raw = raw;
    end

    [tFig, stance] = DLC_plot_swingNstance(data,Center,Radius,thresh);
    close(tFig)
    % quick spatial filter:
    n = 3; % filter step size
    for idx = 1:2 % catch any that don't smooth first time
    for leg = 1:6
        input = stance(leg,:);
        a = smooth(abs(diff(input)),n);
        b = find(a>(1/n)); % more than one transition within the n range
        if ~isempty(b)
          for i = 1:length(b)
            ii = b(i);
            if ii<=n || abs(ii-length(stance))<=n; continue; end
            c = input(ii-n:ii);
            if sum(c)==0||sum(c)==(n+1); continue;
            else
                input(ii) = c(1);
            end
          end
        end
        stance(leg,:) = input;
        % check for holes of 2 frames:
    end
    end

    fig = DLC_plot_swingNstance(stance); % view the updated swing-stance plot
    hold on;
    try vline(laserROI, 'g-');catch;end
    xlabel('time (fps)')
    set(gca,'YTickLabel',[]);
    ylabel('leg')
    fig_title = [FilePath.locations{ifly,4}, ' R' num2str(rep), 'C' num2str(cond)];
    fig_title = strrep(fig_title, '_', '-');
    title({FilePath.structure_name; fig_title})

    % need to add a quick filter for points that skip threshold
    angles(ifly).stance(cond,rep).loc = stance;
    fig_path = [fig_dir, '\' FilePath.locations{ifly,4}]; 
    save_figure(fig, [fig_path, '\R' num2str(rep), 'C' num2str(cond) ' sw-st plot']);
end
% want to save these figures into the folder for that specific fly...

%% Code to save a 3d rotating video of a figure:



ViewZ = [10,60; 50,60; 100,30; 180,30; 280, 10; 360, 10; 410, 10];
VidName = 'C:\Users\evynd\Desktop\3D Rotating Video';
OptionZ.FrameRate=10;OptionZ.Duration=10;OptionZ.Periodic=true;
CaptureFigVid(ViewZ, VidName, OptionZ) 


%% load formerly save fly info:

% Categorize and save plot data to Angles from previously saved sphere fitting data
% select fly & load presaved data if not already in Workspace
if strcmpi(questdlg('Load old data?'), 'Yes')==true
    [file,path] = uigetfile;
    load(fullfile(path, file))
end

ifly = 1;

ntotal = sum(sum(behavior(ifly).walking)); % total walking trials
% pull the list of walking: 
walkLoc = [];
for rep = 1:num.reps
   a(:,1) = find(behavior(ifly).walking(:,rep)==true);
   a(:,2) = rep*ones(length(a),1);
   walkLoc = [walkLoc; a];
   clear a
end
%check that the walk length matches the expected length
if ~(size(walkLoc,1)==ntotal)
    warndlg('mismatched number of walking trials')
    return
end

% Save the data into the 'angles' data structure:
angles(ifly).ballFit = trial;
for idx = 1:ntotal
    if trial(idx).fit == false
        continue
    end
    % pull fly trial information:
    cond = trial(idx).ID(2);
    rep = trial(idx).ID(3);
    % save ball location information into the angles data structure
    offset = abs(trial(idx).offset);
    newCenter = trial(idx).Center-offset;
    angles(ifly).ball.Center{cond,rep} = newCenter;
    angles(ifly).ball.Radius(cond,rep) = trial(idx).Radius;
end

% Find the 'stance' status for each frame: 
LP = 5; % tarsus selection
thresh = 0.02; % error threshold for 'ground'
for idx = 1:ntotal
    cond = walkLoc(idx,1);
    rep = walkLoc(idx,2);
    Center = angles(ifly).ball.Center{cond,rep};
    Radius = angles(ifly).ball.Radius(cond,rep);
    if isempty(Center); continue; end% check for failed ball fit

    % pull the tarsus position data: 
    for leg = 1:6
        % find the raw data from the pose 3d
        pos = [leg_labels{leg} Alphabet(LP) '_x'];
        labels = pose_3d(ifly).labels{cond,rep};
        label_loc = find(strcmpi(pos, labels));
        raw = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);
        % save the raw data
        data(leg).raw = raw;
    end

    [tFig, stance] = DLC_plot_swingNstance(data,Center,Radius,thresh);
    close(tFig)
    
    % quick spatial filter:
    n = 3; % filter step size
    for tt = 1:2 % catch any that don't smooth first time
    for leg = 1:6
        input = stance(leg,:);
        a = smooth(abs(diff(input)),n);
        b = find(a>(1/n)); % more than one transition within the n range
        if ~isempty(b)
          for i = 1:length(b)
            ii = b(i);
            if ii<=n || abs(ii-length(stance))<=n; continue; end
            c = input(ii-n:ii);
            if sum(c)==0||sum(c)==(n+1); continue;
            else
                input(ii) = c(1);
            end
          end
        end
        stance(leg,:) = input;
        % check for holes of 2 frames:
    end
    end

    angles(ifly).stance(cond,rep).loc = stance;
end



clearvars('-except',initial_vars{:})


%% %% FLY BALL AVERAGE: ** doesn't work well -- too many small shifts by vid
% (find the avg ball position across a series of trials -- doesn't fit well
% -- too much movement in the ball position from video to video tracking
% 
% 
% ntotal = sum(sum(behavior(ifly).walking)); % total walking trials
% % pull the list of walking: 
% walkLoc = [];
% for rep = 1:num.reps
%    a(:,1) = find(behavior(ifly).walking(:,rep)==true);
%    a(:,2) = rep*ones(length(a),1);
%    walkLoc = [walkLoc; a];
%    clear a
% end
% %check that the walk length matches the expected length
% if ~(size(walkLoc,1)==ntotal)
%     warndlg('mismatched number of walking trials')
%     return
% end
% 
% for idx = 1:ntotal
%     cond = walkLoc(idx,1);
%     rep = walkLoc(idx,2);
%     disp(['begin ' num2str(idx) '/' num2str(ntotal)])
%     % ---- save trial information ------
%     trial(idx).ID = [ifly,cond,rep];
%     % ----------------------------------
% 
%     % Filter out non-walkers
%     if behavior(ifly).walking(cond,rep)==false
%       continue; % skip other behaviors
%     end
% 
%     % Pull the tarsus position data from the Pose 3d data files:
%     min_vals = [];
%     for leg = 1:6
%         % find the raw data from the pose 3d
%         pos = [leg_labels{leg} Alphabet(LP) '_x'];
%         labels = pose_3d(ifly).labels{cond,rep};
%         label_loc = find(strcmpi(pos, labels));
%         raw = pose_3d(ifly).raw{cond,rep}(:,label_loc:label_loc+2);
% 
%         % save the raw data
%         data(leg).raw = raw(controlROI,:);
%         % find the leg minimum 
%         data(leg).min = min(raw);
%         % offset the individual leg data
%         data(leg).ind = data(leg).raw-data(leg).min;
% 
%         % save leg min for group offset
%         min_vals = [min_vals; min(raw)];
%     end
%     % find the group XYZ value offsets to move whole fly to positive range
%     offset = min(min_vals);
%     for leg = 1:6
%         data(leg).input = data(leg).raw+abs(offset);
%     end
%     % % PLOT THE DATA FOR EACH TARSUS IN SPACE (2 figures)
%     % [fig1, fig2] = DLC_plot_tarsus(data);
% 
%     [cloud_1, cloud_2, cloud_3]  = DLC_getStancePointsWALK(data, false);
% %     cloud = [cloud_1; cloud_2; cloud_3]; % triple up on 'good' data points
%     
%     % ---- save trial information ------
%     trial(idx).cloud_1 = cloud_1;
%     trial(idx).cloud_2 = cloud_2;
%     trial(idx).cloud_3 = cloud_3;
%     trial(idx).data = data;
%     trial(idx).offset = offset; % data shift to positive space
% 
%     close all %clears hidden figures generated by the DLC_getStancePointsWALK function
%     clearvars('-except',var_list{:})
%     idx = idx+1;
% 
% end % trial
% set(0,'DefaultFigureVisible','on');
% disp('Done')
% 
% 
% %% Pool all the data clouds & plot the avg ball: 
% % find the avg ball location across all trials? Use it for the best fit??
% [cloud_1,cloud_2,cloud_3] = deal([]);
% for idx = 1:ntotal
%     cloud_1 = [cloud_1; trial(idx).cloud_1];
%     cloud_2 = [cloud_2; trial(idx).cloud_2];
%     cloud_3 = [cloud_3; trial(idx).cloud_3];
% end
% cloud = [cloud_1; cloud_2; cloud_3];
% 
% 
% % RMSE with minimum optimization function: 
% objective = @(XX) sqrt(mean((pdist2([XX(1),XX(2),XX(3)],cloud)-XX(4)).^2,2));
% x0 = [0,0,0,1.3]; % starting guess [center-x,center-y,center-z,radius]
% [A,b,Aeq,beq] = deal([]);   %empty values
% lb = [-inf,-inf,-inf,1.0];  %lower bounds
% ub = [inf,inf, 0, 2.0];     %upper bounds
% XX = fmincon(objective,x0,A,b,Aeq,beq,lb,ub);
% Center = [XX(1),XX(2),XX(3)];
% Radius = XX(4);
% disp(['Radius: ' num2str(Radius)])
% % 
% % disp number of points that hit the threshold:
% R = Radius*0.03;                % 5-percent threshold
% dist = pdist2(Center, cloud);  % find the euclidian distances
% err = (sum(dist>(Radius-R) & dist<(Radius+R),2)/length(cloud))*100; %percent of points 'on' the sphere
% disp(['Points on ball: ' num2str(err) '%'])
% 
% 
% fig = getfig('',1);
% % scatter3(cloud(:,1), cloud(:,2), cloud(:,3), 20, Color('teal'), 'filled')
% % color code 'on ball vs off ball'
% % ON = teal, OFF = red
% e_dist = pdist2(Center,cloud);
% stance_loc = e_dist<=Radius+R;
% swing_loc = e_dist>Radius+R; 
% scatter3(cloud(stance_loc,1),cloud(stance_loc,2),cloud(stance_loc,3), 20, Color('teal'), 'filled')
% hold on
% scatter3(cloud(swing_loc,1),cloud(swing_loc,2),cloud(swing_loc,3), 20, Color('darkred'), 'filled')
% 
% hold on
% scatter3(Center(1), Center(2), Center(3), 200, 'r', 'filled')
% axis tight
% box off
% set(gca,'visible','off')
% [x,y,z] = sphere;
% s = surf(Radius*x+Center(1), Radius*y+Center(2), Radius*z+Center(3));
% set(s, 'FaceColor', Color('grey'))
% alpha 0.8
% axis vis3d % sets the aspect ratio for 3d rotation
% title(['Point Cluster: ' num2str(ii)])
% 
% 
% 
% % color code 'on ball vs off ball'
% % ON = teal, OFF = red
% e_dist = pdist2(Center,cloud);
% stance_loc = e_dist<=R;
% swing_loc = e_dist>R; 
% scatter3(cloud(stance_loc,1),cloud(stance_loc,2),cloud(stance_loc,3), 20, Color('teal'), 'filled')
% hold on
% scatter3(cloud(swing_loc,1),cloud(swing_loc,2),cloud(swing_loc,3), 20, Color('darkred'), 'filled')

















