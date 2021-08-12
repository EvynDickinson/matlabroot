

% DLC_Step2_basicfigures
% single fly data processing -- requires completion of step1
% 1) predict the fly behavior from FeTi angle of each leg
% 2) fit a sphere to the tarsus position points
% 3) calculate swing or stance for all legs in each video frame
% 4) save the data to 'StepData' in the local and/or online folders



clearvars
clc
close all
warning('off')


% TODO section:
% Set base file paths:
fileroot = 'D:\Tuthill Lab Work\Anipose\';                  %local location for step1 processed Anipose data
googledrive = 'G:\My Drive\Evyn\Data\FicTrac Raw Data\';     %google drive root folder
% ----------------------
folder_date = '6.12.19';                                    %folder to analyze
% ----------------------

local_save = true;                                          %save data in the local PC folder?
online_save = true;                                         %save data to google drive folder?

%% Select fly within the folder_date to load:
pathway = [fileroot folder_date]; % directory
list_flies = dir(pathway);
idx = 1;
for ii = 1:length(list_flies)
    a = list_flies(ii).name; %names of the flies for the day
    if strfind(a,'Fly')==1
        fly_num{idx} = a;
        idx = idx+1;
    else; continue;
    end
end; clear idx
%Select the flies to analyze
idx = listdlg('ListString', fly_num, 'SelectionMode', 'Single');
FlyToLoad = fly_num(idx);
clear a ii idx list_flies

% Load data from the selected fly:
flyID = generate_flyID(folder_date, FlyToLoad{1}(end-2:end));
try
    load([pathway,'\' FlyToLoad{1} '\Analysis\' flyID '_AniposeData.mat'])
catch
    warndlg('No file matching AniposeData found')
    return
end

fig_root = [pathway,'\' FlyToLoad{1} '\Analysis\' flyID];

disp('Fly selected: ')
disp(FlyToLoad)

%% Behavior Prediction:
angles = anipose.angles;
param = anipose.param;
fps = param.Basler_fps;

% make leg behavior and fly behavior predictions
ROI = 1:param.basler_delay*fps; % sort behavior by all frames before the laser
[behavior, figHandle] = DLC_behaviorpredictions(angles,param,ROI);
save_figure(figHandle, [fig_root, 'Behavior Predictions'], '-png');



disp('Finished predicting behavior!')

%% Sphere Fit all possible videos: 
% for each trial, sort by behavior and use the specific function for
% isolating points in 'stance' to use for sphere fitting
pose_3d = anipose.pose_3d;
num = NUM(param);
LP = 5; %leg point = tarsus
rMin = 1.0;     % min radius constraint
rMax = 3.0;     % max radius constraint
thresh = 0.01;  % percent threshold for 'ground'
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
for leg = 1:6
   kolor(leg,:) = Color(leg_colors{leg}); 
end
SZ = 20;        % size of scatter plot points
LW = 1;         % line width for plots

for cond = 1:num.conds
    for rep = 1:num.reps
        % ---- set up blank trial information ------
        trial(cond,rep).cloud = [];
        trial(cond,rep).Center = [];
        trial(cond,rep).Radius = [];
        step.ball.Center{cond,rep} = [];
        step.ball.Radius(cond,rep) = nan;
        
        % pull the data during the control ROI:
        for leg = 1:6
            % find the raw data from the pose 3d
            raw = pose_3d.Leg(LP,leg).data{cond,rep}.all;
            % save the raw data
            data(leg).raw = raw(ROI,:);
            % find the leg minimum 
            data(leg).min = min(raw);
            % offset the individual leg data
            data(leg).ind = data(leg).raw-data(leg).min;
            % save leg min for group offset
            min_vals(leg,:) = min(raw);
        end
        % find the group XYZ value offsets to move whole fly to positive range
        offset = min(min_vals);
        for leg = 1:6
            data(leg).input = data(leg).raw+abs(offset);
        end
        
        % % PLOT THE DATA FOR EACH TARSUS IN SPACE (2 figures)
        % [fig1, fig2] = DLC_plot_tarsus(data);
        trial(cond,rep).data = data;
        trial(cond,rep).offset = offset; % data shift to positive space
        
        % Get a cloud of points to fit the sphere
        switch behavior.class(cond,rep) % fly behavior:
            case 1 %walking
                [cloud_1, cloud_2, cloud_3]  = DLC_getStancePointsWALK(data, false);
                cloud = [cloud_1; cloud_2; cloud_3]; % triple up on 'good' data points
            case 2 %stationary
                % TODO
                continue;
            case 3 %grooming
                % TODO
                continue;
            case 4 %other
                cloud = []; % no process for point selection yet....
                continue;
        end

        % ---- save trial information ------
        trial(cond,rep).cloud = cloud;
        
        % SPHERE FIT
        %RMSE with minimum optimization function: 
        objective = @(XX) sqrt(mean((pdist2([XX(1),XX(2),XX(3)],cloud)-XX(4)).^2,2));
        x0 = [0,0,0,1.3]; % starting guess [center-x,center-y,center-z,radius]
        [A,b,Aeq,beq] = deal([]);   %empty values
        lb = [-inf,-inf,-inf,rMin];  %lower bounds
        ub = [inf,inf, 0, rMax];     %upper bounds
        XX = fmincon(objective,x0,A,b,Aeq,beq,lb,ub);
        Center = [XX(1),XX(2),XX(3)];
        Radius = XX(4);
        disp(['Radius: ' num2str(Radius)])

        % disp number of points that hit the threshold:
        R = Radius*thresh;             % 3-percent threshold
        dist = pdist2(Center, cloud);  % find the euclidian distances
        err = (sum(dist>(Radius-R) & dist<(Radius+R),2)/length(cloud))*100; %percent of points 'on' the sphere
        disp(['Points on ball: ' num2str(err) '%'])
        
        trial(cond,rep).cloud = cloud;
        trial(cond,rep).Center = Center;
        trial(cond,rep).Radius = Radius;

        set(0,'DefaultFigureVisible','on');  %Turn figures back on:
        % --- Figure of ball with full step cycles ---- 
        LW = 1.5;
        
        fig = getfig;
        subplot(1,2,1)
        for leg = 1:6
            input = data(leg).input;
            plot3(input(:,1), input(:,2), input(:,3), 'linewidth', LW, 'color', kolor(leg,:)); hold on
        end
        axis tight
        box off
        set(gca,'visible','off')
        [x,y,z] = sphere;
        s = surf(Radius*x+Center(1), Radius*y+Center(2), Radius*z+Center(3));
        set(s, 'FaceColor', Color('grey'))
        alpha 0.8
        axis vis3d % sets the aspect ratio for 3d rotation

        R = Radius*(1+thresh); % increase threshold
        % find the distance from the center of the sphere to the tarsus and greater
        % than the radius value equates to swing:
        subplot(1,2,2)
        for leg = 1:6
            input = data(leg).input;
            e_dist = pdist2(Center,input);
            stance_loc = e_dist<=R;
            swing_loc = e_dist>R; 
            scatter3(input(stance_loc,1),input(stance_loc,2),input(stance_loc,3), SZ, Color('grey'), 'filled')
            hold on
            scatter3(input(swing_loc,1),input(swing_loc,2),input(swing_loc,3), SZ, kolor(leg,:), 'filled')
        end
        title([strrep(flyID,'_','-') ' Rep: ' num2str(rep) ' Cond: ' num2str(cond)])
        % don't continue until 'ENTER' is pressed:
        currkey = 0;
        while currkey~=1
            pause; % wait for a keypress
            currkey=get(gcf,'CurrentKey'); 
            if strcmpi(currkey, 'return') %|| ~ishghandle(fig)
                currkey = 1;
            else
                currkey = 0;
            end
        end
        close(fig)
        
        answer = questdlg('Fit okay?');      
        switch answer
            case 'Yes'
                trial(cond,rep).fit = true;
                newCenter = Center-abs(offset);
                step.ball.Center{cond,rep} = newCenter;
                step.ball.Radius(cond,rep) = Radius;
            case 'No'
                trial(cond,rep).fit = false;
            case 'Cancel'
                return  % breaks the loop
                
        end
        % Clean up
        close all %clears hidden figures generated by functions
        clear lb ub XX Center Radius x0 objective A b Aeq beq err dist cloud
        clear swing_loc stance_loc e_dist R s x y z data min_vals 
        clear cloud_1 cloud_2 cloud_3 newCenter offset raw
        disp(['logged R' num2str(rep) 'C' num2str(cond)])
    end
end

% immediately save the data! This ^^ is a long process and you don't want
% to do it more than once if possible! 
save([fig_root ' Sphere Fit data'], 'trial');
fprintf('\n\nFinished fitting spheres')

% save the big information:
step.ballFit = trial;



%% OPTIONAL viewing tool for determining the best threshold: 
%if the threshold didn't look good but fit was generally good
answer = questdlg('Check threshold?','Threshold 1');
if strcmpi(answer,'yes')==true
    
    thresh = str2double(inputdlg('Input the best threshold for this fly'));
    idxMax = str2double(inputdlg('Max number of figures to preview?'));
    idx = 0;
    
    set(0,'DefaultFigureVisible','on');
    for cond = 1:num.conds
      for rep = 1:num.reps
        Center = step.ball.Center{cond,rep};
        Radius = step.ball.Radius(cond,rep);
        R = Radius+Radius*thresh;
        if isempty(Center); continue;
        end
        
        if idx>idxMax; continue; end
        
        % PLOT the figure
        fig = getfig('',1);
        for leg = 1:6
            % find the raw data from the pose 3d
            input = pose_3d.Leg(LP,leg).data{cond,rep}.all;
            e_dist = pdist2(Center,input);
            stance_loc = e_dist<=R;
            swing_loc = e_dist>R; 
            scatter3(input(stance_loc,1),input(stance_loc,2),input(stance_loc,3), SZ, Color('grey'), 'filled')
            hold on
            scatter3(input(swing_loc,1),input(swing_loc,2),input(swing_loc,3), SZ, kolor(leg,:), 'filled')
        end
        xlabel('x')
        ylabel('y')
        zlabel('z')
        last_cond = cond;
        last_rep = rep;
        idx = idx+1;
        % wait to open the next figure until 'ENTER' is pressed:
        currkey = 0;
        while currkey~=1
            pause; % wait for a keypress
            currkey=get(gcf,'CurrentKey'); 
            if strcmpi(currkey, 'return')
                currkey = 1;
            else
                currkey = 0;
            end
        end
      end
    end
    % Additional threshold adjustment check
    answer = questdlg('Adjust threshold further?','Threshold 2','Yes','No','Cancel','No');
    if strcmpi(answer,'yes')==true
        thresh = str2double(inputdlg('Input the best threshold for this fly'));
        % plot final figure again with new threshold:
        cond = last_cond; rep = last_rep;
        Center = step.ball.Center{cond,rep};
        Radius = step.ball.Radius(cond,rep);
        R = Radius+Radius*thresh;
        fig = getfig;
        for leg = 1:6
            % find the raw data from the pose 3d
            input = pose_3d.Leg(LP,leg).data{cond,rep}.all;
            e_dist = pdist2(Center,input);
            stance_loc = e_dist<=R;
            swing_loc = e_dist>R; 
            scatter3(input(stance_loc,1),input(stance_loc,2),input(stance_loc,3), SZ, Color('grey'), 'filled')
            hold on
            scatter3(input(swing_loc,1),input(swing_loc,2),input(swing_loc,3), SZ, kolor(leg,:), 'filled')
        end
        xlabel('x')
        ylabel('y')
        zlabel('z')
        %wait for fig to close:
        uiwait(fig)
        answer = questdlg('Adjust threshold further?','Threshold 3','Yes','No','Cancel','No');
        if strcmpi(answer,'yes')==true
            thresh = str2double(inputdlg('Input the best threshold for this fly'));
        end
    elseif strcmpi(answer,'Cancel')==true %abort option
        return
    end
    close all
    disp(['Current threshold: ' num2str(thresh)])
elseif strcmpi(answer,'Cancel')==true
    return % abort program -- stops it. 'no' just skips this section
end

%% Swing-Stance Plots -- required section to categorized each section as swing/stance moving forward
% inputs: data, Center, Radius, threshold, laser logical (optional)
% outputs: figure with swing/stance, logical of each leg and frame
% swing = black, stance = white
LP = 5; % tarsus selection
fig_dir = [fig_root ' swing-stance plots\'];
if ~exist(fig_dir)
    mkdir(fig_dir)
end
disp('Extracting swing-stance locations now')
% Extract and save the step swing/stance information:
for cond = 1:num.conds
  for rep = 1:num.reps
    Center = step.ball.Center{cond,rep};
    Radius = step.ball.Radius(cond,rep);
    step.stance(cond,rep).loc = [];
    if isempty(Center); continue; end
    % pull the tarsus position data: 
    for leg = 1:6
         data(leg).raw = pose_3d.Leg(LP,leg).data{cond,rep}.all;
    end 
    [fig, stance] = DLC_plot_swingNstance(data,Center,Radius,thresh);
    close(fig)
    % quick spatial filter:
    n = 5; % filter step size
    for idx = 1:n-1 % catch any that don't smooth first time
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
    step.stance(cond,rep).loc = stance;
  end
end


answer = questdlg('View swing|stance plots?');
if strcmpi(answer,'yes')==true
    idxMax = str2double(inputdlg('Max number of figures to preview?'));
    idx = 0;
    % PLOT the swing-stance figures: 
    for rep = 1:num.reps
      for cond = 1:num.conds
        % pull information
        Center = step.ball.Center{cond,rep};
        Radius = step.ball.Radius(cond,rep);
        stance = step.stance(cond,rep).loc;
        if isempty(Center); continue; end
        idx = idx+1;
        if idx>idxMax; continue; end
        laserROI = [param.basler_delay*fps,ceil((param.conds_matrix(cond).opto+param.basler_delay)*fps)];
        % figure
        fig = DLC_plot_swingNstance(stance); % view the updated swing-stance plot
        hold on;
        try vline(laserROI, 'g-');catch;end
        xlabel('time (fps)')
        set(gca,'YTickLabel',[]);
        ylabel('leg')
        fig_title = [flyID, ' R' num2str(rep), 'C' num2str(cond)];
        fig_title = strrep(fig_title, '_', '-');
        title({param.cross; fig_title})
        results = save_figure(fig, [fig_dir, 'R' num2str(rep), 'C' num2str(cond) ' sw-st plot']);
        if ~islogical(results) %press 'Cancel' to quit the loop
            return
        end
      end
    end
elseif strcmpi(answer,'Cancel')==true
    return % abort the program
end

disp('Finished fitting steps!')


%% Save all the progress made! 
% save all updated information to 'StepData'

if local_save == true
   save([fig_root ' StepData'], 'step');
end
if online_save == true
    Gdir = [googledrive, folder_date,'\',FlyToLoad{1},'\Analysis'];
    save([Gdir,'\', flyID ' StepData'], 'step')
end   
disp(flyID)
fprintf('\n\n---------------------------\n Saved fly data! All done!\n---------------------------\n')



















