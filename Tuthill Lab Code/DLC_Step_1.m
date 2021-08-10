
% Step one analysis for DLC & anipose analyzed videos
% load and reorganize the angles data to fit the Fictrac data structure
clearvars
clc
close all

%% Load Anipose Data:

a = questdlg('Open previously generated file?');

switch a
    case 'Yes'
        [file,path] = uigetfile;
        load(fullfile(path, file))
        % make a folder for the analysis:
        fig_dir = [output_root, FilePath.structure_name];
        if ~exist(fig_dir, 'dir')
            mkdir(fig_dir)
        end 
    case 'No'
        % rootpath = 'C:\Users\evynd\Desktop\HeadlessFliesDLC\tracking\off_ball-flies\training_videos\';
        % Set the base directory for the data files: 
        fileroot = 'G:\My Drive\Evyn\Data\FicTrac Raw Data\'; %TODO
        output_root = 'C:\matlabroot\DLC\'; %TODO
        FilePath = DLC_select_flies('-inpath', fileroot, '-outpath', output_root, '-export', true);
        % load formerly-generated structures: 

        % make a folder for the analysis:
        fig_dir = [output_root, FilePath.structure_name];
        if ~exist(fig_dir, 'dir')
            mkdir(fig_dir)
        end 
        % fig_dir = [fig_dir, '\' FilePath.structure_name];

        % load structure information (e.g. num.conds, num.reps, etc):
        load([FilePath(1).locations{1,1}, FilePath(1).locations{1,2}, '\' ...
              FilePath(1).locations{1,3}, '/FicTrac Data/Parameters.mat']); 
        num = NUM(param);
        num.flies = length(FilePath.locations);

        % ---- Load the angles data ---- 
        %SLOW -- if your files are on GoogleDrive
        angles = struct; fprintf('\n Finished fly: \n')
        for ifly = 1:num.flies
          tic
          for cond = 1:num.conds
              for rep = 1:num.reps
                metaData = [];
                FolderName = [FilePath.locations{ifly,1}, FilePath.locations{ifly,2}...
                                     '\', FilePath.locations{ifly,3}, '\angles'];

                prep = [FilePath.locations{ifly,4} ' R' num2str(rep) 'C' num2str(cond)...
                        '  ' param.conds_matrix(cond).label, '.csv'];

                csvFile = fullfile(FolderName, prep);
                metaData = readtable(csvFile);
                angles(ifly).labels{cond,rep} = metaData.Properties.VariableNames;
                angles(ifly).raw{cond,rep} = table2array(metaData);
              end
          end
          disp(ifly)
          toc
        end

%         % save struct quickly to not have to sit through that again!
%         save([fig_dir, FilePath.structure_name, ' leg angles'])

        %  ---- Load Joint positions ---- 
        pose_3d = struct; fprintf('\n Finished fly: \n')
        for ifly = 8:num.flies
          tic
          for cond = 1:num.conds
              for rep = 1:num.reps
                metaData = [];
                FolderName = [FilePath.locations{ifly,1}, FilePath.locations{ifly,2}...
                                     '\', FilePath.locations{ifly,3}, '\pose-3d'];

                prep = [FilePath.locations{ifly,4} ' R' num2str(rep) 'C' num2str(cond)...
                        '  ' param.conds_matrix(cond).label, '.csv'];

                csvFile = fullfile(FolderName, prep);
                metaData = readtable(csvFile);
                pose_3d(ifly).labels{cond,rep} = metaData.Properties.VariableNames;
                pose_3d(ifly).raw{cond,rep} = table2array(metaData);
              end
          end
          disp(ifly)
          toc
        end

        save([fig_dir, FilePath.structure_name, ' pose-3d'])
end


% save the names of the initial variables
initial_vars = who; initial_vars{end+1} = 'initial_vars';
% clearvars('-except',initial_vars{:})

%% Reorganize the Angles data and save into structures:
% NOTE: data output from Anipose varies depending on the version of anipose
% some currently have 19 data columns and some have 73 -- autoselect added
fps = param.Basler_fps;
num.frames = param.basler_length*fps+1;
% labels list for visual check (particularly important for loading in new
% data to this:
Joints = {'CoFe', 'FeTi', 'TiTa'};
num.joints = length(Joints);

all_angles = {'L1_CF','L1_FTi','L1_TiTa',...
              'L2_CF','L2_FTi','L2_TiTa',...
              'L3_CF','L3_FTi','L3_TiTa',...
              'R1_CF','R1_FTi','R1_TiTa',...
              'R2_CF','R2_FTi','R2_TiTa',...
              'R3_CF','R3_FTi','R3_TiTa'};
clear input
% ----Angles----
for ifly = 1:num.flies
  % find the correct extraction data columns:
  input(:,1) = angles(ifly).labels{1,1};
  %find the data columns for all angles
  for ii = 1:length(all_angles)
     angle_loc(ii) = find(strcmpi(all_angles(ii), input));
     label_check(ii,1) = input(angle_loc(ii));
  end
  angles(ifly).angle_loc = angle_loc; % from above
  
  %visually compare the joint angle labels:
%   fprintf('\n Please check that the labels align correctly: \n')
%   disp([' Fly: ' FilePath.locations{ifly,4}])
%   disp([all_angles', label_check])
%   
  for cond = 1:num.conds
     for rep = 1:num.reps
       raw = angles(ifly).raw{cond,rep};
       % check that the number of frames is what was expected|replace with
       % nans
       if ~(size(raw,1)==num.frames) %frame numbers don't match=true
           raw = nan(num.frames,size(raw,2));
           fprintf(['\n Missing frames found in ' FilePath.locations{ifly,4} ...
                    ' R' num2str(rep) 'C' num2str(cond) '\n'])
       end
       for leg = 1:6 %six legs on the fly
           ROI = leg*num.joints-(num.joints-1):leg*num.joints;
           loc = angle_loc(ROI);
           %extract the time series data:
           angles(ifly).Leg(leg).data(cond,rep).all = raw(:,loc); 
           a = angles(ifly).labels{cond,rep}(loc);
           angles(ifly).Leg(leg).labels = strrep(a,'_','-'); %swap underscores to dashes
           angles(ifly).Leg(leg).loc = loc;
       end
     end
  end
  clear input
end

clear a ans csvFile file ifly ii leg  loc label_check metaData prep raw rep ROI

% add more fixed variables:
x = -(param.basler_delay):1/fps:param.basler_length-(param.basler_delay);
leg_labels = {'L1', 'L2', 'L3', 'R1', 'R2', 'R3'};


initial_vars = who; initial_vars{end+1} = 'initial_vars';

initial_vars{end+1} = 'group';
initial_vars{end+1} = 'behavior';
initial_vars{end+1} = 'gait';



%% Numerous types of plots below: 



%% Plot single trial: all legs and all three joints
% e.g. 6 suplots : 1/leg
nrows = 3; % 1/leg
ncols = 2; %R||L legs
joint_colors = {'purple', 'darkgoldenrod', 'teal'};
LW = 1;
CCW_offset = num.conds/4;

% ---- INPUT ------
ifly = 3;
rep = 3;
cond = 7;
% -----------------

light_on = param.conds_matrix(cond).opto;
subplot_idx = [1,3,5,2,4,6];
fig = getfig('',1);
for leg = 1:6
    input = [];
    subplot(nrows,ncols,subplot_idx(leg)); hold all
    % extract the joint angle information:
    loc = angles(ifly).Leg(leg).loc;
    input = angles(ifly).raw{cond,rep}(:,loc); %raw data for that leg
    
    for iJ = 1:num.joints %plot all three joints
        y = input(:,iJ);
        plot(x, y, 'color', Color(joint_colors{iJ}), 'linewidth', LW) 
    end
    ylim([0, 200])
    yL = rangeLine(fig); %find max yvalue and get a value back 5% lower
    plot([0, light_on], [yL,yL], 'linewidth', LW+2, 'color', Color(param.LED_light))
    vline(0, 'k:')
    % labels: 
    ylabel(['Leg ' leg_labels{leg} ' angle (\circ)'])
    if subplot_idx(leg) == nrows*ncols-1; xlabel('time (s)'); 
    elseif subplot_idx(leg) == nrows*ncols; xlabel('time (s)');
    end
    fig_title = [FilePath.locations{ifly,4} ' R' num2str(rep), 'C' num2str(cond)];
    fig_title = strrep(fig_title, '_', '-');
    if subplot_idx(leg) == 1; title({fig_title; param.conds_matrix(cond).label}); end
  
end

save_figure(fig, [fig_dir, '\' FilePath.structure_name ' compete 3d ' fig_title]);

clearvars('-except',initial_vars{:})

%% Single condition: all reps overlaid for each joint
% variables: 
t_delay = param.basler_delay*fps;
leg = 1; %1-LF, 2-LM, 3-LR, 4-RF, 5-RM, 6-RR
Joints = {'CoFe', 'FeTi', 'TiTa'};
CCW_offset = num.conds/4;
LW = 1; %linewidth

% ---- INPUT ------
ifly = 1; 
cond = 7;
% -----------------

light_on = param.conds_matrix(cond).opto;

% Plot the data: 
fig = getfig('',1);
for iJ = 1:num.joints %number of joints
    [CW_input, CCW_input] = deal([]);
    subplot(num.joints,1,iJ)
    hold all
    for rep = 1:num.reps
        % plot all CW trials of the condition
        CW_input = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
        CCW_input = angles(ifly).Leg(leg).data(cond+CCW_offset,rep).all(:,iJ);
        plot(x, CW_input, 'color', Color('teal'))
        plot(x, CCW_input, 'color', Color('maroon'))
    end
    ylabel({'Joint angle'; angles(1).Leg(1).labels{iJ}})  
    % plot the opto_stim
    yl = max(ylim); %y axis max number 
    plot([0, light_on], [yl,yl], 'linewidth', LW+2, 'color', Color(param.LED_light))
    fig_name = [FilePath.locations{ifly,4} ' ' param.conds_matrix(cond).label];
    fig_name = strrep(fig_name, '_', '-');
    fig_name = strrep(fig_name, 'cw-', '');
    if iJ == 1; title(fig_name); end
end
xlabel('time (s)')


save_figure(fig, [fig_dir, '\rep overlay ' fig_name]);

clearvars('-except',initial_vars{:})

%% Plot by condition : all reps have subplot with all three joint angles
% e.g. 6 suplots : 1/rep
[nrows, ncols] = subplot_numbers(num.reps*2,2); %two max cols for the CW&CCW
joint_colors = {'purple', 'darkgoldenrod', 'teal'};
LW = 1;
CCW_offset = num.conds/4;

% ---- INPUT ------
ifly = 3;
leg = 1;
cond = 10;
% -----------------

light_on = param.conds_matrix(cond).opto;

fig = getfig('',1);
idx = 1;
for rep = 1:num.reps
  %CW stimulus:
    subplot(nrows,ncols,idx); hold all
    % extract the joint angle information:
    for iJ = 1:num.joints %plot all three joints
        y = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
        plot(x, y, 'color', Color(joint_colors{iJ}), 'linewidth', LW) 
    end
    ylim([0, 200])
    yL = rangeLine(fig);
    plot([0, light_on], [yL,yL], 'linewidth', LW+2, 'color', Color(param.LED_light))
    vline(0, 'k:')
    % labels: 
    ylabel(['Joint angle, rep ' num2str(rep)])
    if rep == num.reps; xlabel('time (s)'); end
    if rep == 1
        fig_name = [FilePath.locations{ifly,4}, ' ' param.conds_matrix(cond).label];
        fig_name = strrep(fig_name, '_', '-');
        title(fig_name); 
    end
    idx = idx+1;
    
    
  % CCW stimulus:
    subplot(nrows,ncols,idx); hold all
    % extract the joint angle information:
    for iJ = 1:num.joints %plot all three joints
        y = angles(ifly).Leg(leg).data(cond+CCW_offset,rep).all(:,iJ);
        plot(x, y, 'color', Color(joint_colors{iJ}), 'linewidth', LW) 
    end
    ylim([0, 200])
    vline(0, 'k:')
    yL = rangeLine(fig); 
    plot([0, light_on], [yL,yL], 'linewidth', LW+2, 'color', Color(param.LED_light))
    % labels: 
    if rep == num.reps; xlabel('time (s)'); end
    if rep == 1; title(param.conds_matrix(cond+CCW_offset).label); end
    idx = idx+1;   
end


fig_name = strrep(fig_name, 'cw-', '');
save_figure(fig, [fig_dir, '\' FilePath.structure_name ' joint angle timecourse ' fig_name]);

clearvars('-except',initial_vars{:})

%% Overlay all the diff joint angles for a given stim for the flies
fig_root = [fig_dir, '\Condition Figs'];
if ~exist(fig_root, 'dir')
    mkdir(fig_root)
end 
Joints = {'Coxa-Femur', 'Femur-Tibia', 'Tibia-Tarsus'};
% ------ Input 1 ------
leg = 1;
% ---------------------
% expected size: 
num.frames = param.basler_length*fps+1; % total frames expected

% Group the data: 
for cond = 1:num.conds
   for iJ = 1:num.joints    
       JA(iJ,cond).data = []; %reset/designate the new data structure
       for ifly = 1:num.flies
           for rep = 1:num.reps    
               input = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
               JA(iJ,cond).data = [JA(iJ,cond).data, input];
           end
       end
   end   
end



% ---- plot the data for a given joint ---- : 
LW = 1; SZ = 50;
% ------ Input 2 ------
iJ = 2; % Femur-tibia
% ---------------------
controlROI = 1:param.basler_delay*fps; %150 = 0.5sec*fps
stimStart = controlROI(end)+1;
condList = 1:7;

for cond = condList
    
    light_off = param.conds_matrix(cond).opto;
    
    % get regions of interest:
    
    stimROI = stimStart:stimStart+(light_off*fps);
    postROI = stimStart+light_off*fps:param.basler_length*fps; 
    % pull averages from regions of interest
    pre = (nanmean(JA(iJ,cond).data(controlROI,:)));
    preAv = mean(pre);
    stim = (nanmean(JA(iJ,cond).data(stimROI,:)));
    stimAv = mean(stim);
    post = (nanmean(JA(iJ,cond).data(postROI,:)));
    postAv = mean(post);

    %     x = -0.5:1/fps:1.5;
    joint_colors = {'purple', 'darkgoldenrod', 'teal'};

    fig = getfig('',1);
    subplot(1,2,1); hold all
        plot(x,JA(iJ,cond).data, 'Color', Color(joint_colors{iJ}), 'linewidth', LW)
        yL = rangeLine(fig); 
        plot([0, light_off], [yL,yL], 'linewidth', LW+2, 'color', Color(param.LED_light))
        % plot a moving average:
        input = mean(JA(iJ,cond).data,2);
        plot(x, input, 'color', 'k', 'linewidth', LW)
        %labels
        xlabel('time (s)')
        ylabel([Joints{iJ} ' joint angle (deg)'])
        title({FilePath.structure_name; param.conds_matrix(cond).label})

    % scatter plot of joint angle per trial:
    n = size(pre,2);
    subplot(1,2,2); hold all
    % pre
        ROI = [1,2];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        scatter(x1, pre, SZ, 'k', 'filled')
        plot(ROI,[preAv, preAv], 'color', 'k', 'linewidth', LW+2)
    % stim
        ROI = [3,4];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        try
            scatter(x1, stim, SZ, Color(joint_colors{iJ}), 'filled')
            plot(ROI,[stimAv, stimAv], 'color', Color(joint_colors{iJ}), 'linewidth', LW+2)
        catch
        end
    % post
        ROI = [5,6];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        scatter(x1, post, SZ, 'k', 'filled')
        plot(ROI,[postAv, postAv], 'color', 'k', 'linewidth', LW+2)
        xlim([condList(1)-1,condList(end)+1])
        xticklabels([])
        xlabel({'Pre          stim           post'; 'Time period'})
        ylabel('avg joint angle (deg)')
        title(['Leg ' leg_labels{leg}])
    % save the figure:
    save_figure(fig, [fig_root, '\joint angles and averages ' param.conds_matrix(cond).label]);
end

clearvars('-except',initial_vars{:})

%% Plot step frequency and joint angle for a single condition (& all reps) of a single fly
joint_colors = {'purple', 'darkgoldenrod', 'teal'};
LW = 1;
  
% ----- Input -----
leg = 1;
% -----------------

% Group the data for all trials of a condition: 
for cond = 1:num.conds
   for iJ = 1:num.joints    
       JA(iJ,cond).data = [];
       for ifly = 1:num.flies
           for rep = 1:num.reps    
               input = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
               JA(iJ,cond).data = [JA(iJ,cond).data, input];
           end
       end
   end   
end


% ----- Input -----
ifly = 3;   % fly selection
leg = 1;    % leg selection
iJ = 2;     % joint selection
cond = 10;   % condition selection
% -----------------

% find step frequency
% use the function peakfinder to pull the peak locations -- however, it
% doesn't have the capability to find max/min outside of the positive/neg
% range, so through the data into another function first to adjust it to
% the peakfinder standards


% FilePath.locations{ifly,4}

fig = getfig('',1); 
for rep = 1:num.reps
        [nrows, ncols] = subplot_numbers(num.reps);
        subplot(nrows, ncols, rep)
        hold all
    % plot angle timecourse    
%         input = JA(iJ,cond).data(:,rep); % use data from JA grouping (not
%         easily organized by fly)
        input = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
        plot(input, 'color', Color(joint_colors{iJ}), 'linewidth', LW)
    % find max and min points in the joint angles
        max_loc = findPeaks('-x',(input), '-extrema', 1);
        max_val = input(round(max_loc));
        min_loc = findPeaks('-x',(input), '-extrema', -1);
        min_val = input(round(min_loc));

    % check that the peaks qualify for the minimum oscillation amplitude:
    min_diff = 10;
    oscillation_diff = max(max_val) - min(min_val);
    if oscillation_diff <= min_diff
        [max_loc, max_val, min_loc, min_val] = deal([]);
        fprintf('\n No peaks found\n')
    end

    % plot flexion points
        scatter(max_loc, max_val, 50, 'r', 'filled');
        scatter(min_loc, min_val, 50, 'b', 'filled');
        ylabel([Joints{iJ} ' joint angle (deg)'])
        xlabel('time (frames)')
        light_on = 150+param.conds_matrix(cond).opto*fps;
        y = rangeLine(fig);
        plot([150,light_on], [y,y], 'color', 'g', 'linestyle', '-', 'linewidth', LW+2)
        ylim([0,180])
    % rolling step frequency:
        yyaxis right
        stepfrequency = 1./(diff(max_loc)/fps);
        plot(max_loc(2:end),stepfrequency, 'color', 'k', 'linewidth', LW)
        ylabel('step frequency (steps/sec)')
    % labels and general
        title({FilePath.structure_name;...
            ['Condition: ' param.conds_matrix(cond).label]; ['Rep: ' num2str(rep)]})

end

fig_name = [FilePath.structure_name ' ' param.conds_matrix(cond).label ' peak selection'];
save_figure(fig, [fig_dir, '\' fig_name]);

clearvars('-except',initial_vars{:})

%% ----- SECTIONS IN PROGRESS BELOW HERE ----- 
%aka variables aren't clearable and need to be processed in order

% % plot the angle data for any trial:
% % ----- Input ------ 
% ifly = 14;  % fly number in structure
% cond = 1;   % condition
% rep = 2;    % repitition
% leg = 1;    % leg on fly
% iJ = 2;     % joint
% % ------------------
% 
% figure;
% plot(angles(ifly).Leg(leg).data(cond, rep).all(:,iJ))

% --------------------------------------------

%% Behavior predictions from the FeTi joint angle of a single leg
% DATA ORGANIZATION SECTION the 'step' information into a structure
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

% Pull the information for each fly:
for ifly = 1:num.flies
  % set blank behavior group:
  [behavior(ifly).walking, behavior(ifly).stationary, behavior(ifly).other] = ...
      deal(false(num.conds,num.reps));
  for leg = 1:6 %number of legs
     for cond = 1:num.conds
        for rep = 1:num.reps
            raw = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
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
            
            % Categorize behavior:
            %walking
            if freq >= min_freq && step_std <= max_std
                angles(ifly).Leg(leg).control_behavior{cond,rep} = 'walking';
                if leg == behavior_leg
                    behavior(ifly).walking(cond,rep) = true;
                end
            end
            %stationary  
            cntl_angles = raw(controlROI);
            a = abs(mean(diff(cntl_angles))) < max_variance; 
            b = (max(cntl_angles)- min(cntl_angles)) < max_range;
            if a && b == true
                angles(ifly).Leg(leg).control_behavior{cond,rep} = 'stationary';
                if leg == behavior_leg
                    behavior(ifly).stationary(cond,rep) = true;
                end
            end
        end
     end    
  end
end

% find 'other' behavior and summarize
for ifly = 1:num.flies
   behavior(ifly).other =...
       ~behavior(ifly).walking & ~behavior(ifly).stationary;
   % Summarize the behavior into name classes
   behavior(ifly).class = nan(num.conds,num.reps);
   loc = behavior(ifly).walking;
   behavior(ifly).class(loc) = 1;
   loc = behavior(ifly).stationary;
   behavior(ifly).class(loc) = 2;
   loc = behavior(ifly).other;
   behavior(ifly).class(loc) = 3;
end

clearvars('-except',initial_vars{:})

%% Plot FeTi joint angle for all trials of a given condition (all flies)
% color coded behavior sorting too
LW = 1;

% ----- Input ------ 
cond = 26;   % condition
iJ = 2;     % Femur-Tibia joint number
leg = 1;    % leg joint angle trace to plot
% ------------------

light_on = (param.basler_delay*fps);
light_off = param.conds_matrix(cond).opto*fps;
light_off = light_off+light_on;

[nrows, ncols] = subplot_numbers(num.flies*num.reps);

fig_title = [FilePath.structure_name, ' cond ' num2str(cond) ' behavior'];
idx = 0;
fig = getfig(fig_title,1);
for ifly = 1:num.flies
  for rep = 1:num.reps
      idx = idx+1;
      subplot(nrows,ncols,idx)
      % pull the data:
      raw = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
      switch behavior(ifly).class(cond,rep)
          case 1 %walking
            pcolor = Color('darkcyan');
            tag = 'W'; %'walking';
          case 2 %stationary
            pcolor = Color('firebrick');
            tag = 'S'; %'stationary';
          case 3 %other
            pcolor = Color('grey');
            tag = 'O'; %'other';
      end

      plot(raw, 'linewidth', LW, 'color', pcolor) %controlROI
      xticklabels([])
      yticklabels([])  
      ylim([0,200])
      xlim([0,length(raw)])
      vline([light_on,light_off], 'k-')
      % write the ID info:
      fly_ID = [FilePath.locations{ifly,4} ' R' num2str(rep) 'C' num2str(cond)];
      fly_ID = strrep(fly_ID, '_', '-');
      xlabel(fly_ID)
      ylabel(tag)
  end
end
%labels:
subplot(nrows,ncols,1)
title(['All trials cond ' num2str(cond)])
    
save_figure(fig, [fig_dir '\' fig_title]);

clearvars('-except',initial_vars{:})

%% Plot all trials for a selected fly with behavior labeled
% color coded behavior sorting too
LW = 1;

% ----- Input ------ 
ifly = 3;   % fly
iJ = 2;      % joint
leg = 1;     % leg choice
% ------------------

B_delay = param.basler_delay;
light_on = (fps*B_delay);
%subplot sizes
ncols = length(param.LED);
nrows = round(num.conds/ncols);
edge_num = 1:ncols:num.conds; %e.g. 1,8,15,22
bottom_num = num.conds-ncols+1:num.conds; %e.g. 22:28

fig_title = [FilePath.structure_name, ' ' FilePath.locations{ifly,4} ' all trials'];

fig = getfig(fig_title,1);
for cond = 1:num.conds
    subplot(nrows,ncols,cond); hold all
    for rep = 1:num.reps
        raw = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
        switch behavior(ifly).class(cond,rep)
            case 1 % walking
                pcolor = Color('darkcyan');
                tag{rep} = 'W'; %'walking';
            case 2 % stationary
                pcolor = Color('firebrick');
                tag{rep} = 'S'; %'stationary';
            case 3 % other
                pcolor = Color('grey');
                tag{rep} = 'O'; %'other';
        end
        plot(raw, 'linewidth', LW, 'color', pcolor) %controlROI
    end    
    % labels
      xticklabels([])
      yticklabels([])  
      ylim([0,180])
      xlim([0,length(raw)])
      light_off = param.conds_matrix(cond).opto*fps+light_on;
      vline(light_on, 'k:')
      y1 = rangeLine(fig);
      plot([light_on, light_off], [y1,y1], 'linewidth', LW+2, 'color', 'g')

      % conditional labels:
      t_add = [tag{1}, ', ', tag{2}, ', ', tag{3}];
      if cond == 1
          title({strrep(FilePath.locations{ifly,4}, '_', '-'); t_add})
      else 
          title(t_add)
      end
      % ylabels
      if sum(cond == edge_num) > 0 
          ylabel([Joints{iJ} ' joint angle'])
      end
      
      % xlabels
      if sum(cond == bottom_num) > 0 
          xlabel('time (s)')
      end
                                                                
end

save_figure(fig, [fig_dir '\' fig_title]);

clearvars('-except',initial_vars{:})

%% compare determined behaviors with those manually labeled:
clear temp
load(['C:\matlabroot\behavior class\' FilePath.cross ' behavior class']);

% Add behavior class fields to the structure:
for ifly = 1:num.flies
    input = group(ifly).behavior;
    group(ifly).walking = strcmpi(input, 'Walking');
    group(ifly).stationary = strcmpi(input, 'Stationary');
    group(ifly).other = strcmpi(input, 'Other');
    group(ifly).grooming = strcmpi(input, 'Grooming');
end

% ------ input ------
iJ = 2; % Femur-tibia
% -------------------

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

%% Overlay of joint angles for all trials for a specific behavior (e.g. walking)
fig_root = [fig_dir, '\Condition Figs'];
if ~exist(fig_root, 'dir')
    mkdir(fig_root)
end 
% ------ Input 1 ------
leg = 1;
state = 'walking';
% ---------------------

% Group the data: 
for cond = 1:num.conds
   for iJ = 1:num.joints    
       JA(iJ,cond).data = []; %reset/designate the new data structure
       for ifly = 1:num.flies
           for rep = 1:num.reps 
               % filter for desired behavior
               if behavior(ifly).(state)(cond,rep) == true
                  input = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
                  JA(iJ,cond).data = [JA(iJ,cond).data, input];
               end
           end
       end
   end   
end

% ---- plot the data for a given joint ---- : 
LW = 1; SZ = 50;
% ------ Input 2 ------
spLoc = 1:2:num.joints*2;
% ---------------------
controlROI = 1:param.basler_delay*fps; %150 = 0.5sec*fps
stimStart = controlROI(end)+1;
condList = 2;

for cond = condList
  fig = getfig('',1);
  
  for iJ = 1:num.joints  
    subplot(num.joints,2,spLoc(iJ)); hold all  
    light_off = param.conds_matrix(cond).opto;
    
    % get regions of interest:
    stimROI = stimStart:stimStart+(light_off*fps);
    postROI = stimStart+light_off*fps:param.basler_length*fps; 
    % pull averages from regions of interest
    pre = (nanmean(JA(iJ,cond).data(controlROI,:)));
    preAv = mean(pre);
    stim = (nanmean(JA(iJ,cond).data(stimROI,:)));
    stimAv = mean(stim);
    post = (nanmean(JA(iJ,cond).data(postROI,:)));
    postAv = mean(post);

    %     x = -0.5:1/fps:1.5;
    joint_colors = {'purple', 'darkgoldenrod', 'teal'};

    
        plot(x,JA(iJ,cond).data, 'Color', Color(joint_colors{iJ}), 'linewidth', LW)
        yL = rangeLine(fig); 
        plot([0, light_off], [yL,yL], 'linewidth', LW+2, 'color', Color(param.LED_light))
        % plot a moving average:
        input = mean(JA(iJ,cond).data,2);
        plot(x, input, 'color', 'k', 'linewidth', LW)
        %labels
        xlabel('time (s)')
        ylabel([Joints{iJ} ' joint angle (\circ)'])
        title({FilePath.structure_name; param.conds_matrix(cond).label})

    % scatter plot of joint angle per trial:
    n = size(pre,2);
    subplot(num.joints,2,spLoc(iJ)+1); hold all
    % pre
        ROI = [1,2];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        scatter(x1, pre, SZ, 'k', 'filled')
        plot(ROI,[preAv, preAv], 'color', 'k', 'linewidth', LW+2)
    % stim
        ROI = [3,4];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        try
            scatter(x1, stim, SZ, Color(joint_colors{iJ}), 'filled')
            plot(ROI,[stimAv, stimAv], 'color', Color(joint_colors{iJ}), 'linewidth', LW+2)
        catch
        end
    % post
        ROI = [5,6];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        scatter(x1, post, SZ, 'k', 'filled')
        plot(ROI,[postAv, postAv], 'color', 'k', 'linewidth', LW+2)
        xlim([0,7]) %
        xticklabels([])
        xlabel({'Pre          stim           post'; 'Time period'})
        ylabel('avg joint angle (\circ)')
        title(['Leg ' leg_labels{leg}])
  end
        
    % save the figure:
    save_figure(fig, [fig_root, '\joint angles and averages ' param.conds_matrix(cond).label]);
end

clearvars('-except',initial_vars{:})

%% Overlay of joint angles for all trials for a specific behavior (e.g. stationary)
fig_root = [fig_dir, '\Condition Figs'];
if ~exist(fig_root, 'dir')
    mkdir(fig_root)
end 
% ------ Input 1 ------
leg = 1;
% state = 'stationary';
state = 'walking';
% ---------------------

% Group the data: 
for cond = 1:num.conds
   for iJ = 1:num.joints    
       JA(iJ,cond).data = []; %reset/designate the new data structure
       for ifly = 1:num.flies
           for rep = 1:num.reps 
               % filter for desired behavior
               if behavior(ifly).(state)(cond,rep) == true
                  input = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
                  JA(iJ,cond).data = [JA(iJ,cond).data, input];
               end
           end
       end
   end   
end

% ---- plot the data for a given joint ---- : 
LW = 1; SZ = 50;
% ------ Input 2 ------
nc = 3;
spLoc = 1:nc:num.joints*nc;
% ---------------------
controlROI = 1:param.basler_delay*fps; %150 = 0.5sec*fps
stimStart = controlROI(end)+1;
condList = 1:7;

for cond = condList
  fig = getfig('',1);
  
  for iJ = 1:num.joints  
%Raw data:
    subplot(num.joints,nc,spLoc(iJ)); hold all  
    light_off = param.conds_matrix(cond).opto;
    
    % get regions of interest:
    stimROI = stimStart:stimStart+(light_off*fps);
    postROI = stimStart+light_off*fps:param.basler_length*fps; 
    % pull averages from regions of interest
    raw = JA(iJ,cond).data;
    pre = (nanmean(raw(controlROI,:)));
    preAv = mean(pre);
    stim = (nanmean(raw(stimROI,:)));
    stimAv = mean(stim);
    post = (nanmean(raw(postROI,:)));
    postAv = mean(post);

    %     x = -0.5:1/fps:1.5;
    joint_colors = {'purple', 'darkgoldenrod', 'teal'};

    
        plot(x,JA(iJ,cond).data, 'Color', Color(joint_colors{iJ}), 'linewidth', LW)
        yL = rangeLine(fig); 
        plot([0, light_off], [yL,yL], 'linewidth', LW+2, 'color', Color(param.LED_light))
        % plot a moving average:
        input = mean(JA(iJ,cond).data,2);
        plot(x, input, 'color', 'k', 'linewidth', LW)
        %labels
        xlabel('time (s)')
        ylabel([Joints{iJ} ' joint angle (\circ)'])
        title({FilePath.structure_name; param.conds_matrix(cond).label})

    
%AVG JOINT ANGLE
    % scatter plot of joint angle per trial:
    n = size(pre,2);
    subplot(num.joints,nc,spLoc(iJ)+1); hold all
    % pre
        ROI = [1,2];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        scatter(x1, pre, SZ, 'k', 'filled')
        plot(ROI,[preAv, preAv], 'color', 'k', 'linewidth', LW+2)
    % stim
        ROI = [3,4];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        try
            scatter(x1, stim, SZ, Color(joint_colors{iJ}), 'filled')
            plot(ROI,[stimAv, stimAv], 'color', Color(joint_colors{iJ}), 'linewidth', LW+2)
        catch
        end
    % post
        ROI = [5,6];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        scatter(x1, post, SZ, 'k', 'filled')
        plot(ROI,[postAv, postAv], 'color', 'k', 'linewidth', LW+2)
        xlim([0,7]) %
        xticklabels([])
        xlabel({'Pre          stim           post'; 'Time period'})
        ylabel('avg joint angle (\circ)')
        title(['Leg ' leg_labels{leg}])
        
        
%ACTIVITY LEVEL
    % scatter plot of total movement as proxy for activity:
    temp = raw(controlROI,:);
    pre = sum(abs(diff(temp)))/size(temp,1); % avg change in joint angle 
    preAv = nanmean(pre);
    temp = raw(stimROI,:);
    stim = sum(abs(diff(temp)))/size(temp,1); % avg change in joint angle
    stimAv = nanmean(stim);
    temp = raw(postROI,:);
    post = sum(abs(diff(temp)))/size(temp,1); % avg change in joint angle 
    postAv = nanmean(post);
    
    n = size(pre,2); % number of trials
    subplot(num.joints,nc,spLoc(iJ)+2); hold all
    % pre
        ROI = [1,2];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        scatter(x1, pre, SZ, 'k', 'filled')
        plot(ROI,[preAv, preAv], 'color', 'k', 'linewidth', LW+2)
    % stim
        ROI = [3,4];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        try
            scatter(x1, stim, SZ, Color(joint_colors{iJ}), 'filled')
            plot(ROI,[stimAv, stimAv], 'color', Color(joint_colors{iJ}), 'linewidth', LW+2)
        catch
        end
    % post
        ROI = [5,6];
        x1 = scatterPoints(ROI, n); %minval, maxval, number of samples
        scatter(x1, post, SZ, 'k', 'filled')
        plot(ROI,[postAv, postAv], 'color', 'k', 'linewidth', LW+2)
        xlim([0,7]) %
        xticklabels([])
        xlabel({'Pre          stim           post'; 'Time period'})
        ylabel('avg \Delta joint angle (\circ)')
        title(['Leg ' leg_labels{leg}])
  end
        
    % save the figure:
    save_figure(fig, [fig_root, '\' state ' joint angles and averages ' param.conds_matrix(cond).label]);
end

clearvars('-except',initial_vars{:})

 %% Stationary Flies: SINGLE condition overlay (all flies)
% Overlay the abs and change in joint angle over time in stationary flies:
LW = 1;
% ----- input ------
cond = 2;
leg = 1;
iJ = 2;
% ------------------
fig = getfig('',1); hold all
for ifly = 1:num.flies
    for rep = 1:num.reps
        % screen for desired behavior
        if behavior(ifly).stationary(cond,rep)==true
            raw = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
        else
            continue; %not stationary, skip to next 'for' cycle
        end
        % plot the data:
        ylim([0,180])
        plot(x, raw, 'color', 'k', 'linewidth', LW)
        % plot the light region:
        y1 = rangeLine(fig);
        plot([0,param.conds_matrix(cond).opto], [y1,y1], 'linewidth', LW+2, 'color', 'g')
    end
end
% labels
vline([0, param.conds_matrix(cond).opto],'k:')
title(strrep(param.conds_matrix(cond).label, '_', '-'))
xlabel('time (s)')
ylabel([Joints{iJ} ' joint angle'])

fig_title = [FilePath.structure_name ' ' param.conds_matrix(cond).label ' stationary all trials'];
save_figure(fig, [fig_dir '\' fig_title]);

clearvars('-except',initial_vars{:})

%% Stationary Flies: ALL conditions (subplots) - avg of trials
% Overlay the abs and change in joint angle over time in stationary flies:
LW = 1;
ncols = length(param.LED);
nrows = num.conds/ncols;

% ----- input ------
leg = 1;
iJ = 2;
% ------------------

fig = getfig('',1); 
for cond = 1:num.conds
    subplot(nrows,ncols,cond)
    hold all
    data = [];
    for ifly = 1:num.flies
        for rep = 1:num.reps
            % screen for desired behavior
            if behavior(ifly).stationary(cond,rep)==true
                raw = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
            else
                continue; %not stationary, skip to next 'for' cycle
            end
            data = [data,raw];
        end
    end
    % plot the data:
%     ylim([0,180])
    xlim([x(1), x(end)])
    % normalize the data to change in joint angle:
    data = data-data(150,:); %normalize by the last control frame
    y.avg = median(data,2);
    y.err = std(data,0,2);
    fill_data = error_fill(x, y.avg, y.err);
    h = fill(fill_data.X, fill_data.Y, 'k', 'EdgeColor','none');
    set(h, 'facealpha', 0.2)
    plot(x, y.avg, 'color', 'k', 'linewidth', LW)

    % plot the light region:
    vline(0,'k:')
    y1 = rangeLine(fig);
    plot([0,param.conds_matrix(cond).opto], [y1,y1], 'linewidth', LW+2, 'color', 'g')
    
    % labels
    title(strrep(param.conds_matrix(cond).label, '_', '-'))
    xlabel('time (s)')
    ylabel([Joints{iJ} ' joint angle'])
end

fig_title = [FilePath.structure_name ' avg change in ' Joints{iJ} ' angle STATIONARY'];
save_figure(fig, [fig_dir '\' fig_title], '-png');

clearvars('-except',initial_vars{:})

%% Overlay one step cycle for all WALKING flies in the control period
% step 1) Isolate the walking trials for a given fly
% step 2) Find the start of each swing cycle (local minimma)
% step 3) Find the subsequent start of a swing cycle
% step 4) Isolate regions of joint angle data that fit
% step 5) Plot and overlay the trajectories
% step 5b)Optional color coordination by speed

% ----- input ------
ifly = 5;   %fly number
leg = 1;    %leg
iJ = 2;     %joint
% ------------------

ROI = 1:(param.basler_delay*fps);
LW = 1;
% build empty matrix for the traces
data = nan(length(ROI),1);
MT = data;

fig = getfig('',1); hold all
for cond = 1:num.conds
   for rep = 1:num.reps
      if behavior(ifly).walking(cond,rep)==false
          continue; % skip other behaviors
      end
      raw = angles(ifly).Leg(leg).data(cond,rep).all(ROI,iJ); % joint angle data
      % find start of swing
      min_loc = angles(ifly).Leg(leg).step(cond,rep).min_loc;
      min_loc(min_loc>ROI(end)) = [];
      min_loc = round(min_loc);
      % find start of stance
      max_loc = angles(ifly).Leg(leg).step(cond,rep).max_loc;
      max_loc(max_loc>ROI(end)) = [];
      max_loc = round(max_loc);
      % isolate steps 
      nSteps = length(min_loc)-1; % number of steps
      if nSteps==0
          continue; % skip trials with no complete steps (should be screened by walking already tho)
      end
      for n = 1:nSteps
        % check for max values between swing points
        step_range = min_loc(n):min_loc(n+1);
%         peak = sum(max_loc>min_loc(n) & max_loc<min_loc(n+1));
%         if peak == 0 
%             continue; %if no peak between steps, skip it
%         end
        step = raw(step_range);
        %save data to matrix
        input = MT;
        input(1:length(step)) = step;
        data = [data, input];
        %plot the step
        plot(step, 'color', Color('grey'), 'linewidth', LW)
      end
   end
end; input = [];
data(:,1) = []; %remove empty first column
% Overlay the avg trace:
input = nanmedian(data,2);
input(isnan(input)) = [];
% only plot to point that 25% of steps are done:
loc = sum(isnan(data),2)>round(size(data,2)*0.25); 
input(loc(1:length(input))) = [];
plot(input, 'color', Color('darkcyan'), 'linewidth', LW+4)
xlim([0,length(input)])
% plot labels
fig_name = strrep([FilePath.locations{ifly,4} ' control step'], '_', '-');
title(fig_name)
xlabel('time (fps)')
ylabel([Joints{iJ} ' angle (deg)'])


fig_title = [FilePath.structure_name ' ' fig_name];
save_figure(fig, [fig_dir '\' fig_title]);

clearvars('-except',initial_vars{:})

%% Avg one step cycle for all WALKING flies in the control period of ALL flies single plot
% step 1) Isolate the walking trials for a given fly
% step 2) Find the start of each swing cycle (local minimma)
% step 3) Find the subsequent start of a swing cycle
% step 4) Isolate regions of joint angle data that fit
% step 5) Plot and overlay the trajectories
% step 5b)Optional color coordination by speed

% ----- input ------
leg = 1;    %leg
iJ = 2;     %joint
% ------------------

ROI = 1:(param.basler_delay*fps);
LW = 1;
% build empty matrix for the traces
data = nan(length(ROI),1);
MT = data;
[nrows,ncols] = subplot_numbers(num.flies);

fig = getfig('',1);
for ifly = 1:num.flies
    data = nan(length(ROI),1);
    MT = data;
    subplot(nrows,ncols,ifly); hold all
    for cond = 1:num.conds
       for rep = 1:num.reps
          if behavior(ifly).walking(cond,rep)==false
              continue; % skip other behaviors
          end
          raw = angles(ifly).Leg(leg).data(cond,rep).all(ROI,iJ); % joint angle data
          % find start of swing
          min_loc = angles(ifly).Leg(leg).step(cond,rep).min_loc;
          min_loc(min_loc>ROI(end)) = [];
          min_loc = round(min_loc);
          % find start of stance
          max_loc = angles(ifly).Leg(leg).step(cond,rep).max_loc;
          max_loc(max_loc>ROI(end)) = [];
          max_loc = round(max_loc);
          % isolate steps 
          nSteps = length(min_loc)-1; % number of steps
          if nSteps==0
              continue; % skip trials with no complete steps (should be screened by walking already tho)
          end
          for n = 1:nSteps
            % check for max values between swing points
            step_range = min_loc(n):min_loc(n+1);
    %         peak = sum(max_loc>min_loc(n) & max_loc<min_loc(n+1));
    %         if peak == 0 
    %             continue; %if no peak between steps, skip it
    %         end
            step = raw(step_range);
            %save data to matrix
            input = MT;
            input(1:length(step)) = step;
            data = [data, input];
            %plot the step (optional)
            plot(step, 'color', Color('grey'), 'linewidth', LW)
          end
       end
    end; input = [];
    data(:,1) = []; %remove empty first column
    % Overlay the avg trace:
    input = nanmean(data,2);
    input(isnan(input)) = [];
    % only plot to point that 25% of steps are done:
    loc = sum(isnan(data),2)>round(size(data,2)*0.25); 
    input(loc(1:length(input))) = [];
    plot(input, 'color', Color('darkcyan'), 'linewidth', LW+2)
    xlim([0,length(input)])
    % plot labels
    fig_name = strrep([FilePath.locations{ifly,4} ' control step'], '_', '-');
    title(fig_name)
    xlabel('time (fps)')
    ylabel([Joints{iJ} ' angle (deg)'])
    clear data input loc
end

fig_title = [FilePath.structure_name ' control steps overlaid all trials'];
save_figure(fig, [fig_dir '\' fig_title]);

clearvars('-except',initial_vars{:})

%% Avg WALKING step cycle for all flies in the control period

% ----- input ------
leg = 1;    %leg
iJ = 2;     %joint
% ------------------

ROI = 1:(param.basler_delay*fps);
LW = 1;

fig = getfig('',1); hold all
for ifly = 1:num.flies
    data = nan(length(ROI),1);
    MT = data;
    for cond = 1:num.conds
       for rep = 1:num.reps
          if behavior(ifly).walking(cond,rep)==false
              continue; % skip other behaviors
          end
          raw = angles(ifly).Leg(leg).data(cond,rep).all(ROI,iJ); % joint angle data
          % find start of swing
          min_loc = angles(ifly).Leg(leg).step(cond,rep).min_loc;
          min_loc(min_loc>ROI(end)) = [];
          min_loc = round(min_loc);
          % find start of stance
          max_loc = angles(ifly).Leg(leg).step(cond,rep).max_loc;
          max_loc(max_loc>ROI(end)) = [];
          max_loc = round(max_loc);
          % isolate steps 
          nSteps = length(min_loc)-1; % number of steps
          if nSteps==0
              continue; % skip trials with no complete steps (should be screened by walking already tho)
          end
          for n = 1:nSteps
            % check for max values between swing points
            step_range = min_loc(n):min_loc(n+1);
            step = raw(step_range);
            step = step-(step(1)); % normalize data
            %save data to matrix
            input = MT;
            input(1:length(step)) = step;
            data = [data, input];
          end
       end
    end; input = [];
    data(:,1) = []; %remove empty first column
    % Overlay the avg trace:
    input = nanmedian(data,2);
    input(isnan(input)) = [];
    % only plot to point that 25% of steps are done:
    loc = sum(isnan(data),2)>round(size(data,2)*0.25); 
    input(loc(1:length(input))) = [];
    plot(input, 'linewidth', LW+2)

    clear data input loc
end

% plot labels
fig_name = [FilePath.structure_name ' control step overlaid all flies'];
title(fig_name)
xlabel('time (fps)')
ylabel([Joints{iJ} ' change in angle (deg)'])
tester = FilePath.locations;
A = cellfun(@(tester) strrep(tester, '_','-'), tester, 'UniformOutput', false); %replace '_'
legend(A{:,4})

save_figure(fig, [fig_dir '\' fig_name]);
clearvars('-except',initial_vars{:})

%% Plot step frequency and joint angle for each condition (& all reps) for a behavior
joint_colors = {'purple', 'darkgoldenrod', 'teal'};
LW = 1;
  
% ----- Input -----
leg = 1;
% -----------------

% Group the data for all trials of a condition: 
for cond = 1:num.conds
   for iJ = 1:num.joints    
       JA(iJ,cond).data = [];
       for ifly = 1:num.flies
           for rep = 1:num.reps    
               input = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
               JA(iJ,cond).data = [JA(iJ,cond).data, input];
           end
       end
   end   
end


% ----- Input -----
ifly = 5;   % fly selection
leg = 1;    % leg selection
iJ = 2;     % joint selection
cond = 4;   % condition selection
% -----------------

% find step frequency
% use the function peakfinder to pull the peak locations -- however, it
% doesn't have the capability to find max/min outside of the positive/neg
% range, so through the data into another function first to adjust it to
% the peakfinder standards


% FilePath.locations{ifly,4}

fig = getfig('',1); 
for rep = 1:num.reps
        [nrows, ncols] = subplot_numbers(num.reps);
        subplot(nrows, ncols, rep)
        hold all
    % plot angle timecourse    
%         input = JA(iJ,cond).data(:,rep); % use data from JA grouping (not
%         easily organized by fly)
        input = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ);
        plot(input, 'color', Color(joint_colors{iJ}), 'linewidth', LW)
    % find max and min points in the joint angles
        max_loc = findPeaks('-x',(input), '-extrema', 1);
        max_val = input(round(max_loc));
        min_loc = findPeaks('-x',(input), '-extrema', -1);
        min_val = input(round(min_loc));

    % check that the peaks qualify for the minimum oscillation amplitude:
    min_diff = 10;
    oscillation_diff = max(max_val) - min(min_val);
    if oscillation_diff <= min_diff
        [max_loc, max_val, min_loc, min_val] = deal([]);
        fprintf('\n No peaks found\n')
    end

    % plot flexion points
        scatter(max_loc, max_val, 50, 'r', 'filled');
        scatter(min_loc, min_val, 50, 'b', 'filled');
        ylabel([Joints{iJ} ' joint angle (deg)'])
        xlabel('time (frames)')
        light_on = 150+param.conds_matrix(cond).opto*fps;
        y = rangeLine(fig);
        plot([150,light_on], [y,y], 'color', 'g', 'linestyle', '-', 'linewidth', LW+2)
        ylim([0,180])
    % rolling step frequency:
        yyaxis right
        stepfrequency = 1./(diff(max_loc)/fps);
        plot(max_loc(2:end),stepfrequency, 'color', 'k', 'linewidth', LW)
        ylabel('step frequency (steps/sec)')
    % labels and general
        title({FilePath.structure_name;...
            ['Condition: ' param.conds_matrix(cond).label]; ['Rep: ' num2str(rep)]})

end

fig_name = [FilePath.structure_name ' ' param.conds_matrix(cond).label ' peak selection'];
save_figure(fig, [fig_dir, '\' fig_name]);

clearvars('-except',initial_vars{:})

%% vvvv WORKING SECTION HERE vvvv





% Position data below: 

%% Overlay swing/stance in walking trials: % WORKNG HERE! START ON WED AGAIN
% step 1) Isolate the walking trials only for a given condition:
% step 2) Find the avg control step
% step 3) Determine if a fly is in swing or stance during opto onset
% step 4) Align 1st stim-step for 'swing' and for 'stance'


% ----- input ------
leg = 1;    %leg
iJ = 2;     %joint
cond = 5;
% ------------------
nrows = 2; % one line
ncols = 3; % control|swing|stance|overlaid|outliers
% time limit post stim on for there to be a step
light_on = (param.basler_delay*fps);
t_lim = (0.5*fps)+light_on; 
t_frames = light_on:light_on+60;

ROI = 1:(param.basler_delay*fps);
LW = 1;
[control, swing, stance] = deal(nan(length(ROI),1));
MT = control;

fig = getfig('',1); 
for ifly = 1:num.flies    
   for rep = 1:num.reps
      if behavior(ifly).walking(cond,rep)==false
          continue; % skip other behaviors
      end
      raw = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ); % joint angle data
      % find start of swing
      min_loc = round(angles(ifly).Leg(leg).step(cond,rep).min_loc);
      % find start of stance
      max_loc = round(angles(ifly).Leg(leg).step(cond,rep).max_loc);
      
      % --- control ---- %
      subplot(nrows,ncols,1); hold all
      swing_strt = round(min_loc(min_loc<ROI(end)));
      stance_strt = round(max_loc(max_loc<ROI(end)));
      nSteps = length(swing_strt)-1; % number of steps
      for n = 1:nSteps
       % check for max values between swing points
        step_range = swing_strt(n):swing_strt(n+1);
        step = raw(step_range);
        step = step-(step(1)); % normalize data
        %save data to matrix
        input = MT;
        input(1:length(step)) = step;
        control = [control, input];
        % plot the trials:
        plot(step, 'color', Color('grey'), 'linewidth', LW)
      end
      
      % determine if fly is in swing or stance at stim start:
      type = swing_strt(end)>stance_strt(end);
      last_swing = find(min_loc==swing_strt(end)); % idx of last swing start in control
      
      % check for flies stopping walking completely...
      post_stim_steps = sum(min_loc>light_on & min_loc<t_lim)-1;
      if post_stim_steps < 2 % min num steps = 2
           switch type
               case true
               % ------ SWING ------ %
                   subplot(nrows,ncols,5); hold all
                   step = raw(t_frames);
                   step = step-(step(1)); % normalize data
                   plot(step, 'color', Color('mediumpurple'), 'linewidth', LW)
               case false
               % ------ SWING ------ %
                   subplot(nrows,ncols,6); hold all
                   step = raw(t_frames);
                   step = step-(step(1)); % normalize data
                   plot(step, 'color', Color('Turquoise'), 'linewidth', LW)
           end
           continue; % jump to next trial
      end

      switch type
          % ------ SWING ------ %
          case true
              subplot(nrows,ncols,2); hold all
              s_start = min_loc(last_swing);
              s_end = min_loc(last_swing+1);
              step_range = s_start:s_end;
              step = raw(step_range);
              step = step-(step(1)); % normalize data
              plot(step, 'color', Color('mediumpurple'), 'linewidth', LW)
              %save data to matrix
              input = MT;
              input(1:length(step)) = step;
              swing = [swing, input];
              
          % ------ STANCE ------ %    
          case false
              subplot(nrows,ncols,3); hold all
              s_start = min_loc(last_swing+1);
              s_end = min_loc(last_swing+2);
              step_range = s_start:s_end;
              step = raw(step_range);
              step = step-(step(1)); % normalize data
              plot(step, 'color', Color('Turquoise'), 'linewidth', LW)
              %save data to matrix
              input = MT;
              input(1:length(step)) = step;
              stance = [stance, input];
      end
   end
end; input = [];

% Plot the avgs and then the overlays:
% CONTROL
subplot(nrows,ncols,1); hold all
control(:,1) = [];
input = nanmedian(control,2);
input(isnan(input)) = [];
% only plot to point that 25% of steps are done:
if size(control,2)> 1
    loc = sum(isnan(control),2)>round(size(control,2)*0.25); 
    input(loc(1:length(input))) = [];
end
plot(input, 'linewidth', LW+2, 'color', 'k')
subplot(nrows,ncols,4); hold all
plot(input, 'linewidth', LW+2, 'color', 'k')

% SWING:
subplot(nrows,ncols,2); hold all
swing(:,1) = [];
input = nanmedian(swing,2);
input(isnan(input)) = [];
% only plot to point that 25% of steps are done:
if size(swing,2)> 1
    loc = sum(isnan(swing),2)>round(size(swing,2)*0.25); 
    input(loc(1:length(input))) = [];
end
plot(input, 'linewidth', LW+2, 'color', Color('Purple'))
subplot(nrows,ncols,4); hold all
plot(input, 'linewidth', LW+2, 'color', Color('Purple'))


% STANCE:
subplot(nrows,ncols,3); hold all
stance(:,1) = [];
input = nanmedian(stance,2);
input(isnan(input)) = [];
% only plot to point that 25% of steps are done:
if size(stance,2)> 1
    loc = sum(isnan(stance),2)>round(size(stance,2)*0.25); 
    input(loc(1:length(input))) = [];
end
plot(input, 'linewidth', LW+2, 'color', Color('Teal'))
subplot(nrows,ncols,4); hold all
plot(input, 'linewidth', LW+2, 'color', Color('Teal'))

% ----- plot labels -----:
%control
subplot(nrows,ncols,1)
title({param.conds_matrix(cond).label;'Control Step'})
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])

%swing
subplot(nrows,ncols,2)
title('First step (swing)')
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])

%stance
subplot(nrows,ncols,3)
title('First step (stance)')
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])

%overlay
subplot(nrows,ncols,4)
title('Avg steps overlaid')
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])

%swing
subplot(nrows,ncols,5)
title('Outliers (swing)')
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])

%stance
subplot(nrows,ncols,6)
title('Outliers (stance)')
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])


fig_title = [FilePath.structure_name ' ' param.conds_matrix(cond).label ' first stim step'];
save_figure(fig, [fig_dir '\' fig_title]);

clearvars('-except',initial_vars{:})
% 
%% Overlay 1st step from swing|stance transitions in walking trials: % 
% pull the swing-stance AND the stance-swing transitions (all flies)

% ----- input ------
leg = 1;    %leg
iJ = 2;     %joint
% ------------------
nrows = 2; % top row = swing-stance, bottom = stance-swing
ncols = 3; % control|swing|stance|overlaid|outliers
% time limit post stim on for there to be a step
light_on = (param.basler_delay*fps);
t_lim = (0.5*fps)+light_on; 
t_frames = light_on:light_on+60;

ROI = 1:(param.basler_delay*fps);
LW = 1;
[control_swing, control_stance, swing, stance] = deal(nan(length(ROI),1));
MT = control_swing;

fig = getfig('',1); 
for ifly = 1:num.flies    
   for rep = 1:num.reps
      if behavior(ifly).walking(cond,rep)==false
          continue; % skip other behaviors
      end
      raw = angles(ifly).Leg(leg).data(cond,rep).all(:,iJ); % joint angle data
      % find start of swing
      min_loc = round(angles(ifly).Leg(leg).step(cond,rep).min_loc);
      % find start of stance
      max_loc = round(angles(ifly).Leg(leg).step(cond,rep).max_loc);
% CONTROL      
      % --- control SWING->STANCE ---- %
      subplot(nrows,ncols,1); hold all
      swing_strt = round(min_loc(min_loc<ROI(end)));
      stance_strt = round(max_loc(max_loc<ROI(end)));
      nSteps = length(swing_strt)-1; % number of steps
      for n = 1:nSteps
       % check for max values between swing points
        step_range = swing_strt(n):swing_strt(n+1);
        step = raw(step_range);
        step = step-(step(1)); % normalize data
        %save data to matrix
        input = MT;
        input(1:length(step)) = step;
        control_stance(:,end+1) = input;
        
        % plot the trials:
        plot(step, 'color', Color('grey'), 'linewidth', LW)
      end
       % --- control STANCE->SWING ---- %
      subplot(nrows,ncols,4); hold all
      nSteps = length(stance_strt)-1; % number of steps
      for n = 1:nSteps
       % check for max values between swing points
        step_range = stance_strt(n):stance_strt(n+1);
        step = raw(step_range);
        step = step-(step(1)); % normalize data
        %save data to matrix
        input = MT;
        input(1:length(step)) = step;
        control_swing(:,end+1) = input;
        % plot the trials:
        plot(step, 'color', Color('grey'), 'linewidth', LW)
      end
% STIM      
      norm_range = round(mean(diff(min_loc))) + 15;
      a = min_loc(min_loc>light_on);
      b = max_loc(max_loc>light_on);
      try a(1) < b(1); catch; continue; end % no steps - exclude
      %in STANCE
      if a(1) < b(1) % stance->swing transition at light-on
          step_range = a(1):a(2);% swing start to next swing start
          pcolor = Color('Turquoise');
          if length(step_range)>norm_range % include standard range:
              step_range = a(1):a(1)+norm_range;
              pcolor = Color('grey');
          end
          % subplot activation:
          subplot(nrows,ncols,2); hold all
          step = raw(step_range);
          step = step-(step(1)); % normalize data
          plot(step, 'color', pcolor, 'linewidth', LW)
          %save data to matrix
          input = MT;
          input(1:length(step)) = step;
          stance(:,end+1) = input;
      %in SWING    
      elseif a(1) > b(1) % swing->stance transition at light-on
          step_range = b(1):b(2);% swing start to next swing start
          pcolor = Color('mediumpurple');
          if length(step_range)>norm_range % include standard range:
              step_range = b(1):b(1)+norm_range;
              pcolor = Color('grey');
          end
          % subplot activation:
          subplot(nrows,ncols,5); hold all
          step = raw(step_range);
          step = step-(step(1)); % normalize data
          plot(step, 'color', pcolor, 'linewidth', LW)
          %save data to matrix
          input = MT;
          input(1:length(step)) = step;
          swing(:,end+1) = input;
      end
   end
end; input = [];

% Plot the avgs and then the overlays:
% CONTROL STANCE
subplot(nrows,ncols,1); hold all
control_stance(:,1) = [];
input = nanmean(control_stance,2);
input(isnan(input)) = [];
% only plot to point that 25% of steps are done:
if size(control_stance,2)> 1
    loc = sum(isnan(control_stance),2)>round(size(control_stance,2)*0.25); 
    input(loc(1:length(input))) = [];
end
plot(input, 'linewidth', LW+2, 'color', 'k')
subplot(nrows,ncols,3); hold all
plot(input, 'linewidth', LW+2, 'color', 'k')
% CONTROL SWING
subplot(nrows,ncols,4); hold all
control_swing(:,1) = [];
input = nanmean(control_swing,2);
input(isnan(input)) = [];
% only plot to point that 25% of steps are done:
if size(control_swing,2)> 1
    loc = sum(isnan(control_swing),2)>round(size(control_swing,2)*0.25); 
    input(loc(1:length(input))) = [];
end
plot(input, 'linewidth', LW+2, 'color', 'k')
subplot(nrows,ncols,6); hold all
plot(input, 'linewidth', LW+2, 'color', 'k')


% SWING:
subplot(nrows,ncols,5); hold all
swing(:,1) = [];
input = nanmedian(swing,2);
input(isnan(input)) = [];
% only plot to point that 25% of steps are done:
if size(swing,2)> 1
    loc = sum(isnan(swing),2)>round(size(swing,2)*0.25); 
    input(loc(1:length(input))) = [];
end
plot(input, 'linewidth', LW+2, 'color', Color('Purple'))
subplot(nrows,ncols,6); hold all
plot(input, 'linewidth', LW+2, 'color', Color('Purple'))


% STANCE:
subplot(nrows,ncols,2); hold all
stance(:,1) = [];
input = nanmedian(stance,2);
input(isnan(input)) = [];
% only plot to point that 25% of steps are done:
if size(stance,2)> 1
    loc = sum(isnan(stance),2)>round(size(stance,2)*0.25); 
    input(loc(1:length(input))) = [];
end
plot(input, 'linewidth', LW+2, 'color', Color('Teal'))
subplot(nrows,ncols,3); hold all
plot(input, 'linewidth', LW+2, 'color', Color('Teal'))

% ----- plot labels -----:
%control - stance
subplot(nrows,ncols,1)
title({param.conds_matrix(cond).label;'Control Stance-Swing'})
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])

%first step in stim, fly was in stance
subplot(nrows,ncols,2)
title('Stim step: Stance')
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])

%stance overlay
subplot(nrows,ncols,3)
title('Stance-Swing step')
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])

%control - swing
subplot(nrows,ncols,4)
title('Stim step: Swing')
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])

%first step in stim, fly was in swing
subplot(nrows,ncols,5)
title('Stim step: Swing')
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])

%swing overlay
subplot(nrows,ncols,6)
title('Swing-Stance step')
xlabel('time (fps)')
ylabel([Joints{iJ} ' \Delta angle (\circ)'])


fig_title = [FilePath.structure_name ' ' param.conds_matrix(cond).label ' isolated first stim step'];
save_figure(fig, [fig_dir '\' fig_title]);

clearvars('-except',initial_vars{:})
% 
%% POSITION BASED DATA ANALYSIS


































