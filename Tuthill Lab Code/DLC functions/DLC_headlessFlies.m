
% Use to look at a series of headless flies on or off the ball (assumes
% that they aren't doing much walking so doesn't focus on gait parameters)

% Use this after you have run DLC_Step1 for importing the data from Anipose

clearvars
clc
close all
warning('off')

% TODO section:
% Set base file paths:
fileroot = 'D:\Tuthill Lab Work\Anipose\';                  %local location for step1 processed Anipose data
googledrive = 'G:\My Drive\Evyn\Data\FicTrac Raw Data\';     %google drive root folder
% ----------------------
folder_date = '1.10.20';                                    %folder to analyze
% ----------------------

local_save = true;                                          %save data in the local PC folder?
online_save = true;                                         %save data to google drive folder?

%% Load data: 

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

param = anipose.param;
num = NUM(param);
angles = anipose.angles;
pose_3d = anipose.pose_3d;
fps = param.Basler_fps;

initial_vars = who; initial_vars{end+1} = 'initial_vars';

initial_vars{end+1} = 'group';
initial_vars{end+1} = 'behavior';
initial_vars{end+1} = 'gait';

%% Behavior predictions based on joint angles: 

% make leg behavior and fly behavior predictions
ROI = 1:param.basler_delay*fps; % sort behavior by all frames before the laser
[behavior, figHandle] = DLC_behaviorpredictions(angles,param,ROI);
save_figure(figHandle, [fig_root, 'Behavior Predictions'], '-png');

disp('Finished predicting behavior!')

%% PLOT joint angle time course
light_on = param.basler_delay*fps;

% ---- input -----
G1.conds = [1:7:28];        % group 1 conditions: control
G2.conds = [1:28];        % group 2 conditions: stim
G1.color = {'Black', 'Darkslategrey', 'darkgrey'};         % group 1 color choice
G2.color = {'Indigo', 'Orangered', 'Teal'};     % group 1 color choice
G1.ROI = 1:601;             % all time
G2.ROI = 1:601;             % all time
for leg = 1:6
   [G1.data,G2.data] = deal([]);
end
% ----------------
G2.conds(G1.conds) = [];

% Group 1 data:
idx = 0;
for n = 1:length(G1.conds) %all user defined conditions
  cond = G1.conds(n); %specific condition for this cycle
  for rep = 1:num.reps
    % check for possible data: * can later add filters here like
    % 'walking only' etc
    if behavior.stationary(cond,rep)==true %only look at stationary flies
      idx = idx+1;
      for leg = 1:6   %all legs
        for iJ = 1:3
          input = angles.Leg(leg).data(cond,rep).all; %(trial,joint,time)
          G1.data(leg,iJ).all(:,idx) = input(G1.ROI,iJ);
        end
      end
    end
  end
end

% Group 2 data:
idx = 0;
for n = 1:length(G2.conds) %all user defined conditions
  cond = G2.conds(n); %specific condition for this cycle
  for rep = 1:num.reps
    % check for possible data: * can later add filters here like
    % 'walking only' etc
    if behavior.stationary(cond,rep)==true %only look at stationary flies
      idx = idx+1;
      for leg = 1:6   %all legs
        for iJ = 1:3
          input = angles.Leg(leg).data(cond,rep).all; %(trial,joint,time)
          G2.data(leg,iJ).all(:,idx) = input(G2.ROI,iJ);
        end
      end
    end
  end
end

% Pull the change in angle & averages: 
for iJ = 1:3
  for leg = 1:6
    % change in joint angle:
    offset = (G1.data(leg,iJ).all(light_on,:)); %joint angle at light on
    G1.data(leg,iJ).change = G1.data(leg,iJ).all - offset;
    offset = (G2.data(leg,iJ).all(light_on,:)); %joint angle at light on
    G2.data(leg,iJ).change = G2.data(leg,iJ).all - offset;
    
    % find the average joint angle:
    G1.data(leg,iJ).avg = nanmean(G1.data(leg,iJ).change,2);
    G2.data(leg,iJ).avg = nanmean(G2.data(leg,iJ).change,2);
    
    % find the err in joint angle:
    G1.data(leg,iJ).err = sem(G1.data(leg,iJ).change,2);
    G2.data(leg,iJ).err = sem(G2.data(leg,iJ).change,2);
  end
end

% Plot the data as a time course:
subplot_loc = [1 3 5 2 4 6]; % organize left legs on leg and right on right

fig = getfig('',1);
for leg = 1:6
    subplot(3,2,subplot_loc(leg))
    hold on; axis normal
    x = -0.5:1/300:1.5;
    
    for iJ = 1:3
        %plot G1 data:
        avg = G1.data(leg,iJ).avg;
        err = G1.data(leg,iJ).err;
        % shaded error regions
        fill_data = error_fill(x, avg, err);
        h = fill(fill_data.X, fill_data.Y, Color(G1.color{iJ}), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(x, avg, 'linewidth', 1, 'color', Color(G1.color{iJ}))
        %plot G2 data:
        avg = G2.data(leg,iJ).avg;
        err = G2.data(leg,iJ).err;
        % shaded error regions
        fill_data = error_fill(x, avg, err);
        h = fill(fill_data.X, fill_data.Y, Color(G2.color{iJ}), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
        plot(x, avg, 'linewidth', 1, 'color', Color(G2.color{iJ}))
    end
    % plot line for laser:
    y = rangeLine(fig);
    plot([0, 0.72], [y,y], 'linewidth',3, 'color', 'g')
    
    % labels
    xlabel('time (s)')
    ylabel('\Delta joint angle (\circ)')
    title(['Leg ' anipose.leg_labels{leg}])
end
 
save_figure(fig, [fig_root, ' joint angles for all legs'],'-png');

clearvars('-except',initial_vars{:})


%% PLOT All joint angles 3D scatter plot *new thing*
light_on = param.basler_delay*fps;

% ---- input -----
G1.conds = [1:7:28];        % group 1 conditions: control
G2.conds = [1:28];        % group 2 conditions: stim
G1.color = 'Black';         % group 1 color choice
G2.color = 'Orangered';     % group 1 color choice
G1.ROI = 1:600;             % all time
G2.ROI = 1:600;             % all time
for leg = 1:6
   [G1.data(leg).all,G2.data(leg).all] = deal([]);
end
% ----------------
G2.conds(G1.conds) = []; 

% Group 1 data:
for n = 1:length(G1.conds) %all user defined conditions
  cond = G1.conds(n); %specific condition for this cycle
  for rep = 1:num.reps
    % check for possible data: * can later add filters here like
    % 'walking only' etc
    if behavior.stationary(cond,rep)==true %only look at stationary flies
      for leg = 1:6   %all legs
        input = angles.Leg(leg).data(cond,rep).all; %(trial,joint,time)
        G1.data(leg).all = [G1.data(leg).all; input(G1.ROI,:)];
      end
    end
  end
end

% Group 2 data:
for n = 1:length(G2.conds) %all user defined conditions
  cond = G2.conds(n); %specific condition for this cycle
  for rep = 1:num.reps
    % check for possible data: * can later add filters here like
    % 'walking only' etc
    if behavior.stationary(cond,rep)==true %only look at stationary flies
      for leg = 1:6   %all legs
        input = angles.Leg(leg).data(cond,rep).all; %(trial,joint,time)
        G2.data(leg).all = [G2.data(leg).all; input(G2.ROI,:)];
      end
    end
  end
end

Joints = anipose.Joints;
SZ = 30;
subplot_loc = [1 3 5 2 4 6]; % organize left legs on leg and right on right

% Plot the data into 3d space: 
fig = getfig('Joint angles by leg',1);
for leg = 1:6
    subplot(3,2,subplot_loc(leg))
    input = G1.data(leg).all;
    scatter3(input(:,1), input(:,2), input(:,3), SZ, Color(G1.color),'filled')
    hold on
    input = G2.data(leg).all;
    scatter3(input(:,1), input(:,2), input(:,3), SZ, Color(G2.color),'filled')
    
    % labels
    grid on
    xlabel(Joints{1}); ylabel(Joints{2}); zlabel(Joints{3});
    title(anipose.leg_labels{leg})
end

save_figure(fig, [fig_root, ' three joint angles scatterplot'],'-png');

clearvars('-except',initial_vars{:})















