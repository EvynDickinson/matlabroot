
% Basic figures for things like step frequency, stride length, swing speed,
% joint angle distributions.
% use data structures preloaded from DLC_Step1 and DLC_Step2

clearvars
clc
close all
warning('off')

% TODO section:
% Set base file paths:
fileroot = 'D:\Tuthill Lab Work\Anipose\';                  %local location for step1 processed Anipose data
googledrive = 'G:\My Drive\Evyn\Data\FicTrac Raw Data\';     %google drive root folder
% ----------------------
folder_date = '6.11.19';                                    %folder to analyze
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
    warndlg('No file matching AniposeData found'); return
end
try
    load([pathway,'\' FlyToLoad{1} '\Analysis\' flyID ' StepData.mat'])
catch
    warndlg('No file matching StepData found'); return
end

fig_root = [pathway,'\' FlyToLoad{1} '\Analysis\' flyID];

disp('Fly selected: ')
disp(FlyToLoad)

param = anipose.param;
num = NUM(param);
angles = anipose.angles;
pose_3d = anipose.pose_3d;
fps = param.Basler_fps;

% save variables for an easy clear between analysis / plotting
initial_vars = who; initial_vars{end+1} = 'initial_vars';

initial_vars{end+1} = 'group';
initial_vars{end+1} = 'behavior';
initial_vars{end+1} = 'gait';
disp('Finished load and setting params')

%% Find step frequency for all feasible trials

% for each trial, pull the step frequency over the full experiement: 
for cond = 1:num.conds
  for rep = 1:num.reps
    % pull the logical of swing|stance per frame
    raw = step.stance(cond,rep).loc'; %flip orientation
    if isempty(raw); continue; end % filter out non-sphere-fit trials
    freq = [];
    for leg = 1:6 %each leg
        input = raw(:,leg);
        stance_strt = find(diff(input)==1);
        swing_strt = find(diff(input)==-1);
        f = 1./([stance_strt(1); diff(stance_strt)]./fps); % step frequency at each stance start
        
        % extrapolate and fill the gaps to have a fps matrix: 
        Y = nan(param.basler_length*fps+1,1);
        X = [1; stance_strt; param.basler_length*fps+1];
        for ii = 1:length(X)-1
            if ii == length(X)-1; val = f(ii-1); 
            else val = f(ii); end
            strt = X(ii);
            fin = X(ii+1);
            Y(strt:fin-1) = val*ones(fin-strt,1);
        end
        Y(end) = Y(end-1);  
        y = smooth(Y,fps/10,'lowess'); 
        
        % save data into struct:
        freq(leg).stance_strt = stance_strt;
        freq(leg).swing_strt = swing_strt;
        freq(leg).freq = y;
        freq(leg).rawFreq = Y;
        freq(leg).inst_freq = [stance_strt,f];
        F_all(:,leg) = y;
    end
    % assign the frequency data to the Gait structure:
    gait.Freq(cond,rep).data = freq;
    gait.Freq(cond,rep).f_all = F_all;
    gait.Freq(cond,rep).f_avg = mean(F_all,2);
  end
end
disp('Finished calculating step frequency')
clearvars('-except',initial_vars{:})

%% PLOT step frequency per SINGLE VIDEO
% note: you should ignore the first bit as there's an edge effect in
% predicting the frequency

% ---- select trial -----
cond = 28;
rep = 3;
% -----------------------

x = -(param.basler_delay):1/fps:param.basler_length-(param.basler_delay);
light_on = (param.conds_matrix(cond).opto);
freq = [];

LW = 1;
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};

fig = getfig('Step Frequency',1); hold all
for leg = 1:6
    y = gait.Freq(cond,rep).data(leg).freq;
    plot(x, y, 'color', Color(leg_colors{leg}), 'linewidth', LW)
    freq(:,leg) = y;
end
% plot the average: 
plot(x, median(freq,2), 'color', 'k', 'linewidth', LW+1)
% laser region lines:
y1 = rangeLine(fig);
plot([0,light_on], [y1,y1],'color',Color(param.LED_light), 'linewidth', 5)
% labels
xlabel('time (s)')
ylabel('Step frequency (hz)')
fig_title = {param.cross; flyID};
fig_title = strrep(fig_title,'_','-');
title(fig_title)
legend([anipose.leg_labels, {'avg'}])

fig = formatFig(fig, false);
separateAxes

% save option
save_figure(fig, [fig_root, 'R' num2str(rep) 'C' num2str(cond) ' Raw step frequency']);
%note: if you want a clean figure for powerpoint and don't want to use
%illustrator to adjust it at all, add: ,'-png' at the end of the
%save_figure function (e.g. save_figure(fighandle, savename, '-png'))
clearvars('-except',initial_vars{:})

%% Calculate the stride length:
LP = 5;
% for each trial, pull the step frequency over the full experiement: 
for cond = 1:num.conds
  for rep = 1:num.reps
    % pull the logical of swing|stance per frame
    raw = step.stance(cond,rep).loc'; %flip orientation
    if isempty(raw); continue; end % filter out non-sphere-fit trials
    stride = []; %reset the stride struct to fill
    for leg = 1:6 %each leg
        [stance_strt,swing_strt,nstrides] = deal([]);
        input = raw(:,leg);
        stance_strt = find(diff(input)==1);
        swing_strt = find(diff(input)==-1);
        % determine which comes first: swing or stance:
        if stance_strt(1)<swing_strt(1) %stance first
            nstrides = min([length(stance_strt),length(swing_strt)]);
            stride(leg).loc = [stance_strt(1:nstrides),swing_strt(1:nstrides)];
        else
            nstrides = min([length(stance_strt), length(swing_strt)-1]);
            stride(leg).loc = [stance_strt(1:nstrides),swing_strt(2:1+nstrides)];
        end
        stride(leg).nstrides = nstrides;
        % pull the tarsus positions
        data = pose_3d.Leg(LP,leg).data{cond,rep}.all;

        % find the tarsus position for each stride:
        stride(leg).start_pos = data(stride(leg).loc(:,1),:);
        stride(leg).end_pos = data(stride(leg).loc(:,2),:);
        
        % find the euclidian distance between the start and end of the stride
        for ii = 1:nstrides
            stride(leg).length(ii,1) = pdist2(stride(leg).start_pos(ii,:),stride(leg).end_pos(ii,:));
        end
    end
    % find the average stride leg
    % assign the frequency data to the Gait structure:
    gait.Stride(cond,rep).data = stride;
  end
end
disp('Finished calculating stride length')

clearvars('-except',initial_vars{:})

%% Plot stride lengths for a given condition
% this is a bit of a strange plot: it plots the stride length (y) as a
% constant value for the duration of each stride, but plots nothing during
% the swing period -- you could redo this with the stride length lasting
% for the entire stance-to-stance cycle, but if you leave a gap you get
% some visual information about the swing duration (time-wise)

% ---- input ----
cond = 28;
rep = 3;
% ---------------

x = -(param.basler_delay):1/fps:param.basler_length-(param.basler_delay);
light_on = (param.conds_matrix(cond).opto);
freq = [];

LW = 3;
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};

fig = getfig('Stride Length',1);
hold all
MT = nan(6,length(x));
for leg = 1:6
   loc = gait.Stride(cond,rep).data(leg).loc;
   X = x(loc); 
   y = gait.Stride(cond,rep).data(leg).length; 
   Y = [y,y];
   % plot line for each stride
   for s = 1:size(Y,1)
       plot(X(s,:),Y(s,:), 'linewidth', LW, 'color', Color(leg_colors{leg}))
       MT(leg,loc(s,1):loc(s,2)) = y(s);
   end
end
% plot the average line: 
avg = smooth(nanmean(MT),10);
plot(x,avg, 'linewidth', LW, 'color', 'k')

% laser region lines:
vline([0,light_on], 'g-')

% labels
xlabel('time (s)')
ylabel('Stride length (pixels?)')
fig_title = {param.cross; [flyID 'R' num2str(rep) 'C' num2str(cond)]};
fig_title = strrep(fig_title,'_','-');
title(fig_title)

fig = formatFig(fig);
separateAxes

% save option
save_figure(fig, [fig_root, 'R' num2str(rep) 'C' num2str(cond) ' Raw stride length']);

clearvars('-except',initial_vars{:})

%% Joint angles sections

%% PLOT *All* legs *single* joint angle histogram
% compare the joint angles across three joints between two different ROIs 
%e.g. control joint angle distribution vs the stimulus period
light_on = param.basler_delay*fps;

% ---- input -----
iJ = 2;                 % joint angle to show
G1.conds = [7,14];      % group 1 conditions to lump together, can be as many as desired
G2.conds = [7,14];      % group 2 conditions to lump together
G1.color = 'Black';     % group 1 color choice
G2.color = 'Orangered'; % group 1 color choice
G1.ROI = 1:light_on;    % control ROI *current
G2.ROI = light_on:(light_on+(0.5*fps))-1; %0.5sec starting at stim on
[G1.data,G2.data] = deal([]);
% ----------------

% Group 1 data:
for n = 1:length(G1.conds) %all user defined conditions
  cond = G1.conds(n); %specific condition for this cycle
  for rep = 1:num.reps
    % check for possible data: * can later add filters here like
    % 'walking only' etc
    if ~isempty(step.ballFit(cond,rep).fit) && step.ballFit(cond,rep).fit==true
      for leg = 1:6   %all legs
        all(:,leg) = angles.Leg(leg).data(cond,rep).all(:,iJ); %(trial,joint,time)
      end
      G1.data = [G1.data; all(G1.ROI,:)];
    end
  end
end

% Group 2 data:
for n = 1:length(G2.conds) %all user defined conditions
  cond = G2.conds(n); %specific condition for this cycle
  for rep = 1:num.reps
    % check for possible data: * can later add filters here like
    % 'walking only' etc
    if ~isempty(step.ballFit(cond,rep).fit) && step.ballFit(cond,rep).fit==true
      for leg = 1:6   %all legs
        all(:,leg) = angles.Leg(leg).data(cond,rep).all(:,iJ); %(trial,joint,time)
      end
      G2.data = [G2.data; all(G2.ROI,:)];
    end
  end
end


% Plot the data:
subplot_loc = [1 3 5 2 4 6]; % organize left legs on leg and right on right
Edges = 0:10:180; % edges for histogram bins
fig = getfig('Joint Angle Histogram',1);
for leg = 1:6   % subplot for each leg
    subplot(3,2,subplot_loc(leg)); hold on
    % plot group 1:
    h = histogram(G1.data(:,leg),Edges);
    h.FaceColor = Color(G1.color); 
    % plot group 2:
    h = histogram(G2.data(:,leg),Edges);
    h.FaceColor = Color(G2.color); 
    
    % plot labels:
    xlabel('Angle')
    ylabel([anipose.Joints{iJ} ' frame count'])
    title(anipose.leg_labels{leg})
    set(gca, 'TickDir', 'out')
    xlim([0,180])
end

% save option ** this will over-right itself each time so add a tag if
% desired ** 
save_figure(fig, [fig_root, ' ' anipose.Joints{iJ} ' joint angle histogram']);

clearvars('-except',initial_vars{:})

%% PLOT All joint angles 3D scatter plot *new thing*
% compare the joint angles across three joints between two different ROIs 
%e.g. control joint angle distribution vs the stimulus period
light_on = param.basler_delay*fps;

% ---- input -----
G1.conds = [7,14];      % group 1 conditions to lump together, can be as many as desired
G2.conds = [7,14];      % group 2 conditions to lump together
G1.color = 'Black';     % group 1 color choice
G2.color = 'Orangered'; % group 1 color choice
G1.ROI = 1:light_on;    % control ROI *current
G2.ROI = light_on:(light_on+(0.5*fps))-1; %0.5sec starting at stim on
for leg = 1:6
   [G1.data(leg).all,G2.data(leg).all] = deal([]);
end
% ----------------

% Group 1 data:
for n = 1:length(G1.conds) %all user defined conditions
  cond = G1.conds(n); %specific condition for this cycle
  for rep = 1:num.reps
    % check for possible data: * can later add filters here like
    % 'walking only' etc
    if ~isempty(step.ballFit(cond,rep).fit) && step.ballFit(cond,rep).fit==true
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
    if ~isempty(step.ballFit(cond,rep).fit) && step.ballFit(cond,rep).fit==true
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

%% PLOT Average joint angle between two ROIs
light_on = param.basler_delay*fps;

% ---- input -----
paired = true;          % paired data? (aka same conds for G1&G2)
G1.conds = [7,14];      % group 1 conditions to lump together, can be as many as desired
G2.conds = [7,14];      % group 2 conditions to lump together
G1.color = 'Black';     % group 1 color choice
G2.color = 'Orangered'; % group 1 color choice
G1.ROI = 1:light_on;    % control ROI *current
G2.ROI = light_on:(light_on+(0.5*fps))-1; %0.5sec starting at stim on
% ----------------

% Group 1 data:
idx = 0;
for n = 1:length(G1.conds) %all user defined conditions
  cond = G1.conds(n); %specific condition for this cycle
  for rep = 1:num.reps
    % check for possible data: * can later add filters here like
    % 'walking only' etc
    if ~isempty(step.ballFit(cond,rep).fit) && step.ballFit(cond,rep).fit==true
      idx = idx+1;
      for leg = 1:6   %all legs
        input = angles.Leg(leg).data(cond,rep).all; %(trial,joint,time)
        for iJ = 1:3  %All joints
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
    if ~isempty(step.ballFit(cond,rep).fit) && step.ballFit(cond,rep).fit==true
      idx = idx+1;
      for leg = 1:6   %all legs
        input = angles.Leg(leg).data(cond,rep).all; %(trial,joint,time)
        for iJ = 1:3  %All joints
          G2.data(leg,iJ).all(:,idx) = input(G2.ROI,iJ);
        end
      end
    end
  end
end

% Pull the averages: 
for iJ = 1:3
  for leg = 1:6
    G1.data(leg,iJ).avg = nanmean(G1.data(leg,iJ).all,1);
    G2.data(leg,iJ).avg = nanmean(G2.data(leg,iJ).all,1);
  end
end

% Plotting information:
%3 graphs -- one per joint
SZ = 30;
LW = 1;
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};

%Plot the data:
fig = getfig('Avg Joint Angle',1);
for iJ = 1:3
    subplot(3,1,iJ)
    hold on
    % plot the avg joint angle+SEM
    for leg = 1:6
        %group 1
        x1 = (leg*2-1)*ones(1,length(G1.data(leg,iJ).avg));
        scatter(x1, G1.data(leg,iJ).avg, SZ, Color(leg_colors{leg}), 'filled')
        %group 2
        x2 = (leg*2)*ones(1,length(G1.data(leg,iJ).avg));
        scatter(x2, G2.data(leg,iJ).avg, SZ, Color(leg_colors{leg}), 'filled')
        % if paired data, plot connecting line by trial
        if paired==true
          for ii = 1:length(x1)
            plot([x1(ii),x2(ii)], [G1.data(leg,iJ).avg(ii), G2.data(leg,iJ).avg(ii)], ...
                'color', Color(leg_colors{leg}), 'linewidth', LW)
          end
        end
    end
    % labels and formatting:
    xlim([0,13])
    ylim([0,180])
    ylabel(['\Delta ' anipose.Joints{iJ} ' angle (\circ)'])
    xlabel('Legs')
end

fig = formatFig(fig, true, [3,1]);
for iJ = 1:3
    subplot(3,1,iJ); separateAxes
end

save_figure(fig, [fig_root, ' Avg joint angle']);

clearvars('-except',initial_vars{:})



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    









