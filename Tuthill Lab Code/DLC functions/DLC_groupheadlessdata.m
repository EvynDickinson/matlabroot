
clear; clc


%% create list of flies to include: 
% Set the base directory for the data files: 
fileroot = 'D:\Tuthill Lab Work\Anipose\';  % local or google
output_root = 'C:\matlabroot\DLC\'; %TODO

a = questdlg('Open previously generated file?');
switch a
    case 'No' 
        FilePath = DLC_select_flies('-inpath', fileroot, '-outpath', output_root, '-export', false);
        % load formerly-generated structures: 

        % make a folder for the analysis:
        fig_dir = [output_root, FilePath.structure_name];
        if ~exist(fig_dir, 'dir')
            mkdir(fig_dir)
        end 
        % Load the individual fly anipose data previously generated with 'DLC_Step1_loaddata.m'
        for ifly = 1:size(FilePath.locations,1)
            flyID = FilePath.locations{ifly,4};
            root_dir = [fileroot, FilePath.locations{ifly,2},'\',FilePath.locations{ifly,3}];
            temp = load([root_dir, '\Analysis\' flyID '_AniposeData']);
            fly(ifly) = temp.anipose;
        end
    case 'Yes'
        [file,path] = uigetfile;
        load(fullfile(path, file))
        % make a folder for the analysis:
        fig_dir = [output_root, FilePath.structure_name];
        if ~exist(fig_dir, 'dir')
            mkdir(fig_dir)
        end 
end

% save the grouped fly data:
save([fig_dir '\' FilePath.structure_name '_anipose']);

% general parameters:
nflies = size(FilePath.locations,1);
fps = fly(1).param.Basler_fps;

disp('Done loading flies')
initial_vars = who; initial_vars{end+1} = 'initial_vars';

initial_vars{end+1} = 'behavior';
initial_vars{end+1} = 'gait';

% Predict behavior for all flies:
for ifly = 1:nflies
    param = fly(ifly).param;
    ROI = 1:param.basler_delay*fps; % sort behavior by all frames before the laser
    [fly(ifly).behavior,f] = DLC_behaviorpredictions(fly(ifly).angles,param,ROI);
    close(f)
end
disp('Behavior sorted')

%% Look at time course for joint angles:
param = fly(1).param;
num = NUM(param);
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
  for ifly = 1:nflies
      for rep = 1:num.reps
        % check for possible data: * can later add filters here like
        % 'walking only' etc
        if fly(ifly).behavior.stationary(cond,rep)==true %only look at stationary flies
          idx = idx+1;
          for leg = 1:6   %all legs
            for iJ = 1:3
              input = fly(ifly).angles.Leg(leg).data(cond,rep).all; %(trial,joint,time)
              G1.data(leg,iJ).all(:,idx) = input(G1.ROI,iJ);
            end
          end
        end
      end
  end
end
% Group 2 data:
idx = 0;
for n = 1:length(G2.conds) %all user defined conditions
  cond = G2.conds(n); %specific condition for this cycle
  for ifly = 1:nflies
      for rep = 1:num.reps
        % check for possible data: * can later add filters here like
        % 'walking only' etc
        if fly(ifly).behavior.stationary(cond,rep)==true %only look at stationary flies
          idx = idx+1;
          for leg = 1:6   %all legs
            for iJ = 1:3
              input = fly(ifly).angles.Leg(leg).data(cond,rep).all; %(trial,joint,time)
              G2.data(leg,iJ).all(:,idx) = input(G2.ROI,iJ);
            end
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
YL = [-10,10];

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
    ylim(YL)
    y = rangeLine(fig);
    plot([0, 0.72], [y,y], 'linewidth',3, 'color', 'g')
    
    % labels
    xlabel('time (s)')
    ylabel('\Delta joint angle (\circ)')
    title(['Leg ' fly(ifly).leg_labels{leg}])
    
end
 
save_figure(fig, [fig_dir, ' joint angle time course all legs'],'-png');

clearvars('-except',initial_vars{:})

%% 



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




