
% GROUP FLIES TOGETHER USING STRUCTURES DESCRIBED IN YOUR FLY SUMMARY
clear; 
clc 
warning('off')

%% create list of flies to include: 
% Set the base directory for the data files: 
fileroot = 'D:\Tuthill Lab Work\Anipose\';  % local or google
output_root = 'C:\matlabroot\DLC\';

a = questdlg('Open previously generated file?','','Yes', 'No', 'Cancel', 'No');
switch a
    case 'No' 
        FilePath = DLC_select_flies('-inpath', fileroot, '-outpath', output_root, '-export', false); %TODO update paths
        % load formerly-generated structures: 

        % make a folder for the analysis:
        fig_dir = [output_root, FilePath.structure_name];
        if ~exist(fig_dir, 'dir')
            mkdir(fig_dir)
        end 
        % Load the individual fly anipose data previously generated with DLC_Step1
        % & DLC_Step2
        fprintf('\nFinished loading: \n')
        for ifly = 1:size(FilePath.locations,1)
            flyID = FilePath.locations{ifly,4};
            root_dir = [fileroot, FilePath.locations{ifly,2},'\',FilePath.locations{ifly,3}];
            temp = load([root_dir, '\Analysis\' flyID '_AniposeData']);
            
            % save anipose data:
            fly(ifly) = temp.anipose;
            
            % save step data
            try %not all flies have data analyzed here
                temp = load([root_dir, '\Analysis\' flyID ' StepData']);
                step(ifly) = temp.step;
            catch
            end
            
            % save sphere fit data:
            try %not all flies have data analyzed here
                temp = load([root_dir, '\Analysis\' flyID ' Sphere Fit data']);
                sphereFit(ifly) = temp;
            catch
            end
            disp(ifly)
            clear temp
        end
        % general parameters:
        nflies = size(FilePath.locations,1);
        fps = fly(1).param.Basler_fps;

        disp('Done loading flies')
        initial_vars = who; initial_vars{end+1} = 'initial_vars';

        % Sort the behavior again of the flies:
        initial_vars{end+1} = 'behavior';
%         initial_vars{end+1} = 'gait';

        % Predict behavior for all flies:
        fprintf('\nSorting behavior: \n')
        for ifly = 1:nflies
            param = fly(ifly).param;
            ROI = 1:param.basler_delay*fps; % sort behavior by all frames before the laser
            [fly(ifly).behavior,f] = DLC_behaviorpredictions(fly(ifly).angles,param,ROI);
            close(f)
            disp(ifly)
        end
        disp('Behavior sorted')
        % save the grouped fly data:
        save([fig_dir '\' FilePath.structure_name '_KinematicData']);
        disp('Finished saving structure')
        
    case 'Yes'
        [file,path] = uigetfile;
        load(fullfile(path, file))
        % make a folder for the analysis:
        fig_dir = [output_root, FilePath.structure_name];
        if ~exist(fig_dir, 'dir')
            mkdir(fig_dir)
        end 
    case 'Cancel'
        return
end

clearvars('-except',initial_vars{:})


%% Heatmap of joint angle change for all legs and joints:
behaviorState = 'all';  % behavior group to isolate e.g. 'walking', 'stationary','all'
color_opt = false;      % false = white background, true = black background formatting
param = fly(1).param;   % assumes all flies have the same params (#conds,etc)
light_on = param.basler_delay*fps;
num = NUM(param);       % pull num reps/conds from a fly
% set colorLimit to nan if you want it automatically selected for you
colorLimit = nan;       % normalizing range for color (e.g. -15 deg:+15 deg range)
stats_text = true;       % display the p-value on the heatmap squares 
fileType = '-png';      % figure save type (e.g. '-png' or '-pdf')
nJoints = 1:3;          % number of joints to cycle through
nLegs = 1:6;            % legs to use

% define the two ROIs for comparison
condList = [4,11];  % what conditions do you want to group?

% auto will be based on the first condition listed...
laser = param.conds_matrix(condList(1)).opto;
switch questdlg(['Auto fit ROI to control vs stim of first cond (' num2str(laser) ' ms)?'])
    case 'Yes' 
        ROI_1 = 1:light_on;
        ROI_2 = light_on+1:light_on+1+round(laser*fps);
    case 'No'
        ROI_1 = 1:150;     % first 0.5 sec (control)
        ROI_2 = 151:200;   % 720 ms of stimulus   
end


% Pull the joint angle data for each 

% Data extraction
for ifly = 1:nflies
%     % skip flies without sphere fits
%     if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(condList)
        cond = condList(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
%             % check for possible data:
%             if isempty(step(ifly).ballFit(cond,rep).fit)...
%                || step(ifly).ballFit(cond,rep).fit==false
%                  continue;
%             end
            idx = idx+1;
            % add data for each leg into the G1 data struct
            for iJ = nJoints
              for leg = nLegs
                % pull control data:
                control = fly(ifly).angles.Leg(leg).data(cond,rep).all(ROI_1,iJ);
                % pull stimulus data:
                stim = fly(ifly).angles.Leg(leg).data(cond,rep).all(ROI_2,iJ);
                G(leg,iJ).fly(ifly).ROIs(idx,1) = nanmean(control);
                G(leg,iJ).fly(ifly).ROIs(idx,2) = nanmean(stim);
                % find the change:
                G(leg,iJ).fly(ifly).all(idx) = nanmean(stim)-nanmean(control);
              end
            end
            
        end
    end
end


% Pool data & find averages within a group:
for leg = nLegs
    for iJ = nJoints
        idx = 0;
        for ifly = 1:length(G(leg,iJ).fly)
            % screen out nonmatching data
            temp = G(leg,iJ).fly(ifly).all;
            if isempty(temp); continue; end
            idx = idx+1;   
            raw(idx) = nanmean(temp);
        end
        % find the cross-fly average
        data(leg,iJ) = nanmean(raw);
        raw = [];
    end
end

% ----------------------------------------
% determine significance using a paired t-test w/mulitple corrections
for leg = nLegs
    for iJ = nJoints
        idx = 0;
        for ifly = 1:length(G(leg,iJ).fly)
            % screen out nonmatching data
            temp = G(leg,iJ).fly(ifly).ROIs;
            if isempty(temp); continue; end
            idx = idx+1;   
            raw(idx,1:2) = nanmean(temp,1); 
        end
        if idx<2 %no stats are possible here given sample size
            disp(['Sample size under stats minimum, leg-' num2str(leg) ' joint-' num2str(iJ)])
        end
        stats.data(leg,iJ).raw = raw; % avg joint angles by fly for ROI1 vs ROI2
        [~, stats.p_val(leg,iJ)] = ttest(raw(:,1),raw(:,2));
        raw = [];
        %labels:
        matLabel{leg,iJ} = [fly(ifly).leg_labels{leg} ' ' fly(ifly).Joints{iJ}];
    end
end
%reshape to linear (columns stack)
p_val = reshape(stats.p_val,numel(stats.p_val),1);
locLabels = reshape(matLabel,numel(matLabel), 1);

% multiple comparisons test:
fprintf('\n Significance vs multiple comparisons\n')
[p_err,p_loc] = sort(p_val);
for idx = 1:numel(p_val)
    q = (p_err((idx)) > 0.05/(length(p_err) +1 - idx));
    r = (p_err((idx)) > (idx/length(p_err))*.05);
    fprintf(['\n ' locLabels{p_loc(idx)} ' significant? ' num2str(~q) ' vs ' num2str(~r) ])
    sig(p_loc(idx)) = ~q;
end
fprintf('\n -- Stats Done -- \n')

stats.sig = reshape(sig,nLegs(end),nJoints(end));
% ----------------------------------------


% make a heatmap with diff colors for flexion vs extension
if color_opt == true
    baseColor = 'black';
else 
    baseColor = 'white';
end

n = length(nJoints)*length(nLegs);
flex = get_color('mediumvioletred', baseColor,n);
extension = get_color(baseColor, 'lime',n);
pkChange = max(abs(data),[],'all');
cMap = [flex; extension];
xlab = fly(1).leg_labels;
ylab = fly(1).Joints;
if isnan(colorLimit)
    C_limit = ceil(pkChange);
elseif pkChange>colorLimit
    warndlg('Max joint range is greater than set limits')
else 
    C_limit = colorLimit;
end
clims = [-C_limit,C_limit]; % this sets the color range!


fig = getfig('Heatmap',1);
    h = imagesc(data');
    caxis(clims);
    colormap(cMap);
    cb = colorbar;
    %label and aesthetics
    cb.FontName = 'Arial';
    cb.FontSize = 10;
    cb.Label.String = 'Joint angle change (\Delta\circ)';
    box on
    ax = gca;
    set(gca,'TickLength', [0 0])
    xticklabels(xlab)
    % adjust the axes:
    curr = string(ax.YAxis.TickLabels); % extract
    curr(1:2:end) = nan; % remove every other one
    curr(2:2:end) = ylab;
    ax.YAxis.TickLabels = curr; % set
    ax.YTickLabelRotation = 90; %rotate to vertical
    ax.FontName = 'Arial';
    %labels
    xlabel('Leg')
    ylabel('Joint')
    title({[behaviorState ': Change in joint angle during activation'];...
          [FilePath.structure_name ' conds: ' num2str(condList)]})
% ----------------------------------------
if color_opt == true
    set(fig, 'color', 'k') 
    labelColor = 'w';
    backColor = 'k';
    ax.LineWidth = 1; %change axes lines to width of 1
    set(ax, 'color', backColor, 'YColor', labelColor, 'XColor', labelColor) 
    set(ax, 'FontName', 'Arial');
    % set the label sizes and color
    labelHandles = findall(ax, 'type', 'text', 'handlevisibility', 'off');
    set(labelHandles,'FontSize', 12, 'color', labelColor)
    cb.Color = labelColor;
else
    labelColor = 'k';
    backColor = 'w';
end
% ----------------------------------------
offset = 0.25;
if stats_text == true 
   for leg = nLegs
       for iJ = nJoints
           texttag = ['p = ' num2str(round(stats.p_val(leg,iJ),4))];
           t = text(leg-offset, iJ, texttag);
           set(t,'Color', labelColor, 'FontName', 'Arial', 'FontSize',10);
           if stats.sig(leg,iJ)==true
               T = text(leg-offset, iJ-offset, '*');
               set(T,'Color', labelColor, 'FontName', 'Arial', 'FontSize',25);
           end
       end
   end
end
% ----------------------------------------

% save figure: 
save_figure(fig, [fig_dir, '\' behaviorState ' flies-joint change heatmap conds-' num2str(condList)],fileType);

clearvars('-except',initial_vars{:})


%% Average change in joint angle by fly:
% plot the avg joint angle in control vs stim periods for each fly in the group
behaviorState = 'stationary';  % behavior group to isolate e.g. 'walking', 'stationary','all'
figTag = '720 ms before and after'; % title tag
color_opt = true;      % false = white background, true = black background formatting
param = fly(1).param;   % assumes all flies have the same params (#conds,etc)
light_on = param.basler_delay*fps;
num = NUM(param);       % pull num reps/conds from a fly
paired = true;          % plots a connected line between G1&G2
fileType = '-png';      % figure save type (e.g. '-png' or '-pdf')

% define the two ROIs for comparison
G1.conds = [7,14];  % 720 ms laser condition
G2.conds = [7,14];  % 720 ms laser condition
nJoints = 1:3;          % number of joints to cycle through
nLegs = 1:6;            % legs to use

% auto will be based on the first condition listed...
switch questdlg('Auto fit ROI to G1=control vs G2=stim?')
    case 'Yes' 
        G1.ROI = 1:light_on;
        G2.ROI = light_on+1:light_on+1+round(param.conds_matrix(G2.conds(1)).opto*fps);
    case 'No'
        G1.ROI = 151:366;     % first 0.5 sec (control)
        G2.ROI = 151:366;   % 720 ms of stimulus
end


% PULL DATA:

% Group 1 data extraction
for ifly = 1:nflies
%     % skip flies without sphere fits
%     if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(G1.conds)
        cond = G1.conds(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
%             % check for possible data:
%             if isempty(step(ifly).ballFit(cond,rep).fit)...
%                || step(ifly).ballFit(cond,rep).fit==false
%                  continue;
            idx = idx+1;
            % add data for each leg into the G1 data struct
            for iJ = nJoints
              for leg = nLegs
                input = fly(ifly).angles.Leg(leg).data(cond,rep).all(G1.ROI,iJ);
                G1.f(ifly).data(leg,iJ).all(:,idx) = input;
              end
            end
            
        end
    end
end
% Group 2 data extraction
for ifly = 1:nflies
%     % skip flies without sphere fits
%     if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(G2.conds)
        cond = G2.conds(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
%             % check for possible data:
%             if isempty(step(ifly).ballFit(cond,rep).fit)...
%                || step(ifly).ballFit(cond,rep).fit==false
%                  continue;
            idx = idx+1;
            % add data for each leg into the G1 data struct
            for iJ = nJoints
              for leg = nLegs
                input = fly(ifly).angles.Leg(leg).data(cond,rep).all(G2.ROI,iJ);
                G2.f(ifly).data(leg,iJ).all(:,idx) = input;
              end
            end
        end
    end
end

% average the trials with the ROI groups:
idx = 0;
for ifly = 1:length(G1.f)
    % screen out nonmatching data
    if isempty(G1.f(ifly).data); continue; end
    idx = idx+1;
    for leg = nLegs
        for iJ = nJoints
            input = G1.f(ifly).data(leg,iJ).all;
            trialMean = mean(input);
            flyMean = mean(trialMean);
            % save data to group:
            G1.data(leg,iJ).avg(idx) = flyMean;
        end
    end
end
idx = 0;
for ifly = 1:length(G2.f)
    % screen out nonmatching data
    if isempty(G2.f(ifly).data); continue; end
    idx = idx+1;
    for leg = nLegs
        for iJ = nJoints
            input = G2.f(ifly).data(leg,iJ).all;
            trialMean = mean(input);
            flyMean = mean(trialMean);
            % save data to group:
            G2.data(leg,iJ).avg(idx) = flyMean;
        end
    end
end

% PLOT A FIGURE:
%plot avg joint angle for each fly within the cond groups for each leg
%3 graphs -- one per joint
SZ = 30;
LW = 1;
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};


fig = getfig('Joint Angle Change',1);
for iJ = nJoints
    subplot(nJoints(end),1,iJ)
    hold on
    % Plotting:
    for leg = nLegs
        x1 = (leg*2-1)*ones(1,length(G1.data(leg,iJ).avg));
        scatter(x1, G1.data(leg,iJ).avg, SZ, Color(leg_colors{leg}), 'filled')
        %group 2
        x2 = (leg*2)*ones(1,length(G2.data(leg,iJ).avg));
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
%     ylim([0,180])
    ylabel(['\Delta ' fly(1).Joints{iJ} ' angle (\circ)'])
    xlabel('Legs')
    ax = gca;
    set(ax, 'XTick', ((nLegs.*2)-0.5))
    set(ax, 'XTickLabels', fly(1).leg_labels)
end
subplot(nJoints(end),1,1)
title({FilePath.structure_name; 'Change in joint angle during activation'; [behaviorState ': ' figTag]});


fig = formatFig(fig, color_opt, [nJoints(end),1]);
save_figure(fig, [fig_dir '\Change in joint angle during activation ' behaviorState '-' figTag], fileType);

clearvars('-except',initial_vars{:})


%% Timecourse figure for change in joint angles for each leg:
% Plot the avg joint angle with err over time
behaviorState = 'walking';  % behavior group to isolate e.g. 'walking', 'stationary','all'
shading = true;             % shading vs no shading
color_opt = false;          % false = white background, true = black background formatting
param = fly(1).param;       % assumes all flies have the same params (#conds,etc)
num = NUM(param);           % pull num reps/conds from a fly
fileType = '-png';          % figure save type (e.g. '-png' or '-pdf')

% define the two ROIs for comparison
condList = [7,14];          % 720 ms laser condition 
nJoints = 1:3;              % number of joints to cycle through
nLegs = 1:6;                % legs to use

% soft parameters
light_on = param.basler_delay;
laser = param.conds_matrix(condList(1)).opto;
rec = param.basler_length;
time = linspace(-light_on,(rec-light_on),(fps*rec+1));
Jcolors = {'purple', 'darkgoldenrod', 'teal'}; % colors for the joints
% G1.color = {'Black', 'Darkslategrey', 'darkgrey'};         % group 1 color choice
% G2.color = {'Indigo', 'Orangered', 'Teal'};     % group 1 color choice
LW = 1;

% Data extraction
for ifly = 1:nflies
%     % skip flies without sphere fits
%     if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(condList)
        cond = condList(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
%             % check for possible data:
%             if isempty(step(ifly).ballFit(cond,rep).fit)...
%                || step(ifly).ballFit(cond,rep).fit==false
%                  continue;
%             end
            idx = idx+1;
            % add data for each leg into the G1 data struct
            for leg = nLegs
              % pull joint angles over time data:
              input = fly(ifly).angles.Leg(leg).data(cond,rep).all;
              for iJ = nJoints
                % find the change:
                G(leg,iJ).fly(ifly).all(:,idx) = input(:,iJ);
              end
            end
        end
    end
end

% Pool data & find fly averages:
for leg = nLegs
    for iJ = nJoints
        idx = 0; [raw, changeRaw] = deal([]);
        for ifly = 1:length(G(leg,iJ).fly)
            % screen out nonmatching data
            temp = G(leg,iJ).fly(ifly).all;
            if isempty(temp); continue; end
            idx = idx+1;   
            raw(:,idx) = mean(temp,2);
        end
        % find the cross-fly average
        % control avg: 
%         offset = mean(raw(1:light_on*fps,:));
%         changeRaw = raw - offset;
        data(leg,iJ).all = raw;
        data(leg,iJ).avg = nanmean(raw,2);
        data(leg,iJ).err = nanstd(raw,0,2);
    end
end

% Plot the data for each joint angle by leg:
ncol = 2; % FOR BOTH SIDES OF THE BODY
subplot_list = [1:ncol:nLegs(end),2:ncol:nLegs(end)];
nrow = ceil(nLegs(end)/ncol);
% show joint name in appropriate color
if color_opt == true
    labelColor = 'white';
else
    labelColor = 'black';
end
Joints = ['\color{' labelColor '}Joints:'];
for iJ = nJoints
    % format: '\color[rbg]{0 0 0}TEXT'
    a = sprintf('%s{%f %f %f}%s', '\color[rgb]', Color(Jcolors{iJ}),fly(1).Joints{iJ});
    Joints = [Joints ' ' a];
end


fig = getfig('Change in joint angle over time', 1);
for leg = nLegs
    subplot(nrow,ncol,subplot_list(leg)); hold on
    
    for iJ = nJoints
        y = data(leg,iJ).avg;
        err = data(leg,iJ).err;
        if shading==true
            fill_data = error_fill(time, y, err);
            h = fill(fill_data.X, fill_data.Y, Color(Jcolors{iJ}), 'EdgeColor','none');
            set(h, 'facealpha', 0.4)
        end
        plot(time,y+err, 'linewidth', 0.5, 'color', Color(Jcolors{iJ}))
        plot(time,y-err, 'linewidth', 0.5, 'color', Color(Jcolors{iJ}))
        plot(time,y, 'linewidth', LW, 'color', Color(Jcolors{iJ}))
    end
    ylim([0,180])
    ax = gca;
    ax.YTick = [0 90 180];
    % Add the laser indicator
    y1 = rangeLine(fig);
    plot([0,laser],[y1,y1], 'linewidth', 3, 'color', Color(param.LED_light))
    v_line(0, labelColor, ':')
    title(['Leg ' fly(1).leg_labels{leg}])
    % ylabels on left-most column only
    if ~(mod(subplot_list(leg),2)==0) % odd number subplot->label yaxis
        ylabel('Joint Angle (\circ)')
    end
    % xlabels on bottom row only
    if subplot_list(leg)>(nLegs(end)-ncol)
        xlabel('Time (s)')
    end
end

fig = formatFig(fig,color_opt,[nrow,ncol]);

% add joint labels to the first graph:
subplot(nrow,ncol,subplot_list(1));
toptitle = ['\color{' labelColor '}Leg ' fly(1).leg_labels{1}];
title({toptitle; Joints})

% save figure: 
save_figure(fig, [fig_dir, '\' behaviorState ' flies joint angles ' num2str(laser) ' opto'],fileType);

clearvars('-except',initial_vars{:})


%% Histogram of joint angles
% plots histogram of all frames within two ROIs for 2 groups
% plots the fly joint angle mean-of-means above histo with significance star
% ---- input -----
iJ = 2;                 % joint selection
behaviorState = 'walking';  % behavior group to isolate e.g. 'walking', 'stationary','all'
% figTag = '720 ms before and after'; % title tag %TODO
color_opt = false;      % false = white background, true = black background formatting
cleanAxes = true;       % set all axes to same limits & standarize tick marks
fileType = '-png';      % figure save type (e.g. '-png' or '-pdf')
G1.conds = [7,14];      % group 1 conditions to lump together, can be as many as desired
G2.conds = [7,14];      % group 2 conditions to lump together
G1.color = 'SlateGrey';     % group 1 color choice
G2.color = 'blue'; % group 1 color choice
nLegs = 1:6;            % leg selection
fileType = '-png';      % figure save type (e.g. '-png' or '-pdf')
% ----------------
param = fly(1).param;   % assumes all flies have the same params (#conds,etc)
light_on = param.basler_delay*fps;
num = NUM(param);       % pull num reps/conds from a fly
% auto will be based on the first condition listed for each group
switch questdlg('Auto fit ROI to G1=control vs G2=stim?')
    case 'Yes' 
        G1.ROI = 1:light_on;
        G2.ROI = light_on+1:light_on+1+round(param.conds_matrix(G2.conds(1)).opto*fps);
    case 'No'
        G1.ROI = 1:light_on;    % control ROI *current
        G2.ROI = light_on:(light_on+(0.5*fps))-1; %0.5sec starting at stim on
end


% Group 1 data extraction
for ifly = 1:nflies
%     % skip flies without sphere fits
%     if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(G1.conds)
        cond = G1.conds(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep), behaviorState)
              continue;
            end
%             % check for possible data:
%             if isempty(step(ifly).ballFit(cond,rep).fit)...
%                || step(ifly).ballFit(cond,rep).fit==false
%                  continue;
            idx = idx+1;
            % add data for each leg into the G1 data struct
            for leg = nLegs
              input = fly(ifly).angles.Leg(leg).data(cond,rep).all(G1.ROI,iJ);
              G1.f(ifly).data(leg).all(:,idx) = input;
            end
        end
    end
end
% Group 2 data extraction
for ifly = 1:nflies
%     % skip flies without sphere fits
%     if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(G2.conds)
        cond = G2.conds(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
%             % check for possible data:
%             if isempty(step(ifly).ballFit(cond,rep).fit)...
%                || step(ifly).ballFit(cond,rep).fit==false
%                  continue;
            idx = idx+1;
            % add data for each leg into the G1 data struct
            for leg = nLegs
              input = fly(ifly).angles.Leg(leg).data(cond,rep).all(G2.ROI,iJ);
              G2.f(ifly).data(leg).all(:,idx) = input;
            end
        end
    end
end

% average the trials with the ROI groups:
idx = 0;
for ifly = 1:length(G1.f)
    % screen out nonmatching data
    if isempty(G1.f(ifly).data); continue; end
    idx = idx+1;
    for leg = nLegs
        input = G1.f(ifly).data(leg).all;
        trialMean = mean(input);
        flyMean = mean(trialMean);
        % save data to group:
        G1.data(leg).avg(idx) = flyMean;
        if idx==1
            G1.data(leg).all = input;
        else
            loc = size(G1.data(leg).all,2);
            width = size(input,2);
            G1.data(leg).all(:,loc+1:loc+width) = input;
        end
    end
end
idx = 0;
for ifly = 1:length(G1.f)
    % screen out nonmatching data
    if isempty(G2.f(ifly).data); continue; end
    idx = idx+1;
    for leg = nLegs
        input = G2.f(ifly).data(leg).all;
        trialMean = mean(input);
        flyMean = mean(trialMean);
        % save data to group:
        G2.data(leg).avg(idx) = flyMean;
        if idx==1
            G2.data(leg).all = input;
        else
            loc = size(G2.data(leg).all,2);
            width = size(input,2);
            G2.data(leg).all(:,loc+1:loc+width) = input;
        end
    end
end

% Check for significance in the groups:
% show the p-value of a t-test between G1&G2 --  for paired group (aka same
% conds)
[~,p_val] = deal([]);
for leg = nLegs
    try
        [~, p_val(leg)] = ttest(G1.data(leg).avg,G2.data(leg).avg);
    catch
        p_val(leg) = nan;
    end
end
% multiple comparisons test:
fprintf('\n Significance vs multiple comparisons\n')
[p_err,p_loc] = sort(p_val);
for idx = nLegs
    q = (p_err((idx)) > 0.05/(length(p_err) +1 - idx));
    r = (p_err((idx)) > (idx/length(p_err))*.05);
    fprintf(['\nLeg ' fly(1).leg_labels{p_loc(idx)} ' significant? ' num2str(~q) ' vs ' num2str(~r) ])
    sig(p_loc(idx)) = ~q;
end
fprintf('\nStats Done\n')


% -------------PLOTS-------------

subplot_loc = [1 3 5 2 4 6]; % organize left legs on leg and right on right
Edges = 0:10:180; % edges for histogram bins
LW = 1;
SZ = 100;
err_tranceparency = 0.6; %transparency on the error plot lines

Joints = fly(1).Joints;
leg_labels = fly(1).leg_labels;
if color_opt==true
    baseColor = 'white';
    if strcmpi(G1.color,'black')
        G1.color = baseColor;
    end
    if strcmpi(G2.color,'black')
        G2.color = baseColor;
    end 
else
    baseColor = 'black';
end


% Make joint angle histograms:
fig = getfig('Joint Angle Histogram',1);
for leg = 1:6   % subplot for each leg
    subplot(3,2,subplot_loc(leg))
    hold on
    % plot group 1:
    h = histogram(G1.data(leg).all,Edges);
    h.FaceColor = Color(G1.color); 
    % plot group 2:
    h = histogram(G2.data(leg).all,Edges);
    h.FaceColor = Color(G2.color); 
    
    % sort tick marks to make space for average bars
    ax = gca;
    yL = ylim;
    y = yL(2)*1.10;
    
    % G1 avg and err line:
    avg1 = mean(G1.data(leg).avg);
    err = std(reshape(G1.data(leg).all,numel(G1.data(leg).all),1));
    scatter(avg1,y, SZ, Color(G1.color), 'filled','s')
    P = plot([avg1-err,avg1+err],[y,y], 'linewidth', LW, 'color', Color(G1.color));
    P.Color(4) = err_tranceparency;
    % G2 avg and err line:
    avg2 = mean(G2.data(leg).avg);
    err = std(reshape(G2.data(leg).all,numel(G2.data(leg).all),1));
    scatter(mean(G2.data(leg).avg),y, SZ, Color(G2.color), 'filled','s')
    p = plot([avg2-err,avg2+err],[y,y], 'linewidth', LW, 'color', Color(G2.color));
    p.Color(4) = err_tranceparency;
    
    if sig(leg)==true
        y = yL(2)*1.22;
        x = avg1-((avg1-avg2)/2);
        t = scatter(x,y,50,'*');
        t.MarkerFaceColor = 'none';
        t.MarkerEdgeColor =  Color(baseColor); 
%         t.LineWidth = 1;
    end
    
    % peep the ylimits:
    finalYlim(leg,:) = ylim;
    
    % plot labels:
    if subplot_loc(leg)>4
        xlabel('Angle (\circ)')
    end
    ylabel([Joints{iJ} ' (count)'])
    title(['Leg ' leg_labels{leg} ' | p=' num2str(p_val(leg))]); 
    set(gca, 'TickDir', 'out')
    xlim([0,180])
end

% tidy axes:
if cleanAxes==true
    ymax = max(finalYlim(:,2));
    K = 1:ceil(ymax/3);
    D = [K(rem(ymax,K)==0) ymax];
    [~,idx] = min(abs(D-3));
    new_yticks = linspace(0,ymax,D(idx));
    for leg = nLegs
        subplot(3,2,leg)
        ax = gca;
        ylim([0,ymax])
        ax.YTick = new_yticks;
        ax.XTick = [0,90,180];
    end
end

% save options:
fig = formatFig(fig, color_opt, [3,2]);
save_figure(fig, [fig_dir '\Histogram of ' Joints{iJ} ' - ' behaviorState], fileType);
 
clearvars('-except',initial_vars{:})


%% Step Frequency over time by leg (currently only for walking since no other sphere fits happen)

behaviorState = 'walking';  % behavior group to isolate e.g. 'walking', 'stationary','all'
shading = true;            % shading vs no shading
allTrials = false;           % plot lines for each trial, not just the avg
color_opt = false;          % false = white background, true = black background formatting
norm_axes = true;           % standardize the axes across legs?
param = fly(1).param;       % assumes all flies have the same params (#conds,etc)
num = NUM(param);           % pull num reps/conds from a fly
fileType = '-png';          % figure save type (e.g. '-png' or '-pdf')

% define the two ROIs for comparison
condList = [7,14];          % 720 ms laser condition 
nLegs = 1:6;                % legs to use

% soft parameters
light_on = param.basler_delay;
laser = param.conds_matrix(condList(1)).opto;
rec = param.basler_length;
time = linspace(-light_on,(rec-light_on),(fps*rec+1));
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
LW = 1;


% Calculate step frequency
for ifly = 1:nflies
    % skip flies without sphere fits
    if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(condList)
        cond = condList(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
            
            % check for poor fit data:
            if isempty(step(ifly).ballFit(cond,rep).fit)...
               || step(ifly).ballFit(cond,rep).fit==false
                 continue;
            end
            
            % Print stats about each trial that is selected: (error catching)
            state = fly(ifly).behavior.behavior{cond,rep};
            id = ['fly-' num2str(ifly) ', cond-' num2str(cond) ', rep-' num2str(rep) ' : '];
            disp([id state])
            
            idx = idx+1;
            % pull the logical of swing|stance per frame
            raw = step(ifly).stance(cond,rep).loc'; %flip orientation
            freq = [];
            for leg = nLegs %each leg
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
            gait(ifly).Freq(idx).data = freq;
            gait(ifly).Freq(idx).f_all = F_all;
            gait(ifly).Freq(idx).f_avg = mean(F_all,2);
        end
    end
end

% Pool step frequency data across flies
% note: you should ignore the first bit as there's an edge effect in
% predicting the frequency            
idx = 0; freq = [];
for ifly = 1:nflies
    ntrials = length(gait(ifly).Freq);
    for ii = 1:ntrials
      idx = idx+1;  
      for leg = nLegs
        f = gait(ifly).Freq(ii).f_all(:,leg);
        if isempty(f); continue; end
        freq(leg).data(:,idx) = f;
      end
    end
end

% Find the averages
for leg = nLegs
    freq(leg).avg = mean(freq(leg).data,2);
    freq(leg).err = sem(freq(leg).data,2);
end

% Plot all traces overlaid:
xlimits = [-0.4,1.5];

nrow = 3;
ncol = 2;
sbpt = [1:2:6,2:2:6];
fig = getfig('Step Frequency',1);
for leg = nLegs
    subplot(nrow,ncol,sbpt(leg))
    hold on
    y = freq(leg).avg;
    err = freq(leg).err;
    % plot individual traces
    if allTrials==true
       plot(time, freq(leg).data, 'color', Color('grey'), 'linewidth', 0.5)
    end
    
    % plot error region
    if shading==true
        fill_data = error_fill(time, y, err);
        h = fill(fill_data.X, fill_data.Y, Color(leg_colors{leg}), 'EdgeColor','none');
        set(h, 'facealpha', 0.4)
    end
    % plot error edges and avg line *not optional
    plot(time, y, 'color', Color(leg_colors{leg}), 'linewidth', LW)
    plot(time, y+err, 'color', Color(leg_colors{leg}), 'linewidth', 0.5)
    plot(time, y-err, 'color', Color(leg_colors{leg}), 'linewidth', 0.5)
    
    % labels etc
    xlabel('time (s)')
    ylabel('Step frequency (hz)')
    % laser region lines:
    xlim(xlimits)
    if norm_axes==false
        y1 = rangeLine(fig);
        plot([0,light_on], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
    end
    %pull ylimits to later normalize
    ylimits(leg,:) = ylim;

end

if norm_axes==true
    ub = max(ylimits(:,2));
    lb = min(ylimits(:,1));
    for leg = nLegs
        subplot(nrow,ncol,leg)
        ylim([lb,ub])
        % plot stim region
        y1 = rangeLine(fig,2);
        plot([0,light_on], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
    end
end

fig = formatFig(fig,color_opt,[nrow,ncol]); 
 
save_figure(fig, [fig_dir, '\' behaviorState ' step frequency for conds ' num2str(condList)], fileType);


clearvars('-except',initial_vars{:})


%% Average step frequency (legs combined) compare two groups

behaviorState = 'walking';  % behavior group to isolate e.g. 'walking', 'stationary','all'
shading = true;            % shading vs no shading
allTrials = false;           % plot lines for each trial, not just the avg
color_opt = true;          % false = white background, true = black background formatting
norm_axes = true;           % smae axes range for both time course figures
paired = true;              % line between control-stim points comparison
param = fly(1).param;       % assumes all flies have the same params (#conds,etc)
num = NUM(param);           % pull num reps/conds from a fly
fileType = '-png';          % figure save type (e.g. '-png' or '-pdf')

% define the two ROIs for comparison
G1.conds = [1,8];           % control 0ms groups
G2.conds = [5,12];          % laser 720 ms condition 
nLegs = 1:6;                % legs to use
G1.color = 'black';         % control color
G2.color = 'teal';          % 'stim' group color
if color_opt==true
    if strcmpi(G1.color, 'black')
        G1.color = 'white';
    end
    if strcmpi(G2.color, 'black')
        G2.color = 'white';
    end
end

% soft parameters
light_on = param.basler_delay;
G1.laser = param.conds_matrix(G1.conds(1)).opto;
G2.laser = param.conds_matrix(G2.conds(1)).opto;
G1.light_off = param.conds_matrix(G1.conds(1)).opto;
G2.light_off = param.conds_matrix(G2.conds(1)).opto;
rec = param.basler_length;
time = linspace(-light_on,(rec-light_on),(fps*rec+1));
LW = 1;
SZ = 50;

% auto will be based on the first condition listed...
switch questdlg({'Auto fit G1 stim ROI laser length from G2?';...
                 '(If G1 isn''t the 0 ms condition, then no)'})
    case 'Yes' 
        control_ROI = 30:round(light_on*fps); %strt 30 to avoid edge effect
        stim_ROI = round(light_on*fps)+1:round(light_on*fps)+1+round(param.conds_matrix(G2.conds(1)).opto*fps);
    case 'No' %manually select prestim and stim regions
        control_ROI = 1:150;     % first 0.5 sec (control)
        stim_ROI = 151:366;   % 720 ms of stimulus
end

% Graph parameters
xlimits = [-0.4,1.6];
nrow = 2;
ncol = 3;
G1sp = [1,2];   %subplot idx for G1 timecourse
G2sp = [4,5];   %subplot idx for G2 timecourse
scsb = [3,6];   %subplot idx for change in freq scatter
plotIdx(1).idx = G1sp;
plotIdx(2).idx = G2sp;
plotIdx(3).idx = scsb;

% G1: Calculate step frequency
disp('G1 selected trials')
for ifly = 1:nflies
    % skip flies without sphere fits
    if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(G1.conds)
        cond = G1.conds(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
            
            % check for poor fit data:
            if isempty(step(ifly).ballFit(cond,rep).fit)...
               || step(ifly).ballFit(cond,rep).fit==false
                 continue;
            end
            
            % Print stats about each trial that is selected: (error catching)
            state = fly(ifly).behavior.behavior{cond,rep};
            id = ['fly-' num2str(ifly) ', cond-' num2str(cond) ', rep-' num2str(rep) ' : '];
            disp([id state])
            
            idx = idx+1;
            % pull the logical of swing|stance per frame
            raw = step(ifly).stance(cond,rep).loc'; %flip orientation
            freq = [];
            for leg = nLegs %each leg
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
            G1.gait(ifly).Freq(idx).data = freq;
            G1.gait(ifly).Freq(idx).f_all = F_all;
            G1.gait(ifly).Freq(idx).f_avg = mean(F_all,2);
        end
    end
end

% G1: Calculate step frequency
disp('G2 selected trials')
for ifly = 1:nflies
    % skip flies without sphere fits
    if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(G2.conds)
        cond = G2.conds(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
            
            % check for poor fit data:
            if isempty(step(ifly).ballFit(cond,rep).fit)...
               || step(ifly).ballFit(cond,rep).fit==false
                 continue;
            end
            
            % Print stats about each trial that is selected: (error catching)
            state = fly(ifly).behavior.behavior{cond,rep};
            id = ['fly-' num2str(ifly) ', cond-' num2str(cond) ', rep-' num2str(rep) ' : '];
            disp([id state])
            
            idx = idx+1;
            % pull the logical of swing|stance per frame
            raw = step(ifly).stance(cond,rep).loc'; %flip orientation
            freq = [];
            for leg = nLegs %each leg
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
            G2.gait(ifly).Freq(idx).data = freq;
            G2.gait(ifly).Freq(idx).f_all = F_all;
            G2.gait(ifly).Freq(idx).f_avg = mean(F_all,2);
        end
    end
end

% G1: Pool step frequency data across flies           
[flycount,idx] = deal(0); G1.freq = [];
for ifly = 1:nflies
    ntrials = length(G1.gait(ifly).Freq);
    if isempty(ntrials); continue; end
    flycount = flycount+1;
    for ii = 1:ntrials
      idx = idx+1;  
      for leg = nLegs
        f = G1.gait(ifly).Freq(ii).f_all(:,leg);
        G1.freq(leg).data(:,idx) = f;
      end
    end
end
disp(['G1 N-flies : ' num2str(flycount)])
disp(['G1 N-trials : ' num2str(idx)])
G1.ntrials = idx;

% G2: Pool step frequency data across flies           
[flycount,idx] = deal(0); G2.freq = [];
for ifly = 1:nflies
    ntrials = length(G2.gait(ifly).Freq);
    if isempty(ntrials); continue; end
    flycount = flycount+1;
    for ii = 1:ntrials
      idx = idx+1;  
      for leg = nLegs
        f = G2.gait(ifly).Freq(ii).f_all(:,leg);
        G2.freq(leg).data(:,idx) = f;
      end
    end
end
disp(['G2 N-flies : ' num2str(flycount)])
disp(['G2 N-trials : ' num2str(idx)])
G2.ntrials = idx;


% G1: Find the averages (avg legs, then across flies)
temp = [];
for n = 1:G1.ntrials
  for leg = nLegs
    temp(:,leg) = G1.freq(leg).data(:,n);
  end
  G1.trialmean(:,n) = mean(temp,2);
end
G1.avg = mean(G1.trialmean,2);
G1.err = std(G1.trialmean,0,2);
%pull control vs stim avg for scatter plot
G1.control = mean(G1.trialmean(control_ROI,:));
G1.stim = mean(G1.trialmean(stim_ROI,:));

% G2: Find the averages (avg legs, then across flies)
temp = [];
for n = 1:G2.ntrials
  for leg = nLegs
    temp(:,leg) = G2.freq(leg).data(:,n);
  end
  G2.trialmean(:,n) = mean(temp,2);
end
G2.avg = mean(G2.trialmean,2);
G2.err = std(G2.trialmean,0,2);
%pull control vs stim avg for scatter plot
G2.control = mean(G2.trialmean(control_ROI,:));
G2.stim = mean(G2.trialmean(stim_ROI,:));


% PLOT HERE:
fig = getfig('Step Frequency two groups',1);
% G1 time course
subplot(nrow,ncol,G1sp)
    hold on
    y = G1.avg;
    err = G1.err;
    % plot individual traces
    if allTrials==true
       plot(time, G1.trialmean, 'color', Color('grey'), 'linewidth', 0.5)
    end

    % plot error region
    if shading==true
        fill_data = error_fill(time, y, err);
        h = fill(fill_data.X, fill_data.Y, Color(G1.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.4)
    end
    % plot error edges and avg line *not optional
    plot(time, y, 'color', Color(G1.color), 'linewidth', LW)
    plot(time, y+err, 'color', Color(G1.color), 'linewidth', 0.5)
    plot(time, y-err, 'color', Color(G1.color), 'linewidth', 0.5)

    % labels etc
    xlabel('time (s)')
    ylabel('Step frequency (hz)')
    title(['Group 1 : conds ' num2str(G1.conds)])
    % laser region lines:
    xlim(xlimits)
    if norm_axes==false
        y1 = rangeLine(fig);
        plot([0,G1.light_off], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
    end
    %pull ylimits to later normalize
    ylimits(1,:) = ylim;

% G2 time course
subplot(nrow,ncol,G2sp)
    hold on
    y = G2.avg;
    err = G2.err;
    % plot individual traces
    if allTrials==true
       plot(time, G2.trialmean, 'color', Color('grey'), 'linewidth', 0.5)
    end

    % plot error region
    if shading==true
        fill_data = error_fill(time, y, err);
        h = fill(fill_data.X, fill_data.Y, Color(G2.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.4)
    end
    % plot error edges and avg line *not optional
    plot(time, y, 'color', Color(G2.color), 'linewidth', LW)
    plot(time, y+err, 'color', Color(G2.color), 'linewidth', 0.5)
    plot(time, y-err, 'color', Color(G2.color), 'linewidth', 0.5)

    % labels etc
    xlabel('time (s)')
    ylabel('Step frequency (hz)')
    title(['Group 2 : conds ' num2str(G2.conds)])
    % laser region lines:
    xlim(xlimits)
    if norm_axes==false
        y1 = rangeLine(fig);
        plot([0,G2.light_off], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
    end
    %pull ylimits to later normalize
    ylimits(2,:) = ylim;
    
%standardize axes for time course figures if opted
    if norm_axes==true
        ub = max(ylimits(:,2));
        lb = min(ylimits(:,1));
        subplot(nrow,ncol,[1,2])
            ylim([lb,ub])
            % plot stim region
            y1 = rangeLine(fig,2);
            plot([0,G1.light_off], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
        subplot(nrow,ncol,[4,5])
            ylim([lb,ub])
            % plot stim region
            y1 = rangeLine(fig,2);
            plot([0,G2.light_off], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
    end

% control vs stim ROI both groups:
subplot(nrow,ncol,scsb); hold on
    % G1
    x1 = 1*ones(1,G1.ntrials);
    x2 = 2*ones(1,G1.ntrials);
    scatter(x1, G1.control, SZ, Color(G1.color), 'filled')
    scatter(x2, G1.stim, SZ, Color(G1.color), 'filled')
    if paired==true
      for n = 1:G1.ntrials
        plot([x1(n),x2(n)], [G1.control(n), G1.stim(n)], ...
            'color', Color(G1.color), 'linewidth', LW)
      end
    end
    % G2
    x1 = 3*ones(1,G2.ntrials);
    x2 = 4*ones(1,G2.ntrials);
    scatter(x1, G2.control, SZ, Color(G2.color), 'filled')
    scatter(x2, G2.stim, SZ, Color(G2.color), 'filled')
    if paired==true
      for n = 1:G2.ntrials
        plot([x1(n),x2(n)], [G2.control(n), G2.stim(n)], ...
            'color', Color(G2.color), 'linewidth', LW)
      end
    end
    % labels etc
    xlim([0,5])
    title('Control vs stim')
    ylabel('Mean step frequency (Hz)')
    xlabel('Groups')
    ax = gca;
    set(ax, 'XTick', [1.5,3.5])
    set(ax, 'XTickLabels', {'G1', 'G2'})
    title('Step freq control vs. stim')
   
fig = formatFig(fig,color_opt,[nrow,ncol],plotIdx); 

save_figure(fig, [fig_dir, '\' behaviorState ' step frequency between groups'], fileType);

clearvars('-except',initial_vars{:})


%% Stride distance over time for each leg %% LOOK HERE
euclidian_dist = true;     % euclidian or total step distance
behaviorState = 'walking';  % behavior group to isolate e.g. 'walking', 'stationary','all'
shading = true;            % shading vs no shading
allTrials = false;           % plot lines for each trial, not just the avg
color_opt = true;          % false = white background, true = black background formatting
norm_axes = true;           % smae axes range for both time course figures
paired = true;              % line between control-stim points comparison
print_trials = true;        % print the ID of all trials being selected *good for errors
param = fly(1).param;       % assumes all flies have the same params (#conds,etc)
num = NUM(param);           % pull num reps/conds from a fly
fileType = '-png';          % figure save type (e.g. '-png' or '-pdf')
LP = 5;                     % leg position for movement trace (e.g. tarsus stride dist)

% define the two ROIs for comparison
G1.conds = [7,14];           % control 0ms groups
nLegs = 1:6;                % legs to use
% 

if color_opt==true
    baseColor = 'white';
    l_color = 'black';
else 
    baseColor = 'black';
    l_color = 'white';
end

% soft parameters
light_on = param.basler_delay;
G1.laser = param.conds_matrix(G1.conds(1)).opto;
G1.light_off = param.conds_matrix(G1.conds(1)).opto;
rec = param.basler_length;
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
leg_labels = fly(1).leg_labels;

% Graph parameters
time = linspace(-light_on,(rec-light_on),(fps*rec+1));
xlimits = [-0.4,1.6];   % cuts of beginning of timecouse (edge effect)
nrow = 3;               % number of rows (3 legs per side)
ncol = 2;               % number of colums (L & R)
LW = 1;                 % linewidth
SZ = 50;                % scatter point size
sbpt = [1,3,5,2,4,6];   % subplot locations for the legs

% % auto will be based on the first condition listed...
% switch questdlg('Autofit ROIs to control vs laser on?')
%     case 'Yes' 
%         control_ROI = 30:round(light_on*fps); %strt 30 to avoid edge effect
%         stim_ROI = round(light_on*fps)+1:round(light_on*fps)+1+round(param.conds_matrix(G1.conds(1)).opto*fps);
%     case 'No' %manually select prestim and stim regions
%         control_ROI = 1:150;     % first 0.5 sec (control)
%         stim_ROI = 151:366;   % 720 ms of stimulus
% end

% Calculate step frequency
disp('Selected trials')
for ifly = 1:nflies
    % skip flies without sphere fits
    if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(G1.conds)
        cond = G1.conds(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
            
            % check for poor fit data:
            if isempty(step(ifly).ballFit(cond,rep).fit)...
               || step(ifly).ballFit(cond,rep).fit==false
                 continue;
            end
            
            % Print stats about each trial that is selected: (error catching)
            if print_trials==true
                state = fly(ifly).behavior.behavior{cond,rep};
                id = ['fly-' num2str(ifly) ', cond-' num2str(cond) ', rep-' num2str(rep) ' : '];
                disp([id state])
            end
            
            idx = idx+1;
            % pull the logical of swing|stance per frame
            raw = step(ifly).stance(cond,rep).loc'; %flip orientation    
            
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
                data = fly(ifly).pose_3d.Leg(LP,leg).data{cond,rep}.all;
                if euclidian_dist==true
                    % find the tarsus position for each stride:
                    stride(leg).start_pos = data(stride(leg).loc(:,1),:);
                    stride(leg).end_pos = data(stride(leg).loc(:,2),:);
                    % find the euclidian distance between the start and end of the stride
                    for ii = 1:nstrides
                        stride(leg).length(ii,1) = pdist2(stride(leg).start_pos(ii,:),stride(leg).end_pos(ii,:));
                    end
                else % total distance traveled:
                    stride(leg).start_pos = data(stride(leg).loc(:,1),:);
                    stride(leg).end_pos = data(stride(leg).loc(:,2),:);
                    for ii = 1:nstrides
                        frames = stride(leg).loc(ii,1):stride(leg).loc(ii,2);
                        step_pos = data(frames,:);
                        temp = pdist2(step_pos,step_pos);
                        for tt = 1:length(frames)-1
                            d(tt) = temp(tt,tt+1);
                        end
                        stride(leg).length(ii,1) = sum(d);
                    end
                end
                % extrapolate and fill the graps between stance (i.e. swing)
                [X,Y,y] = deal([]);
                Y = nan(param.basler_length*fps+1,1);
                X = [1; stride(leg).loc(:,1); param.basler_length*fps+1];
                for ii = 1:length(X)-1
                    f = stride(leg).length;
                    if ii == length(X)-1; val = f(ii-1); 
                    else val = f(ii); end
                    strt = X(ii);
                    fin = X(ii+1);
                    Y(strt:fin-1) = val*ones(fin-strt,1);
                end
                Y(end) = Y(end-1);  
                y = smooth(Y,fps/10,'lowess'); 
                stride(leg).filled = y;
            end
            % find the average stride leg
            % assign the frequency data to the Gait structure:
            G1.gait(ifly).Stride(idx).data = stride;
        end
    end
end

% Pool step frequency data across flies           
idx = 0; G1.dist = [];
for ifly = 1:nflies
    ntrials = length(G1.gait(ifly).Stride);
    if isempty(ntrials); continue; end
    for ii = 1:ntrials
      idx = idx+1;  
      for leg = nLegs
        f = G1.gait(ifly).Stride(ii).data(leg).filled;
        G1.all(leg,ifly).data(:,idx) = f;
        G1.dist(leg).data(:,idx) = f;
      end
    end
end

% Find averages
for leg = nLegs
    for ifly = 1:nflies
        if isempty(G1.all); continue; end
        %avg stride length for each fly
        G1.flyav(leg).data(:,ifly) = mean(G1.all(leg,ifly).data,2);
%         G1.flyav(leg).control(ifly) = mean(G1.flyav(leg).data(control_ROI,ifly));
%         G1.flyav(leg).stim(ifly) = mean(G1.flyav(leg).data(stim_ROI,ifly));
    end
    G1.avg(:,leg) = mean(G1.flyav(leg).data,2);
    G1.err(:,leg) = sem(G1.flyav(leg).data,2);
%     % t-test for each leg:
%     [~,G1.p(leg)] = ttest(G1.flyav(leg).control,G1.flyav(leg).stim);
end

% ------ PLOT HERE -------:
fig = getfig('Step distance across legs',1);
for leg = nLegs
% G1 time course
subplot(nrow,ncol,sbpt(leg))
    hold on
    y = G1.avg(:,leg);
    err = G1.err(:,leg);
    % plot individual traces
    if allTrials==true
       plot(time, G1.flyav(leg).data, 'color', Color('grey'), 'linewidth', 0.5)
    end

    % plot error region
    if shading==true
        fill_data = error_fill(time, y, err);
        h = fill(fill_data.X, fill_data.Y, Color(leg_colors{leg}), 'EdgeColor','none');
        set(h, 'facealpha', 0.4)
    end
    % plot error edges and avg line *not optional
    plot(time, y, 'color', Color(leg_colors{leg}), 'linewidth', LW)
    plot(time, y+err, 'color', Color(leg_colors{leg}), 'linewidth', 0.5)
    plot(time, y-err, 'color', Color(leg_colors{leg}), 'linewidth', 0.5)

    % labels etc
    xlabel('time (s)')
    ylabel('Step distance (au)')
    title(['Leg ' leg_labels{leg}])
    % laser region lines:
    xlim(xlimits)
    if norm_axes==false
        y1 = rangeLine(fig);
        plot([0,G1.light_off], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
    end
    %pull ylimits to later normalize
    ylimits(leg,:) = ylim;
end

%standardize axes for time course figures if opted
if norm_axes==true
    ub = max(ylimits(:,2));
    lb = min(ylimits(:,1));
    for leg = nLegs
    subplot(nrow,ncol,leg)
        ylim([lb,ub])
        % plot stim region
        y1 = rangeLine(fig,2);
        plot([0,G1.light_off], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
        l = legend({'SEM', 'avg'}); set(l,'Color',Color(l_color));
        [l.EdgeColor,l.TextColor] = deal(Color(baseColor));
    end
end

% Other labels:
subplot(nrow,ncol,1)
title({['Conds ' num2str(G1.conds)]; ['Leg ' num2str(leg_labels{1})]})
   
fig = formatFig(fig,color_opt,[nrow,ncol]); 

save_figure(fig, [fig_dir, '\' behaviorState ' step distance by leg-conds ' num2str(G1.conds)], fileType);

clearvars('-except',initial_vars{:})


%% Change in stride distance for each leg
euclidian_dist = true;     % euclidian or total step distance
behaviorState = 'walking';  % behavior group to isolate e.g. 'walking', 'stationary','all'
shading = true;            % shading vs no shading
allTrials = false;           % plot lines for each trial, not just the avg
color_opt = true;          % false = white background, true = black background formatting
norm_axes = true;           % same axes range for both time course figures
paired = true;              % line between control-stim points comparison
print_trials = true;        % print the ID of all trials being selected *good for errors
param = fly(1).param;       % assumes all flies have the same params (#conds,etc)
num = NUM(param);           % pull num reps/conds from a fly
fileType = '-png';          % figure save type (e.g. '-png' or '-pdf')
LP = 5;                     % leg position for movement trace (e.g. tarsus stride dist)

% define the two ROIs for comparison
G1.conds = [7,14];           % control 0ms groups
nLegs = 1:6;                % legs to use
if color_opt==true
    baseColor = 'white';
    l_color = 'black';
else 
    baseColor = 'black';
    l_color = 'white';
end

% soft parameters
light_on = param.basler_delay;
G1.laser = param.conds_matrix(G1.conds(1)).opto;
G1.light_off = param.conds_matrix(G1.conds(1)).opto;
rec = param.basler_length;
leg_colors = {'blue', 'red', 'orange', 'purple', 'green', 'cyan'};
leg_labels = fly(1).leg_labels;

% auto will be based on the first condition listed...
switch questdlg({'Auto fit control ROI to pre-laser and stim';...
                 ' to laser duration? (0ms defaults to 0.5s)'})
    case 'Yes' 
        control_ROI = 30:round(light_on*fps); %strt 30 to avoid edge effect
        stim_ROI = round(light_on*fps)+1:round(light_on*fps)+1+round(param.conds_matrix(G1.conds(1)).opto*fps);
        if G1.light_off==0
            stim_ROI = round(light_on*fps)+1:round(light_on*fps)+1+round(0.5*fps);
        end
    case 'No' %manually select prestim and stim regions
        control_ROI = 1:150;     % first 0.5 sec (control)
        stim_ROI = 151:366;   % 720 ms of stimulus
end


% Graph parameters
% time = linspace(-light_on,(rec-light_on),(fps*rec+1));
% xlimits = [-0.4,1.6];   % cuts of beginning of timecouse (edge effect)
% nrow = 3;               % number of rows (3 legs per side)
% ncol = 2;               % number of colums (L & R)
LW = 1;                 % linewidth
SZ = 50;                % scatter point size
% sbpt = [1,3,5,2,4,6];   % subplot locations for the legs

% Calculate step frequency
disp('Selected trials')
for ifly = 1:nflies
    % skip flies without sphere fits
    if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(G1.conds)
        cond = G1.conds(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
            
            % check for poor fit data:
            if isempty(step(ifly).ballFit(cond,rep).fit)...
               || step(ifly).ballFit(cond,rep).fit==false
                 continue;
            end
            
            % Print stats about each trial that is selected: (error catching)
            if print_trials==true
                state = fly(ifly).behavior.behavior{cond,rep};
                id = ['fly-' num2str(ifly) ', cond-' num2str(cond) ', rep-' num2str(rep) ' : '];
                disp([id state])
            end
            
            idx = idx+1;
            % pull the logical of swing|stance per frame
            raw = step(ifly).stance(cond,rep).loc'; %flip orientation    
            
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
                data = fly(ifly).pose_3d.Leg(LP,leg).data{cond,rep}.all;
                if euclidian_dist==true
                    % find the tarsus position for each stride:
                    stride(leg).start_pos = data(stride(leg).loc(:,1),:);
                    stride(leg).end_pos = data(stride(leg).loc(:,2),:);
                    % find the euclidian distance between the start and end of the stride
                    for ii = 1:nstrides
                        stride(leg).length(ii,1) = pdist2(stride(leg).start_pos(ii,:),stride(leg).end_pos(ii,:));
                    end
                else % total distance traveled:
                    stride(leg).start_pos = data(stride(leg).loc(:,1),:);
                    stride(leg).end_pos = data(stride(leg).loc(:,2),:);
                    for ii = 1:nstrides
                        frames = stride(leg).loc(ii,1):stride(leg).loc(ii,2);
                        step_pos = data(frames,:);
                        temp = pdist2(step_pos,step_pos);
                        for tt = 1:length(frames)-1
                            d(tt) = temp(tt,tt+1);
                        end
                        stride(leg).length(ii,1) = sum(d);
                    end
                end
                % extrapolate and fill the graps between stance (i.e. swing)
                [X,Y,y] = deal([]);
                Y = nan(param.basler_length*fps+1,1);
                X = [1; stride(leg).loc(:,1); param.basler_length*fps+1];
                for ii = 1:length(X)-1
                    f = stride(leg).length;
                    if ii == length(X)-1; val = f(ii-1); 
                    else val = f(ii); end
                    strt = X(ii);
                    fin = X(ii+1);
                    Y(strt:fin-1) = val*ones(fin-strt,1);
                end
                Y(end) = Y(end-1);  
                y = smooth(Y,fps/10,'lowess'); 
                stride(leg).filled = y;
            end
            % find the average stride leg
            % assign the frequency data to the Gait structure:
            G1.gait(ifly).Stride(idx).data = stride;
        end
    end
end

% Pool step frequency data across flies           
idx = 0; G1.dist = [];
for ifly = 1:nflies
    ntrials = length(G1.gait(ifly).Stride);
    if isempty(ntrials); continue; end
    for ii = 1:ntrials
      idx = idx+1;  
      for leg = nLegs
        f = G1.gait(ifly).Stride(ii).data(leg).filled;
        G1.all(leg,ifly).data(:,idx) = f;
        G1.dist(leg).data(:,idx) = f;
      end
    end
end

% Find averages
for leg = nLegs
    for ifly = 1:nflies
        if isempty(G1.all); continue; end
        %avg stride length for each fly
        G1.flyav(leg).data(:,ifly) = mean(G1.all(leg,ifly).data,2);
        G1.flyav(leg).control(ifly) = mean(G1.flyav(leg).data(control_ROI,ifly));
        G1.flyav(leg).stim(ifly) = mean(G1.flyav(leg).data(stim_ROI,ifly));
    end
    G1.avg(:,leg) = mean(G1.flyav(leg).data,2);
    G1.err(:,leg) = sem(G1.flyav(leg).data,2);
    % t-test for each leg:
    [~,G1.p(leg)] = ttest(G1.flyav(leg).control,G1.flyav(leg).stim);
end

% STATS: determine significance using a paired t-test w/mulitple corrections
p_val = G1.p';
fprintf('\n Significance vs multiple comparisons\n')
[p_err,p_loc] = sort(p_val);
for idx = 1:numel(p_val)
    q = (p_err((idx)) > 0.05/(length(p_err) +1 - idx));
    r = (p_err((idx)) > (idx/length(p_err))*.05);
    fprintf(['\n ' leg_labels{p_loc(idx)} ' significant? ' num2str(~q) ' vs ' num2str(~r) ])
    sig(p_loc(idx)) = ~q;
end
fprintf('\n -- Stats Done -- \n')
G1.sig = sig;

fig = getfig('Change in stride distance',1);
    hold on
    % Plotting:
    for leg = nLegs
        y1 = G1.flyav(leg).control;
        x1 = (leg*2-1)*ones(length(y1),1);
        scatter(x1, y1, SZ, Color(leg_colors{leg}), 'filled')
        y2 = G1.flyav(leg).stim;
        x2 = (leg*2)*ones(1,length(y2));
        scatter(x2, y2, SZ, Color(leg_colors{leg}), 'filled')
        % if paired data, plot connecting line by trial
        if paired==true
          for ii = 1:length(x1)
            plot([x1(ii),x2(ii)], [y1(ii), y2(ii)], ...
                'color', Color(leg_colors{leg}), 'linewidth', LW)
          end
        end
    end
    % labels and formatting:
    xlim([0,13])
    ylabel('Avg step distance (au)')
    xlabel('Legs')
    ax = gca;
    set(ax, 'XTick', ((nLegs.*2)-0.5))
    set(ax, 'XTickLabels', fly(1).leg_labels)    
    
    % stats 
    ax = gca;
    ub = ax.YTick;
    rng = range([ub(1),ub(end)])*0.05; %significance star offset
    pk = ub(end)+mean(diff(ub))+rng;
    for leg = nLegs
        if G1.sig(leg)==true
            x = [leg*2-1,leg*2];
            plot(x, [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
            s = scatter(leg*2-0.5,pk+rng, 50, '*');
            s.CData = Color(baseColor);
        end
    end
    
title({FilePath.structure_name; 'Avg stride distance';...
      [behaviorState ' flies | conds: ' num2str(G1.conds)]});

fig = formatFig(fig,color_opt); 
save_figure(fig, [fig_dir, '\' behaviorState ' change in step distance by leg-conds ' num2str(G1.conds)], fileType);

clearvars('-except',initial_vars{:})


%% Stride distance (combined across legs)
euclidian_dist = true;     % euclidian or total step distance
behaviorState = 'walking';  % behavior group to isolate e.g. 'walking', 'stationary','all'
shading = true;            % shading vs no shading
allTrials = false;           % plot lines for each trial, not just the avg
color_opt = true;          % false = white background, true = black background formatting
norm_axes = true;           % smae axes range for both time course figures
paired = true;              % line between control-stim points comparison
print_trials = true;        % print the ID of all trials being selected *good for errors
param = fly(1).param;       % assumes all flies have the same params (#conds,etc)
num = NUM(param);           % pull num reps/conds from a fly
fileType = '-png';          % figure save type (e.g. '-png' or '-pdf')
LP = 5;                     % leg position for movement trace (e.g. tarsus stride dist)

% define the two ROIs for comparison
G1.conds = [1,8];           % control 0ms groups
G2.conds = [7,14];          % laser 720 ms condition 
nLegs = 1:6;                % legs to use
G1.color = 'black';         % control color
G2.color = 'teal';          % 'stim' group color
if color_opt==true
    baseColor = 'white';
    l_color = 'black';
    if strcmpi(G1.color, 'black')
        G1.color = baseColor;
    end
    if strcmpi(G2.color, 'black')
        G2.color = baseColor;
    end
else 
    baseColor = 'black';
    l_color = 'white';
end

% soft parameters
light_on = param.basler_delay;
G1.laser = param.conds_matrix(G1.conds(1)).opto;
G2.laser = param.conds_matrix(G2.conds(1)).opto;
G1.light_off = param.conds_matrix(G1.conds(1)).opto;
G2.light_off = param.conds_matrix(G2.conds(1)).opto;
rec = param.basler_length;
time = linspace(-light_on,(rec-light_on),(fps*rec+1));
LW = 1;
SZ = 50;

% auto will be based on the first condition listed...
switch questdlg({'Auto fit G1 stim ROI laser length from G2?';...
                 '(If G1 isn''t the 0 ms condition, then no)'})
    case 'Yes' 
        control_ROI = 30:round(light_on*fps); %strt 30 to avoid edge effect
        stim_ROI = round(light_on*fps)+1:round(light_on*fps)+1+round(param.conds_matrix(G2.conds(1)).opto*fps);
    case 'No' %manually select prestim and stim regions
        control_ROI = 1:150;     % first 0.5 sec (control)
        stim_ROI = 151:366;   % 720 ms of stimulus
end

% Graph parameters
xlimits = [-0.4,1.6];
nrow = 2;
ncol = 3;
G1sp = [1,2];   %subplot idx for G1 timecourse
G2sp = [4,5];   %subplot idx for G2 timecourse
scsb = [3,6];   %subplot idx for change in freq scatter
plotIdx(1).idx = G1sp;
plotIdx(2).idx = G2sp;
plotIdx(3).idx = scsb;


% G1: Calculate step frequency
disp('G1 selected trials')
for ifly = 1:nflies
    % skip flies without sphere fits
    if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(G1.conds)
        cond = G1.conds(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
            
            % check for poor fit data:
            if isempty(step(ifly).ballFit(cond,rep).fit)...
               || step(ifly).ballFit(cond,rep).fit==false
                 continue;
            end
            
            % Print stats about each trial that is selected: (error catching)
            if print_trials==true
                state = fly(ifly).behavior.behavior{cond,rep};
                id = ['fly-' num2str(ifly) ', cond-' num2str(cond) ', rep-' num2str(rep) ' : '];
                disp([id state])
            end
            
            idx = idx+1;
            % pull the logical of swing|stance per frame
            raw = step(ifly).stance(cond,rep).loc'; %flip orientation    
            
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
                data = fly(ifly).pose_3d.Leg(LP,leg).data{cond,rep}.all;
                if euclidian_dist==true
                    % find the tarsus position for each stride:
                    stride(leg).start_pos = data(stride(leg).loc(:,1),:);
                    stride(leg).end_pos = data(stride(leg).loc(:,2),:);
                    % find the euclidian distance between the start and end of the stride
                    for ii = 1:nstrides
                        stride(leg).length(ii,1) = pdist2(stride(leg).start_pos(ii,:),stride(leg).end_pos(ii,:));
                    end
                else % total distance traveled:
                    stride(leg).start_pos = data(stride(leg).loc(:,1),:);
                    stride(leg).end_pos = data(stride(leg).loc(:,2),:);
                    for ii = 1:nstrides
                        frames = stride(leg).loc(ii,1):stride(leg).loc(ii,2);
                        step_pos = data(frames,:);
                        temp = pdist2(step_pos,step_pos);
                        for tt = 1:length(frames)-1
                            d(tt) = temp(tt,tt+1);
                        end
                        stride(leg).length(ii,1) = sum(d);
                    end
                end
                % extrapolate and fill the graps between stance (i.e. swing)
                [X,Y,y] = deal([]);
                Y = nan(param.basler_length*fps+1,1);
                X = [1; stride(leg).loc(:,1); param.basler_length*fps+1];
                for ii = 1:length(X)-1
                    f = stride(leg).length;
                    if ii == length(X)-1; val = f(ii-1); 
                    else val = f(ii); end
                    strt = X(ii);
                    fin = X(ii+1);
                    Y(strt:fin-1) = val*ones(fin-strt,1);
                end
                Y(end) = Y(end-1);  
                y = smooth(Y,fps/10,'lowess'); 
                stride(leg).filled = y;
            end
            % find the average stride leg
            % assign the frequency data to the Gait structure:
            G1.gait(ifly).Stride(idx).data = stride;
        end
    end
end


% G2: Calculate step frequency
disp('G1 selected trials')
for ifly = 1:nflies
    % skip flies without sphere fits
    if isempty(step(ifly).ball); continue; end
    idx = 0;
    for n = 1:length(G2.conds)
        cond = G2.conds(n);
        for rep = 1:num.reps
            % check for a behavior match:
            if strcmpi(behaviorState,'all')
            elseif ~strcmpi(fly(ifly).behavior.behavior(cond,rep),behaviorState)
              continue;
            end
            
            % check for poor fit data:
            if isempty(step(ifly).ballFit(cond,rep).fit)...
               || step(ifly).ballFit(cond,rep).fit==false
                 continue;
            end
            
            % Print stats about each trial that is selected: (error catching)
            if print_trials==true
                state = fly(ifly).behavior.behavior{cond,rep};
                id = ['fly-' num2str(ifly) ', cond-' num2str(cond) ', rep-' num2str(rep) ' : '];
                disp([id state])
            end
            
            idx = idx+1;
            % pull the logical of swing|stance per frame
            raw = step(ifly).stance(cond,rep).loc'; %flip orientation    
            
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
                data = fly(ifly).pose_3d.Leg(LP,leg).data{cond,rep}.all;
                if euclidian_dist==true
                    % find the tarsus position for each stride:
                    stride(leg).start_pos = data(stride(leg).loc(:,1),:);
                    stride(leg).end_pos = data(stride(leg).loc(:,2),:);
                    % find the euclidian distance between the start and end of the stride
                    for ii = 1:nstrides
                        stride(leg).length(ii,1) = pdist2(stride(leg).start_pos(ii,:),stride(leg).end_pos(ii,:));
                    end
                else % total distance traveled:
                    stride(leg).start_pos = data(stride(leg).loc(:,1),:);
                    stride(leg).end_pos = data(stride(leg).loc(:,2),:);
                    for ii = 1:nstrides
                        frames = stride(leg).loc(ii,1):stride(leg).loc(ii,2);
                        step_pos = data(frames,:);
                        temp = pdist2(step_pos,step_pos);
                        for tt = 1:length(frames)-1
                            d(tt) = temp(tt,tt+1);
                        end
                        stride(leg).length(ii,1) = sum(d);
                    end
                end
                % extrapolate and fill the graps between stance (i.e. swing)
                [X,Y,y] = deal([]);
                Y = nan(param.basler_length*fps+1,1);
                X = [1; stride(leg).loc(:,1); param.basler_length*fps+1];
                for ii = 1:length(X)-1
                    f = stride(leg).length;
                    if ii == length(X)-1; val = f(ii-1); 
                    else val = f(ii); end
                    strt = X(ii);
                    fin = X(ii+1);
                    Y(strt:fin-1) = val*ones(fin-strt,1);
                end
                Y(end) = Y(end-1);  
                y = smooth(Y,fps/10,'lowess'); 
                stride(leg).filled = y;
            end
            % find the average stride leg
            % assign the frequency data to the Gait structure:
            G2.gait(ifly).Stride(idx).data = stride;
        end
    end
end

% G1: Pool step frequency data across flies           
[flycount,idx] = deal(0); G1.dist = [];
for ifly = 1:nflies
    ntrials = length(G1.gait(ifly).Stride);
    if isempty(ntrials); continue; end
    flycount = flycount+1;
    for ii = 1:ntrials
      idx = idx+1;  
      for leg = nLegs
        f = G1.gait(ifly).Stride(ii).data(leg).filled;
        G1.dist(leg).data(:,idx) = f;
      end
    end
end
disp(['G1 N-flies : ' num2str(flycount)])
disp(['G1 N-trials : ' num2str(idx)])
G1.ntrials = idx;


% G2: Pool step frequency data across flies           
[flycount,idx] = deal(0); G2.dist = [];
for ifly = 1:nflies
    ntrials = length(G2.gait(ifly).Stride);
    if isempty(ntrials); continue; end
    flycount = flycount+1;
    for ii = 1:ntrials
      idx = idx+1;  
      for leg = nLegs
        f = G2.gait(ifly).Stride(ii).data(leg).filled;
        G2.dist(leg).data(:,idx) = f;
      end
    end
end
disp(['G2 N-flies : ' num2str(flycount)])
disp(['G2 N-trials : ' num2str(idx)])
G2.ntrials = idx;


% G1: Find the averages (avg legs, then across flies)
temp = [];
for n = 1:G1.ntrials
  for leg = nLegs
    temp(:,leg) = G1.dist(leg).data(:,n);
  end
  G1.trialmean(:,n) = mean(temp,2);
end
G1.avg = mean(G1.trialmean,2);
G1.err = std(G1.trialmean,0,2);
%pull control vs stim avg for scatter plot
G1.control = mean(G1.trialmean(control_ROI,:));
G1.stim = mean(G1.trialmean(stim_ROI,:));

% G2: Find the averages (avg legs, then across flies)
temp = [];
for n = 1:G2.ntrials
  for leg = nLegs
    temp(:,leg) = G2.dist(leg).data(:,n);
  end
  G2.trialmean(:,n) = mean(temp,2);
end
G2.avg = mean(G2.trialmean,2);
G2.err = std(G2.trialmean,0,2);
%pull control vs stim avg for scatter plot
G2.control = mean(G2.trialmean(control_ROI,:));
G2.stim = mean(G2.trialmean(stim_ROI,:));


% STATS: determine significance using a paired t-test w/mulitple corrections
[~,G1.p] = ttest(G1.control, G1.stim);
[~,G2.p] = ttest(G2.control, G2.stim);
p_val = [G1.p;G2.p];

% multiple comparisons test:
fprintf('\n Significance vs multiple comparisons\n')
[p_err,p_loc] = sort(p_val);
for idx = 1:numel(p_val)
    q = (p_err((idx)) > 0.05/(length(p_err) +1 - idx));
    r = (p_err((idx)) > (idx/length(p_err))*.05);
    fprintf(['\n Group ' num2str(p_loc(idx)) ' significant? ' num2str(~q) ' vs ' num2str(~r) ])
    sig(p_loc(idx)) = ~q;
    switch p_loc(idx)
        case 1
            G1.sig = ~q;
        case 2
            G2.sig = ~q;
    end
end
fprintf('\n -- Stats Done -- \n')

% PLOT HERE:
fig = getfig('Step distance two groups',1);
% G1 time course
subplot(nrow,ncol,G1sp)
    hold on
    y = G1.avg;
    err = G1.err;
    % plot individual traces
    if allTrials==true
       plot(time, G1.trialmean, 'color', Color('grey'), 'linewidth', 0.5)
    end

    % plot error region
    if shading==true
        fill_data = error_fill(time, y, err);
        h = fill(fill_data.X, fill_data.Y, Color(G1.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.4)
    end
    % plot error edges and avg line *not optional
    plot(time, y, 'color', Color(G1.color), 'linewidth', LW)
    plot(time, y+err, 'color', Color(G1.color), 'linewidth', 0.5)
    plot(time, y-err, 'color', Color(G1.color), 'linewidth', 0.5)

    % labels etc
    xlabel('time (s)')
    ylabel('Step distance (au)')
    title(['Group 1 : conds ' num2str(G1.conds)])
    % laser region lines:
    xlim(xlimits)
    if norm_axes==false
        y1 = rangeLine(fig);
        plot([0,G1.light_off], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
    end
    %pull ylimits to later normalize
    ylimits(1,:) = ylim;

% G2 time course
subplot(nrow,ncol,G2sp)
    hold on
    y = G2.avg;
    err = G2.err;
    % plot individual traces
    if allTrials==true
       plot(time, G2.trialmean, 'color', Color('grey'), 'linewidth', 0.5)
    end

    % plot error region
    if shading==true
        fill_data = error_fill(time, y, err);
        h = fill(fill_data.X, fill_data.Y, Color(G2.color), 'EdgeColor','none');
        set(h, 'facealpha', 0.4)
    end
    % plot error edges and avg line *not optional
    plot(time, y, 'color', Color(G2.color), 'linewidth', LW)
    plot(time, y+err, 'color', Color(G2.color), 'linewidth', 0.5)
    plot(time, y-err, 'color', Color(G2.color), 'linewidth', 0.5)

    % labels etc
    xlabel('time (s)')
    ylabel('Step distance (au)')
    title(['Group 2 : conds ' num2str(G2.conds)])
    % laser region lines:
    xlim(xlimits)
    if norm_axes==false
        y1 = rangeLine(fig);
        plot([0,G2.light_off], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
    end
    %pull ylimits to later normalize
    ylimits(2,:) = ylim;
    
%standardize axes for time course figures if opted
    if norm_axes==true
        ub = max(ylimits(:,2));
        lb = min(ylimits(:,1));
        subplot(nrow,ncol,[1,2])
            ylim([lb,ub])
            % plot stim region
            y1 = rangeLine(fig,2);
            plot([0,G1.light_off], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
            l = legend({'std', 'avg'}); set(l,'Color',Color(l_color));
            [l.EdgeColor,l.TextColor] = deal(Color(baseColor));
            
        subplot(nrow,ncol,[4,5])
            ylim([lb,ub])
            % plot stim region
            y1 = rangeLine(fig,2);
            plot([0,G2.light_off], [y1,y1],'color', Color(param.LED_light), 'linewidth', 5)
            l = legend({'std', 'avg'}); set(l,'Color',Color(l_color));
            [l.EdgeColor,l.TextColor] = deal(Color(baseColor));
    end

% control vs stim ROI both groups:
subplot(nrow,ncol,scsb); hold on
    % G1
    x1 = 1*ones(1,G1.ntrials);
    x2 = 2*ones(1,G1.ntrials);
    scatter(x1, G1.control, SZ, Color(G1.color), 'filled')
    scatter(x2, G1.stim, SZ, Color(G1.color), 'filled')
    if paired==true
      for n = 1:G1.ntrials
        plot([x1(n),x2(n)], [G1.control(n), G1.stim(n)], ...
            'color', Color(G1.color), 'linewidth', LW)
      end
    end
    % G2
    x1 = 3*ones(1,G2.ntrials);
    x2 = 4*ones(1,G2.ntrials);
    scatter(x1, G2.control, SZ, Color(G2.color), 'filled')
    scatter(x2, G2.stim, SZ, Color(G2.color), 'filled')
    if paired==true
      for n = 1:G2.ntrials
        plot([x1(n),x2(n)], [G2.control(n), G2.stim(n)], ...
            'color', Color(G2.color), 'linewidth', LW)
      end
    end
    % labels etc
    xlim([0,5])
    title('Control vs stim')
    ylabel('Mean step distance (au)')
    xlabel('Groups')
    ax = gca;
    set(ax, 'XTick', [1.5,3.5])
    set(ax, 'XTickLabels', {'G1', 'G2'})
    title({'Step dist control vs. stim';...
          ['P=' num2str(G1.p) ' | P=' num2str(G2.p)]})
    % stats 
    y_lims = ylim; ax = gca;
    ub = ax.YTick;
    rng = range([ub(1),ub(end)])*0.05; %significance star offset
    pk = ub(end)+mean(diff(ub))+rng;
    if G1.sig==true
        plot([1,2], [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
        s = scatter(1.5,pk+rng, 50, '*');
        s.CData = Color(baseColor);
    end
    if G2.sig==true
        plot([3,4], [pk,pk], 'color', Color(baseColor), 'linewidth', LW+1)
        s = scatter(3.5,pk+rng, 50, '*');
        s.CData = Color(baseColor);
    end
   
fig = formatFig(fig,color_opt,[nrow,ncol],plotIdx); 

save_figure(fig, [fig_dir, '\' behaviorState ' step distance between groups'], fileType);

clearvars('-except',initial_vars{:})


%%  Extraneous



























