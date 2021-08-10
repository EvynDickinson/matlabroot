%% Load Fly Structure (already built with "MIDSTIM_Step_4_groupflies.m")
clc; close all; clear all

[filename, directory] = uigetfile('*.mat', 'Select your Fly Structure');
fullname = fullfile(directory, filename); 
load(fullname);
try fly = FLY; catch
end
if isfield(fly(1), 'parameters') 
    for kk = 1:length(fly)
        fly(kk).param = fly(kk).parameters;
    end
    fly = rmfield(fly,'parameters');
end
kk = 1;
parameters = fly(1).param;
fprintf(['\n Loaded ' filename '\n'])
setGlobalx(parameters)

% load labels & condition parameters
[num, Fictrac, Matlab, labels] = NUM(fly(1).param);  
num.fly = length(fly);
Type = {'Stim', 'Control'};

% Create a Figures Folder for this structure:
% figures_dir = [directory 'MN line figures\' filename(1:end-4) '\'];
figures_dir = [directory 'Interneuron Lines/' filename(1:end-4) '/'];
if ~isfolder(figures_dir)
    mkdir(figures_dir)
end 

cond_figures_dir = [figures_dir 'Conditions/'];
if ~isfolder(cond_figures_dir)
    mkdir(cond_figures_dir)
end 
clc
% ELIMINATE OUTLIERS
fly = eliminate_speed_outliers(fly, 5);

% LOAD BEHAVIOR CLASSIFICATION DATA:
beep
switch questdlg('Load behavior classification data?')
    case 'Yes'
        load([directory, '/behavior class/', structure_name,  ' behavior class'])
end

% TRAJECTORY PATHS AND HEATMAPS
ROI = 1:60;
beep
switch questdlg('Load former trajectory info?')
    case 'Yes' 
        load([figures_dir, 'Trajectory Data'])
    case {'No', 'Cancel'}
        return
end


% USE THE 'GROUP' STRUCTURE FROM 'Behavior_Categorization.m'

% Convert behvariors into numbers and then a filter:
for kk = 1:num.fly
    for cond = 1:num.conds
        for rep = 1:num.reps
            state = group(kk).behavior{cond, rep};
            phase = group(kk).phase{cond, rep};
            if ~ischar(state)
                state = cell2mat(state);
            end
            if ~ischar(phase)
                phase = cell2mat(phase);
            end
            switch state
                case {'stationary', 'Stationary'}
                    STATE = 1;
                case {'walking', 'Walking'}
                    STATE = 2;
                case {'grooming', 'Grooming'}
                    STATE = 3;
                case {'other', 'Other'}
                    STATE = 4;
            end
            switch phase
                case {'stance', 'Stance'}
                    PHASE = 1;
                case {'swing', 'Swing'}
                    PHASE = 2;
            end 
            group(kk).STATE(cond, rep) = STATE;
            group(kk).PHASE(cond, rep) = PHASE;
        end
    end
    group(kk).walking = (group(kk).STATE==2);
    group(kk).stationary = (group(kk).STATE==1);
    group(kk).grooming = (group(kk).STATE==3);
    group(kk).other = (group(kk).STATE==4);
    group(kk).stance = (group(kk).PHASE==1);
    group(kk).swing = (group(kk).PHASE==2);
end


switch questdlg('Load joint angle info?') % JA_Data
    case 'Yes' 
        load([figures_dir, structure_name, ' joint angle data'])
    case {'No', 'Cancel'}
        return
end


%% Average Trajectory for walkers and non walkers: (finished)
clear stim ctrl cw ccw
p.color_idx = Color('LimeGreen'); 
% p.color_idx = Color('Red'); 
BClass = {'Walking', 'Nonwalking'};
for classification = 0:1 %1 = selecting only for non walkers; 0=walkers
    p.wide = 3;
    p.thin = 1;
    p.dotwide = 40;
    p.dotthin = 20;
    p.cap = 0;
    x_lim = [-.6, 0.6]; 
    y_lim = [-.1, 2];

    fig = getfig('Destination Trials');
    % STIM LENGTH TRIAL:

    for cond = 2:7 %WALKING
        idx = cond-1;
        stim(cond).x.data = []; stim(cond).y.data = [];
        ctrl(cond).x.data = []; ctrl(cond).y.data = [];
        ROI = ceil(parameters.light_length(cond)*parameters.fps);
        for kk = 1:num.fly % concatenate all the data into one vector
            for rep = 1:num.reps
                % stim X
                    cw = traj(cond).x(kk, rep).data(ROI);
                    if group(kk).walking(cond,rep)==classification
    %                     fly(kk).Control.mean_speed(cond).data(rep) < min_speed
                        cw = NaN;
                    end
                    ccw = traj(cond+7).x(kk, rep).data(ROI)*-1;
                    if group(kk).walking(cond+7,rep)==classification
    %                     fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
                        ccw = NaN;
                    end
                stim(cond).x.data = [stim(cond).x.data, cw, ccw];
                % stim Y
                    cw = traj(cond).y(kk, rep).data(ROI);
                    if group(kk).walking(cond,rep)==classification
    %                     fly(kk).Control.mean_speed(cond).data(rep) < min_speed
                        cw = NaN;
                    end
                    ccw = traj(cond+7).y(kk, rep).data(ROI);
                    if group(kk).walking(cond+7,rep)==classification
    %                     fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
                        ccw = NaN;
                    end
                stim(cond).y.data = [stim(cond).y.data, cw, ccw];    
                % control X
                cntrlcond = 1;
                    cw = traj(cntrlcond).x(kk, rep).data(ROI);
                    if group(kk).walking(cntrlcond,rep)==classification
    %                     fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
                        cw = NaN;
                    end
                    ccw = traj(cntrlcond+7).x(kk, rep).data(ROI)*-1;
                    if group(kk).walking(cntrlcond+7,rep)==classification
    %                     fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
                        ccw = NaN;
                    end
                ctrl(cond).x.data = [ctrl(cond).x.data, cw, ccw];
                % control Y
                    cw = traj(cntrlcond).y(kk, rep).data(ROI);
                    if group(kk).walking(cntrlcond,rep)==classification
    %                     fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
                        cw = NaN;
                    end
                    ccw = traj(cntrlcond+7).y(kk, rep).data(ROI);
                    if group(kk).walking(cntrlcond+7,rep)==classification
    %                     fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
                        ccw = NaN;
                    end
                ctrl(cond).y.data = [ctrl(cond).y.data, cw, ccw]; 
            end    
        end    
        % Stats:
        for ax = ['x', 'y']
            stim(cond).(ax).avg = nanmean(stim(cond).(ax).data);
            ctrl(cond).(ax).avg = nanmean(ctrl(cond).(ax).data); 
            stim(cond).(ax).err = sem(stim(cond).(ax).data,2);
            ctrl(cond).(ax).err = sem(ctrl(cond).(ax).data,2);
        end       
        % PLOT THE DATA:     
        figure(fig); subplot(1,2,1); hold all
        %stim data -- wide bars and point in light dependent color
           errorbar(stim(cond).x.avg, stim(cond).y.avg,... %coordinates
                 stim(cond).y.err, stim(cond).y.err,... %vertical error bar
                 stim(cond).x.err, stim(cond).x.err,... %horizontal error bar
                 'LineStyle','none', 'color', p.color_idx,...
                 'LineWidth', p.wide, 'CapSize', p.cap);
           scatter(stim(cond).x.avg, stim(cond).y.avg, p.dotwide, p.color_idx, 'filled')
        %control data -- wide bars|dot first in black and then thin in light dependent color
                 % WIDE
           errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
                 ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
                 ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
                 'LineStyle','none', 'color', 'k',...
                 'LineWidth', p.wide, 'CapSize', p.cap);
           scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotwide, 'k', 'filled')
    %              % THIN
    %        errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
    %              ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
    %              ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
    %              'LineStyle','none', 'color', p.color_idx(cond-1,:),...
    %              'LineWidth', p.thin, 'CapSize', p.cap);
    %        scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotthin, p.color_idx(idx,:), 'filled')     

        xlim(x_lim); ylim(y_lim)
        vline(0, 'r-');hline(0, 'r-')
        getaxes(fig, 10);
    end

    for cond = 16:21 %TURNING
        idx = cond-15;
        stim(cond).x.data = []; stim(cond).y.data = [];
        ctrl(cond).x.data = []; ctrl(cond).y.data = [];
        ROI = ceil(parameters.light_length(idx+1)*parameters.fps);
        for kk = 1:num.fly % concatenate all the data into one vector
            for rep = 1:num.reps
                % stim X
                    cw = traj(cond).x(kk, rep).data(ROI);
                    if group(kk).walking(cond,rep)==classification
    %                     fly(kk).Control.mean_speed(cond).data(rep) < min_speed
                        cw = NaN;
                    end
                    ccw = traj(cond+7).x(kk, rep).data(ROI)*-1;
                    if group(kk).walking(cond+7,rep)==classification
    %                     fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
                        ccw = NaN;
                    end
                stim(cond).x.data = [stim(cond).x.data, cw, ccw];
                % stim Y
                    cw = traj(cond).y(kk, rep).data(ROI);
                    if group(kk).walking(cond,rep)==classification
    %                     fly(kk).Control.mean_speed(cond).data(rep) < min_speed
                        cw = NaN;
                    end
                    ccw = traj(cond+7).y(kk, rep).data(ROI);
                    if group(kk).walking(cond+7,rep)==classification
    %                     fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
                        ccw = NaN;
                    end
                stim(cond).y.data = [stim(cond).y.data, cw, ccw];    
                % control X
                    cntrlcond = 15;
                    cw = traj(cntrlcond).x(kk, rep).data(ROI);
                    if group(kk).walking(cntrlcond,rep)==classification
    %                     fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
                        cw = NaN;
                    end
                    ccw = traj(cntrlcond+7).x(kk, rep).data(ROI)*-1;
                    if group(kk).walking(cntrlcond+7,rep)==classification
    %                     fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
                        ccw = NaN;
                    end
                ctrl(cond).x.data = [ctrl(cond).x.data, cw, ccw];
                % control Y
                    cw = traj(cntrlcond).y(kk, rep).data(ROI);
                    if group(kk).walking(cntrlcond,rep)==classification
    %                     fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
                        cw = NaN;
                    end
                    ccw = traj(cntrlcond+7).y(kk, rep).data(ROI);
                    if group(kk).walking(cntrlcond+7,rep)==classification
    %                     fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
                        ccw = NaN;
                    end
                ctrl(cond).y.data = [ctrl(cond).y.data, cw, ccw]; 
            end    
        end    
        % Stats:
        for ax = ['x', 'y']
            stim(cond).(ax).avg = nanmean(stim(cond).(ax).data);
            ctrl(cond).(ax).avg = nanmean(ctrl(cond).(ax).data); 
            stim(cond).(ax).err = sem(stim(cond).(ax).data,2);
            ctrl(cond).(ax).err = sem(ctrl(cond).(ax).data,2);
        end       
        % PLOT THE DATA:     
        figure(fig); subplot(1,2,2); hold all
        %stim data -- wide bars and point in light dependent color
           errorbar(stim(cond).x.avg, stim(cond).y.avg,... %coordinates
                 stim(cond).y.err, stim(cond).y.err,... %vertical error bar
                 stim(cond).x.err, stim(cond).x.err,... %horizontal error bar
                 'LineStyle','none', 'color', p.color_idx,...
                 'LineWidth', p.wide, 'CapSize', p.cap);
           scatter(stim(cond).x.avg, stim(cond).y.avg, p.dotwide, p.color_idx, 'filled')
        %control data -- wide bars|dot first in black and then thin in light dependent color
                 % WIDE
           errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
                 ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
                 ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
                 'LineStyle','none', 'color', 'k',...
                 'LineWidth', p.wide, 'CapSize', p.cap);
           scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotwide, 'k', 'filled')

        xlim(x_lim); ylim(y_lim)
        vline(0, 'r-');hline(0, 'r-')
        getaxes(fig, 10);
    end

    figure(fig)
    subplot(1,2,1)
        title({parameters.cross; 'Walking Cue'})
        xlabel('position (cm)'); ylabel('position (cm)')
    subplot(1,2,2)
        title({['Group: ' BClass{classification+1}]; 'Turning Cue'})          

    switch classification
        case 0
        save_figure(fig, [figures_dir, structure_name ' Average Trajectory walking selective']);
        case 1
        save_figure(fig, [figures_dir, structure_name ' Average Trajectory nonwalking selective']);
    end
end


%% Time-course velocity and rotvel figures comparing between types of stance (finished)
% for all light lengths together
% min_speed = 0.3;
x_data = -2:1/parameters.fps:2;
p.color_idx = [Color('black'); Color('DarkGreen', 'LimeGreen', 6)];
p.color_idx = [p.color_idx; p.color_idx];
p.lstyle = '-';
p.thin = 1;
% istate = {'stationary', 'walking', 'grooming', 'other'};
% iphase = {'stance', 'swing'};
istate = {'stationary', 'walking'};
iphase = {'stance', 'swing'};

for ii = 1:length(istate)
  for jj = 1:length(iphase)
    phase = iphase{jj};
    state = istate{ii};

    % Create data structure
    vari = {'speed', 'rotvelocity'};
    for pp = 1:2 % speed|rotvel
        switch pp
            case 1 % velocity
                adj = 1;
            case 2 % rotational velocity
                adj = -1;
        end
        param = vari{pp};
        for tt = 1:2 % Stim|Control
            idx = 0;
            for cond = [1:7, 15:21] % condition
                idx = idx+1;
                S(idx).(Type{tt}).(param).data = [];
                S(idx).(Type{tt}).(param).arena = [];
                for kk = 1:num.fly % fly number
                    cw = (fly(kk).(Type{tt}).(param)(cond).data);
    %                    filter =  fly(kk).Control.mean_speed(cond).data <= min_speed;                
                        filter = group(kk).(state)(cond,:);
                        cw(:,~filter) = NaN;
                        filter = group(kk).(phase)(cond,:);
                        cw(:,~filter) = NaN;
                    ccw = (fly(kk).(Type{tt}).(param)(cond+7).data);
    %                     filter = fly(kk).Control.mean_speed(cond+7).data <= min_speed;
                        filter = group(kk).(state)(cond+7,:);
                        ccw(:,~filter) = NaN;
                        filter = group(kk).(phase)(cond+7,:);
                        ccw(:,~filter) = NaN;
                    S(idx).(Type{tt}).(param).data = ...
                        [S(idx).(Type{tt}).(param).data, nanmean(cw, 2).*adj, nanmean(ccw,2)];
                end
                % Stats:
                S(idx).(Type{tt}).(param).avg = nanmean(S(idx).(Type{tt}).(param).data,2);
                S(idx).(Type{tt}).(param).err = sem(S(idx).(Type{tt}).(param).data,2);    
                S(idx).(Type{tt}).(param).arena = fly(kk).(Type{tt}).data(cond,1).raw(:,labels.arena_X);

            end
        end
        % Concatenate the Stim and Control data
        for idx = 1:14
            S(idx).all_raw.(param).avg = [S(idx).Control.(param).avg(1:60); S(idx).Stim.(param).avg]; 
                %adjust the values to be change in 'param' so they align to zero
                S(idx).all_raw.(param).offset = S(idx).Stim.(param).avg(1);
                S(idx).all.(param).avg = S(idx).all_raw.(param).avg-S(idx).all_raw.(param).offset;
            S(idx).all.(param).err = [S(idx).Control.(param).err(1:60); S(idx).Stim.(param).err];

            S(idx).all.(param).arena = [S(idx).Control.(param).arena(1:60); S(idx).Stim.(param).arena];
        end
    end

    % Plot the data: 
    fig = getfig('Time Course', 1);
    for pp = 1:2
        switch pp % Speed or Rotational Velocity
            case 1 
                param = 'speed'; adj2=1;
                idx = 0;
            case 2
                param = 'rotvelocity'; adj2=1;
                idx = 2;
        end

    % WALKING:
        subplot(2,2,idx+1); 
        hold all
        for cond = 1:7
        y_data = S(cond).all_raw.(param).avg;
        y_err = S(cond).all.(param).err;

        fill_data = error_fill(x_data, y_data.*adj2, y_err.*adj2);
            h = fill(fill_data.X, fill_data.Y, p.color_idx(cond,:), 'EdgeColor','none');
            set(h, 'facealpha', 0.2)
        plot(x_data, y_data.*adj2, 'LineStyle',p.lstyle, 'color', p.color_idx(cond,:),...
             'LineWidth', p.thin);
        end
        hold all
    % hline(0, 'k')
        vline(0, 'k:')
    %     plot(x_data, S(cond).all.(param).arena.*6, 'LineStyle',p.lstyle, 'color', 'r',...
    %          'LineWidth', p.thin+3)
         ylabel(param)
         title({structure_name; 'WALKING'; ['state: ' state ' LT1: ' phase]});%
         fig = getaxes(fig, 10);

    % TURNING:
        subplot(2,2,idx+2);   
        hold all
        for cond = 8:14
        y_data = S(cond).all_raw.(param).avg;
        y_err = S(cond).all.(param).err;

        fill_data = error_fill(x_data, y_data.*adj2, y_err.*adj2);
            h = fill(fill_data.X, fill_data.Y, p.color_idx(cond,:), 'EdgeColor','none');
            set(h, 'facealpha', 0.2)
        plot(x_data, y_data.*adj2, 'LineStyle',p.lstyle, 'color', p.color_idx(cond,:),...
             'LineWidth', p.thin);
        end
        hold all
    % hline(0, 'k')
        vline(0, 'k:')
    %     plot(x_data, S(cond).all.(param).arena.*1, 'LineStyle',p.lstyle, 'color', 'r',...
    %          'LineWidth', p.thin+3)  
         ylabel(param)
         title('TURNING')
        fig = getaxes(fig, 10);
    end


    % Find the min and max to set as axes scales:
    for idx = [1,3]
        yl = [];
        subplot(2,2,idx)
            yl = [yl, ylim];

        subplot(2,2,idx+1)
            yl = [yl, ylim];

        subplot(2,2,idx)
            ylim([min(yl), max(yl)]);
            vline(0,'k')
        subplot(2,2,idx+1)
            ylim([min(yl), max(yl)]);
            vline(0,'k')
    end
    for idx = 3:4
    subplot(2,2,idx)
    hline(0, 'k:')
    end

    save_figure(fig, [figures_dir, structure_name, ' speed and rotvel condition overlays ' state ' ' phase]);

  end
end


%% Percent time in behaviors across the three reps & within each fly in the structure (finished)

behavecode = {'stationary, walking, grooming, other'};
coloropts = {'gray', 'mediumblue', 'darkcyan', 'coral'};

BehaveDist.byfly = [];
BehaveDist.data = [];
for kk = 1:length(group)
    temp = [group(kk).STATE(:,1); group(kk).STATE(:,2); group(kk).STATE(:,3)];
    BehaveDist.byfly = [BehaveDist.byfly,temp];
    BehaveDist.data = [BehaveDist.data; group(kk).STATE];
end

% Percent time in state by fly:
for kk = 1:num.fly
    for ii = 1:4 %behavior states
        BehaveDist.sum.fly(ii,kk) = sum(BehaveDist.byfly(:,kk)==ii);
        BehaveDist.percent.fly(ii,kk) = ...
            (BehaveDist.sum.fly(ii,kk)/(num.reps*num.conds))*100;
    end
    labeldist{kk} = num2str(round(fly(kk).param.total_distance_traveled));
end
fig = getfig('Name', 1);

b = bar(BehaveDist.percent.fly', 'stacked');
for ii = 1:4
    b(ii).FaceColor = 'flat';
    b(ii).CData = Color(coloropts{ii});
end
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', labeldist)
xlabel({'Total distance traveled by fly (cm)'; 'Fly Number'})
ylabel({'% Time in activity'; '(Grey-stationary, blue-walking, green-grooming, orange-other)'})
title({structure_name; ' Distribution of activities by fly'})
save_figure(fig, [figures_dir, structure_name, ' fly behavior dist']);

% Percent time in state by rep all flies:
for rep = 1:num.reps
    for ii = 1:4
        BehaveDist.sum.rep(ii,rep) = sum(BehaveDist.data(:, rep)==ii);
        BehaveDist.percent.rep(ii,rep) = ...
                    BehaveDist.sum.rep(ii,rep)/(num.conds*num.fly)*100;
    end
end
fig = getfig;
b = bar(BehaveDist.percent.rep', 'stacked');
for ii = 1:3
    b(ii).FaceColor = 'flat';
    b(ii).CData = Color(coloropts{ii});
end
xlabel('Trial|Rep')
ylabel({'% Time in activity'; '(Grey-stationary, blue-walking, green-grooming, orange-other)'})
title({structure_name; ' Distribution of activities across time'})
save_figure(fig, [figures_dir, structure_name, ' rep behavior dist']);
  

%% Time-Course SPEED|ROT VEL Analysis grouped by initial state:(finished)
param = {'speed', 'rotvelocity'};

for condition = 1:2 %speed|rotvelocity
   
    fig = getfig('',1);
    coloropts = {'black', 'mediumblue', 'darkcyan', 'coral'}; %stationary, walking, grooming, other
    idx = 0;
    for icond = 1:num.conds %control for now
        idx = idx+1;
        for ii = 1:4
            InputData(ii).x = [];
            InputData(ii).y = [];
            InputData(ii).Ind = 1;
            InputData(ii).Color = Color(coloropts{ii});
        end

        for kk = 1:num.fly
          for rep = 1:num.reps
            % filter:
            filter = group(kk).STATE(icond,rep);
            % data:
            x = (-2:1/30:2)';
            y = [fly(kk).Control.(param{condition})(icond).data(2:end, rep);...
                 fly(kk).Stim.(param{condition})(icond).data(:, rep)];
            InputData(filter).x(:,InputData(filter).Ind) = x;
            InputData(filter).y(:,InputData(filter).Ind) = y;
            % set the next position
            InputData(filter).Ind = InputData(filter).Ind+1; 
          end
        end

        subplot(4,7,idx)
        fig = TimeCourseFigure(InputData, fig);
        getaxes(fig, 10);
        vline(0, 'k')
        title(getcond(icond))
    end
    subplot(4,7,1)
    xlabel('Time (s)')
    switch condition
        case 1 
            ylabel('cm/s')
        case 2
            ylabel('deg/s')
    end
    title({structure_name; getcond(1)})

    switch condition
        case 1
            save_figure(fig, [figures_dir, structure_name,...
                ' Behavior sorted time course Speed']);
        case 2
            save_figure(fig, [figures_dir, structure_name,...
                ' Behavior sorted time course Rot Velocity']);
    end

end


%% Combined CW and CCW state-dependent time course graphs:(finished)

param = {'speed', 'rotvelocity'};
for condition = 1:2 %speed|rotvelocity

    fig = getfig('',1);
    coloropts = {'black', 'mediumblue', 'darkcyan', 'coral'}; %stationary, walking, grooming, other
    idx = 0;
    for icond = [1:7, 15:21] %control for now
        idx = idx+1;
        for ii = 1:4
            InputData(ii).x = [];
            InputData(ii).y = [];
            InputData(ii).Ind = 1;
            InputData(ii).Color = Color(coloropts{ii});
        end

        for tt = 1:2
            adj = 1;
            switch tt % CW and CCW
                case 1
                    cond = icond;
                    if condition == 2
                        adj = -1;
                    end
                case 2
                    cond = icond+7;
            end
            for kk = 1:num.fly
              for rep = 1:num.reps
                % filter:
                filter = group(kk).STATE(cond,rep);
                % data:
                x = (-2:1/30:2)';
                y = [fly(kk).Control.(param{condition})(cond).data(2:end, rep);...
                     fly(kk).Stim.(param{condition})(cond).data(:, rep)];
                InputData(filter).x(:,InputData(filter).Ind) = x;
                InputData(filter).y(:,InputData(filter).Ind) = y.*adj;
                % set the next position
                InputData(filter).Ind = InputData(filter).Ind+1; 

              end
            end
        end

        subplot(2,7,idx)
        fig = TimeCourseFigure(InputData, fig);
        getaxes(fig, 10);
        vline(0, 'k')
        title(getcond(icond))
    end
    subplot(2,7,1)
    xlabel('Time (s)')
    switch condition
        case 1 
            ylabel('cm/s')
        case 2
            ylabel('deg/s')
    end
    title({structure_name; getcond(1)})

    switch condition
        case 1
            save_figure(fig, [figures_dir, structure_name,...
                ' Combined behavior sorted time course Speed']);
        case 2
            save_figure(fig, [figures_dir, structure_name,...
                ' Combined behavior sorted time course Rot Velocity']);
    end


end

%% Speed-dependent time-course of running flies (finished)

clear InputData
maxframes = num.fly*num.reps;
nframes = 0.2*num.fps;
BehaveDist.AllSpeeds = [];
% get speeds during the 200ms pre-stim:
for icond = 1:num.conds
    BehaveDist.Speeds(icond,1:maxframes) = NaN;
    temp = [];
    for kk = 1:num.fly
        filter = group(kk).walking(icond,:);
        y = mean(fly(kk).Control.speed(icond).data(end-nframes:end,:));
        fly(kk).behavior.screenspeed(icond, :) = y;
        fly(kk).behavior.group = group(kk);
        temp = [temp, y(filter)];
    end
    BehaveDist.Speeds(icond,1:length(temp)) = temp;
    BehaveDist.AllSpeeds = [BehaveDist.AllSpeeds, temp];
end

% histogram to divide by initial speeds:
nbins = 4;
% figure;
% histogram(BehaveDist.AllSpeeds,nbins); title(['Speed Distribution within ' structure_name])
% ylabel('Count'); xlabel('Mean speed during 200ms pre-stim (cm/s)')
% get the bin edges:
[~,E] = discretize(BehaveDist.AllSpeeds,nbins);
maxspeed = max(BehaveDist.AllSpeeds);

% Make the figures for both flies:

param = {'speed', 'rotvelocity'};
for condition = 1:2
    % Throw data into the order needed for generating a figure
    fig = getfig('Speed Based Analysis', 1);
    idx = 0;
    for icond = [1:7, 15:21]
        idx = idx+1;
        for tt = 1:2 % CW|CCW
            adj = 1;
            switch tt % CW and CCW
                case 1
                    cond = icond;
                    if condition == 2
                        adj = -1;
                    end
                case 2
                    cond = icond+7;
            end
            coloropts = {'Thistle', 'Orchid', 'BlueViolet', 'Indigo'};
%             coloropts = Color('Plum','Indigo', nbins);
            for ii = 1:nbins
                InputData(ii).Ind = 1;
                InputData(ii).x = (-2:1/num.fps:2)';
                InputData(ii).y = [];
                InputData(ii).Color = Color(coloropts{ii});
%                 InputData(ii).Color = coloropts(ii,:);
            end
            for kk = 1:num.fly
                for irep = 1:num.reps
                   meanspeed = fly(kk).behavior.screenspeed(cond,irep);
                    % create a sort filter for the speeds:
                    ii = discretize(meanspeed,E);
                    %filter data for walking only:
                    if group(kk).STATE(cond,irep) == 2
                        InputData(ii).y(:,InputData(ii).Ind) = ...
                            ([fly(kk).Control.(param{condition})(cond).data(2:61,irep);...
                             fly(kk).Stim.(param{condition})(cond).data(:,irep)]).*adj;
                        InputData(ii).Ind = InputData(ii).Ind+1;
                    end
                end
            end
        end
        % Generate a figure
        subplot(2,7,idx)
        fig = TimeCourseFigure(InputData, fig,1);
        getaxes(fig, 10);
        switch condition
            case 1 %speed
                ylim([0, maxspeed+0.25])
            case 2 %rotvel
        end
        vline(0, 'g')
        title(getcond(icond))
        vline((fly(1).param.conds_matrix(icond).opto), 'g')
    end
    subplot(2,7,1)
    xlabel('Time (s)')
    switch condition
        case 1 
            ylabel('cm/s')
        case 2
            ylabel('deg/s')
    end
    title({structure_name; getcond(1)})
    switch condition
        case 1
            save_figure(fig, [figures_dir, structure_name,...
                ' Speed dependent change in walking Speed']);
        case 2
            save_figure(fig, [figures_dir, structure_name,...
                ' Speed dependent change in walking Rot Velocity']);
    end
end
    


%% Loaded Vs Unloaded Time Course -- speed and rotational velocity (finished)
% Variables: 
clear InputData
param = {'speed', 'rotvelocity'};
coloropts = [Color('black'); Color('grey'); Color('saddlebrown'); Color('goldenrod')];

for condition = 1:2 % speed|rotvelocity
    % Throw data into the order needed for generating a figure
    fig = getfig('', 1);
    
    idx = 0;
    for icond = [1:7, 15:21] %combine CW&CCW
        idx = idx+1; % condition count (i.e. 1-7)      
        % Set the Input Data Structure:
        for ii = 1:4 % all types of loading
            InputData(ii).Ind = 1;
            InputData(ii).x = (-2:1/num.fps:2)';
            InputData(ii).y = [];
            InputData(ii).Color = coloropts(ii,:);
        end
        for tt = 1:2 % CW|CCW
            adj = 1;
            switch tt % CW and CCW
                case 1
                    cond = icond;
                    if condition == 2
                        adj = -1;
                    end
                case 2
                    cond = icond+7;
            end
            for kk = 1:num.fly %for all flies within the condition:
                for irep = 1:num.reps
                   istate = group(kk).STATE(icond,irep);
                   iphase = group(kk).PHASE(icond,irep);
                    % send data to appropriate grouping based on state and phase
                   switch iphase %loaded or unloaded
                       case 1 %loaded
                           if istate == 1 || istate == 3 %stationary or grooming
                               ii = 1;
                           elseif istate == 2 %walking
                               ii = 3;
                           end
                       case 2 %unloaded
                           if istate == 1 || istate == 3 %stationary or grooming
                               ii = 2;
                           elseif istate == 2 %walking
                               ii = 4;
                           end
                   end
                    % load the data into the InputData structure's
                    % appropriate field:
                    if group(kk).STATE(cond,irep) < 4
                        InputData(ii).y(:,InputData(ii).Ind) = ...
                            ([fly(kk).Control.(param{condition})(cond).data(2:61,irep);...
                             fly(kk).Stim.(param{condition})(cond).data(:,irep)]).*adj;
                        InputData(ii).Ind = InputData(ii).Ind+1;
                    end
                end
            end
        end
        % Generate a figure
        subplot(2,7,idx)
        fig = TimeCourseFigure(InputData, fig, 1);
        getaxes(fig, 10);
%         switch condition
%             case 1 %speed
%                 ylim([0, maxspeed+0.25])
%             case 2 %rotvel
%         end
        vline(0, 'g')
        title(getcond(icond))
        vline((fly(1).param.conds_matrix(icond).opto), 'g')
    end
    subplot(2,7,1)
    xlabel('Time (s)')
    switch condition
        case 1 
            ylabel('cm/s')
        case 2
            ylabel('deg/s')
    end
    title({structure_name; getcond(1)})
    switch condition
        case 1
            save_figure(fig, [figures_dir, structure_name,...
                ' Loaded vs Unloaded time course Speed']);
        case 2
            save_figure(fig, [figures_dir, structure_name,...
                ' Loaded vs Unloaded time course Rot Velocity']);
    end
end
    
beep

















