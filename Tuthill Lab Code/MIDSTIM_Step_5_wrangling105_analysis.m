 tic

%% Load Fly Structure (already built with "MIDSTIM_Step_4_groupflies.m")
clc; close all; clear all

min_speed = 0.3;


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
figures_dir = [directory 'Interneuron Lines\' filename(1:end-4) '\'];
if ~isfolder(figures_dir)
    mkdir(figures_dir)
end 

cond_figures_dir = [figures_dir 'Conditions\'];
if ~isfolder(cond_figures_dir)
    mkdir(cond_figures_dir)
end 
clc
% ELIMINATE OUTLIERS
fly = eliminate_speed_outliers(fly, 5);

% LOAD BEHAVIOR CLASSIFICATION DATA:
switch questdlg('Load behavior classification data?')
    case 'Yes'
        load([directory, '\behavior class\' structure_name, ' behavior class'])
end

%%
% TRAJECTORY PATHS AND HEATMAPS
ROI = 1:60;
switch questdlg('Load former trajectory info?')
    case 'Yes' 
        load([figures_dir, 'Trajectory Data'])
    case 'No'
    [trajectory_fig, heatmap_fig, traj] = MIDSTIM_gen_trajectories(fly, ROI);
    % Save the trajectory data! 
    save([figures_dir, 'Trajectory Data'], 'traj');
    save_figure(trajectory_fig, [figures_dir, structure_name, ' Flat Path Trajectories']);
    save_figure(heatmap_fig, [figures_dir, structure_name, ' Flat Path Heat Map']);
    close all
    case 'Cancel'
        return 
end

% figure(trajectory_fig) 
% for ii = 1:21
%     subplot(3,7,ii)
%     xlim([-4,4])
%     ylim([-2, 7])smooth
% end


%% Average Trajectory:
clear stim ctrl cw ccw
p.color_idx = Color('LimeGreen'); 
% p.color_idx = Color('Red'); 

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
                if fly(kk).Control.mean_speed(cond).data(rep) < min_speed
                    cw = NaN;
                end
                ccw = traj(cond+7).x(kk, rep).data(ROI)*-1;
                if fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
                    ccw = NaN;
                end
            stim(cond).x.data = [stim(cond).x.data, cw, ccw];
            % stim Y
                cw = traj(cond).y(kk, rep).data(ROI);
                if fly(kk).Control.mean_speed(cond).data(rep) < min_speed
                    cw = NaN;
                end
                ccw = traj(cond+7).y(kk, rep).data(ROI);
                if fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
                    ccw = NaN;
                end
            stim(cond).y.data = [stim(cond).y.data, cw, ccw];    
            % control X
            cntrlcond = 1;
                cw = traj(cntrlcond).x(kk, rep).data(ROI);
                if fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
                    cw = NaN;
                end
                ccw = traj(cntrlcond+7).x(kk, rep).data(ROI)*-1;
                if fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
                    ccw = NaN;
                end
            ctrl(cond).x.data = [ctrl(cond).x.data, cw, ccw];
            % control Y
                cw = traj(cntrlcond).y(kk, rep).data(ROI);
                if fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
                    cw = NaN;
                end
                ccw = traj(cntrlcond+7).y(kk, rep).data(ROI);
                if fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
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
    getaxes(fig, 10)
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
                if fly(kk).Control.mean_speed(cond).data(rep) < min_speed
                    cw = NaN;
                end
                ccw = traj(cond+7).x(kk, rep).data(ROI)*-1;
                if fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
                    ccw = NaN;
                end
            stim(cond).x.data = [stim(cond).x.data, cw, ccw];
            % stim Y
                cw = traj(cond).y(kk, rep).data(ROI);
                if fly(kk).Control.mean_speed(cond).data(rep) < min_speed
                    cw = NaN;
                end
                ccw = traj(cond+7).y(kk, rep).data(ROI);
                if fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
                    ccw = NaN;
                end
            stim(cond).y.data = [stim(cond).y.data, cw, ccw];    
            % control X
                cntrlcond = 15;
                cw = traj(cntrlcond).x(kk, rep).data(ROI);
                if fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
                    cw = NaN;
                end
                ccw = traj(cntrlcond+7).x(kk, rep).data(ROI)*-1;
                if fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
                    ccw = NaN;
                end
            ctrl(cond).x.data = [ctrl(cond).x.data, cw, ccw];
            % control Y
                cw = traj(cntrlcond).y(kk, rep).data(ROI);
                if fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
                    cw = NaN;
                end
                ccw = traj(cntrlcond+7).y(kk, rep).data(ROI);
                if fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
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
             % THIN
%        errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
%              ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
%              ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
%              'LineStyle','none', 'color', p.color_idx(idx,:),...
%              'LineWidth', p.thin, 'CapSize', p.cap);
%        scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotthin, p.color_idx(idx,:), 'filled')     

    xlim(x_lim); ylim(y_lim)
    vline(0, 'r-');hline(0, 'r-')
    getaxes(fig, 10)
end

figure(fig)
subplot(1,2,1)
    title({parameters.cross; 'Walking'})
    xlabel('cm'); ylabel('cm')
subplot(1,2,2)
    title('Turning')          

save_figure(fig, [figures_dir, structure_name ' Average Trajectory']);
beep  

% %% Average Trajectory:
% clear stim ctrl cw ccw
% p.color_idx = Color('DarkGreen', 'GreenYellow', 6);
% p.wide = 3;
% p.thin = 1;
% p.dotwide = 40;
% p.dotthin = 20;
% p.cap = 0;
% x_lim = [-.6, 0.6]; 
% y_lim = [-.1, 2];
% 
% fig = getfig('Destination Trials');
% % STIM LENGTH TRIAL:
% 
% for cond = 2:7 %WALKING
%     idx = cond-1;
%     stim(cond).x.data = []; stim(cond).y.data = [];
%     ctrl(cond).x.data = []; ctrl(cond).y.data = [];
%     ROI = ceil(parameters.light_length(cond)*parameters.fps);
%     for kk = 1:num.fly % concatenate all the data into one vector
%         for rep = 1:num.reps
%             % stim X
%                 cw = traj(cond).x(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cond).data(rep) < min_speed
%                     cw = NaN;
%                 end
%                 ccw = traj(cond+7).x(kk, rep).data(ROI)*-1;
%                 if fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
%                     ccw = NaN;
%                 end
%             stim(cond).x.data = [stim(cond).x.data, cw, ccw];
%             % stim Y
%                 cw = traj(cond).y(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cond).data(rep) < min_speed
%                     cw = NaN;
%                 end
%                 ccw = traj(cond+7).y(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
%                     ccw = NaN;
%                 end
%             stim(cond).y.data = [stim(cond).y.data, cw, ccw];    
%             % control X
%             cntrlcond = 1;
%                 cw = traj(cntrlcond).x(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
%                     cw = NaN;
%                 end
%                 ccw = traj(cntrlcond+7).x(kk, rep).data(ROI)*-1;
%                 if fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
%                     ccw = NaN;
%                 end
%             ctrl(cond).x.data = [ctrl(cond).x.data, cw, ccw];
%             % control Y
%                 cw = traj(cntrlcond).y(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
%                     cw = NaN;
%                 end
%                 ccw = traj(cntrlcond+7).y(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
%                     ccw = NaN;
%                 end
%             ctrl(cond).y.data = [ctrl(cond).y.data, cw, ccw]; 
%         end    
%     end    
%     % Stats:
%     for ax = ['x', 'y']
%         stim(cond).(ax).avg = nanmean(stim(cond).(ax).data);
%         ctrl(cond).(ax).avg = nanmean(ctrl(cond).(ax).data); 
%         stim(cond).(ax).err = sem(stim(cond).(ax).data,2);
%         ctrl(cond).(ax).err = sem(ctrl(cond).(ax).data,2);
%     end       
%     % PLOT THE DATA:     
%     figure(fig); subplot(1,2,1); hold all
%     %stim data -- wide bars and point in light dependent color
%        errorbar(stim(cond).x.avg, stim(cond).y.avg,... %coordinates
%              stim(cond).y.err, stim(cond).y.err,... %vertical error bar
%              stim(cond).x.err, stim(cond).x.err,... %horizontal error bar
%              'LineStyle','none', 'color', p.color_idx(cond-1,:),...
%              'LineWidth', p.wide, 'CapSize', p.cap);
%        scatter(stim(cond).x.avg, stim(cond).y.avg, p.dotwide, p.color_idx(idx,:), 'filled')
%     %control data -- wide bars|dot first in black and then thin in light dependent color
%              % WIDE
%        errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
%              ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
%              ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
%              'LineStyle','none', 'color', 'k',...
%              'LineWidth', p.wide, 'CapSize', p.cap);
%        scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotwide, 'k', 'filled')
%              % THIN
%        errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
%              ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
%              ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
%              'LineStyle','none', 'color', p.color_idx(cond-1,:),...
%              'LineWidth', p.thin, 'CapSize', p.cap);
%        scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotthin, p.color_idx(idx,:), 'filled')     
% 
%     xlim(x_lim); ylim(y_lim)
%     vline(0, 'r-');hline(0, 'r-')
%     getaxes(fig, 10)
% end
% 
% 
% for cond = 16:21 %TURNING
%     idx = cond-15;
%     stim(cond).x.data = []; stim(cond).y.data = [];
%     ctrl(cond).x.data = []; ctrl(cond).y.data = [];
%     ROI = ceil(parameters.light_length(idx+1)*parameters.fps);
%     for kk = 1:num.fly % concatenate all the data into one vector
%         for rep = 1:num.reps
%             % stim X
%                 cw = traj(cond).x(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cond).data(rep) < min_speed
%                     cw = NaN;
%                 end
%                 ccw = traj(cond+7).x(kk, rep).data(ROI)*-1;
%                 if fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
%                     ccw = NaN;
%                 end
%             stim(cond).x.data = [stim(cond).x.data, cw, ccw];
%             % stim Y
%                 cw = traj(cond).y(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cond).data(rep) < min_speed
%                     cw = NaN;
%                 end
%                 ccw = traj(cond+7).y(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cond+7).data(rep) < min_speed
%                     ccw = NaN;
%                 end
%             stim(cond).y.data = [stim(cond).y.data, cw, ccw];    
%             % control X
%                 cntrlcond = 15;
%                 cw = traj(cntrlcond).x(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
%                     cw = NaN;
%                 end
%                 ccw = traj(cntrlcond+7).x(kk, rep).data(ROI)*-1;
%                 if fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
%                     ccw = NaN;
%                 end
%             ctrl(cond).x.data = [ctrl(cond).x.data, cw, ccw];
%             % control Y
%                 cw = traj(cntrlcond).y(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cntrlcond).data(rep) < min_speed
%                     cw = NaN;
%                 end
%                 ccw = traj(cntrlcond+7).y(kk, rep).data(ROI);
%                 if fly(kk).Control.mean_speed(cntrlcond+7).data(rep) < min_speed
%                     ccw = NaN;
%                 end
%             ctrl(cond).y.data = [ctrl(cond).y.data, cw, ccw]; 
%         end    
%     end    
%     % Stats:
%     for ax = ['x', 'y']
%         stim(cond).(ax).avg = nanmean(stim(cond).(ax).data);
%         ctrl(cond).(ax).avg = nanmean(ctrl(cond).(ax).data); 
%         stim(cond).(ax).err = sem(stim(cond).(ax).data,2);
%         ctrl(cond).(ax).err = sem(ctrl(cond).(ax).data,2);
%     end       
%     % PLOT THE DATA:     
%     figure(fig); subplot(1,2,2); hold all
%     %stim data -- wide bars and point in light dependent color
%        errorbar(stim(cond).x.avg, stim(cond).y.avg,... %coordinates
%              stim(cond).y.err, stim(cond).y.err,... %vertical error bar
%              stim(cond).x.err, stim(cond).x.err,... %horizontal error bar
%              'LineStyle','none', 'color', p.color_idx(idx,:),...
%              'LineWidth', p.wide, 'CapSize', p.cap);
%        scatter(stim(cond).x.avg, stim(cond).y.avg, p.dotwide, p.color_idx(idx,:), 'filled')
%     %control data -- wide bars|dot first in black and then thin in light dependent color
%              % WIDE
%        errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
%              ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
%              ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
%              'LineStyle','none', 'color', 'k',...
%              'LineWidth', p.wide, 'CapSize', p.cap);
%        scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotwide, 'k', 'filled')
%              % THIN
%        errorbar(ctrl(cond).x.avg, ctrl(cond).y.avg,... %coordinates
%              ctrl(cond).y.err, ctrl(cond).y.err,... %vertical error bar
%              ctrl(cond).x.err, ctrl(cond).x.err,... %horizontal error bar
%              'LineStyle','none', 'color', p.color_idx(idx,:),...
%              'LineWidth', p.thin, 'CapSize', p.cap);
%        scatter(ctrl(cond).x.avg, ctrl(cond).y.avg, p.dotthin, p.color_idx(idx,:), 'filled')     
% 
%     xlim(x_lim); ylim(y_lim)
%     vline(0, 'r-');hline(0, 'r-')
%     getaxes(fig, 10)
% end
% 
% figure(fig)
% subplot(1,2,1)
%     title({parameters.cross; 'Walking'})
%     xlabel('cm'); ylabel('cm')
% subplot(1,2,2)
%     title('Turning')          
% 
% save_figure(fig, [figures_dir, structure_name ' Average Trajectory']);

%% Tuning Curves:
% min_speed = 0.5;
p.color_idx = Color('LimeGreen', 'DarkGreen',2);
p.wide = 3;
p.thin = 2;
p.dotwide = 50;
p.dotthin = 20;
p.cap = 0;
x_lim = [-.6, 0.6]; 
y_lim = [-.1, 2];
adj = 1; %for rotational velocity neg 1

x_data = parameters.light_length(2:7);

fig = getfig('Tuning Curves');
for pp = 1:2 
    switch pp
        case 1 % Speeeeeed
            param = 'speed';
            adj = 1; 
        case 2 % Rotational velocity
            adj = -1;
            param = 'rotvelocity';
    end

%***WALKING***
    for tt = 1:2
      for cond = 2:7 
        stim(cond).(param).data = [];
        ctrl(cond).(param).data = [];
        ROI = ceil(parameters.light_length(cond)*parameters.fps);
        switch tt
            case 1 %Stim data
                strt = 1;
                stp = ROI;
                p.fill = p.color_idx(2,:);
                p.fill2 = 'k';
                p.lstyle = '-';
            case 2 %Stim Recovery data
                strt = ROI;
                stp = ROI+parameters.fps;
                p.fill = 'w';
                p.fill2 = 'w';
                p.lstyle = ':';
        end
        for kk = 1:num.fly % concatenate all the data into one vector
            % Stim 
            cw = mean(fly(kk).Stim.(param)(cond).data(strt:stp,:));
                filter = fly(kk).Control.mean_speed(cond).data <= min_speed;
                cw(filter) = NaN;
            ccw = mean(fly(kk).Stim.(param)(cond+7).data(strt:stp,:));
                filter = fly(kk).Control.mean_speed(cond+7).data <= min_speed;
                ccw(filter) = NaN;
            stim(cond).(param).data = [stim(cond).(param).data, cw.*adj, ccw];

            % Control
            cntrlcond = 1;
            cw = mean(fly(kk).Stim.(param)(cntrlcond).data(strt:stp,:));
                filter = fly(kk).Control.mean_speed(cntrlcond).data <= min_speed;
                cw(filter) = NaN;
            ccw = mean(fly(kk).Stim.(param)(cntrlcond+7).data(strt:stp,:));
                filter = fly(kk).Control.mean_speed(cntrlcond+7).data <= min_speed;
                ccw(filter) = NaN;
            ctrl(cond).(param).data = [ctrl(cond).(param).data, cw.*adj, ccw];
        end  

        % Stats:
        stim(cond).(param).avg = nanmean(stim(cond).(param).data); 
        ctrl(cond).(param).avg = nanmean(ctrl(cond).(param).data); 
        stim(cond).(param).err = sem(stim(cond).(param).data,2);
        ctrl(cond).(param).err = sem(ctrl(cond).(param).data,2);

        % Restructure the data:
        s.stim.avg(cond-1) = stim(cond).(param).avg;
        s.ctrl.avg(cond-1) = ctrl(cond).(param).avg;
        s.stim.err(cond-1) = stim(cond).(param).err;
        s.ctrl.err(cond-1) = ctrl(cond).(param).err;
      end
      % PLOT THE DATA:     
        figure(fig); subplot(2,2,pp); hold all
        %stim data
        fill_data = error_fill(x_data, s.stim.avg, s.stim.err);
            h = fill(fill_data.X, fill_data.Y, p.color_idx(1,:), 'EdgeColor','none');
            set(h, 'facealpha', 0.2)
        plot(x_data, s.stim.avg, 'LineStyle',p.lstyle, 'color', p.color_idx(2,:),...
             'LineWidth', p.thin);
        scatter(x_data, s.stim.avg, p.dotwide, p.color_idx(2,:),...
             'MarkerFaceColor', p.fill);
        %control data
        fill_data = error_fill(x_data, s.ctrl.avg, s.ctrl.err);
            h = fill(fill_data.X, fill_data.Y, Color('grey'), 'EdgeColor','none');
            set(h, 'facealpha', 0.2)
        plot(x_data, s.ctrl.avg, 'LineStyle', p.lstyle, 'color', 'k',...
             'LineWidth', p.thin);
        scatter(x_data, s.ctrl.avg, p.dotwide, 'k',...
             'MarkerFaceColor', p.fill2);
        % xlim(x_lim); ylim(y_lim)     
        title({structure_name; 'WALKING'})
        ylabel(param)
        fig = getaxes(fig, 15);
    end
    
%***TURNING***
    for tt = 1:2
      for cond = 16:21 
        idx = cond-14;
        stim(cond).(param).data = [];
        ctrl(cond).(param).data = [];
        ROI = ceil(parameters.light_length(idx)*parameters.fps);
          switch tt
            case 1 %Stim data
                strt = 1;
                stp = ROI;
                p.fill = p.color_idx(2,:);
                p.fill2 = 'k';
                p.lstyle = '-';
            case 2 %Stim Recovery data
                strt = ROI;
                stp = ROI+parameters.fps;
                p.fill = 'w';
                p.fill2 = 'w';
                p.lstyle = ':';
        end
        for kk = 1:num.fly % concatenate all the data into one vector
            % Stim 
            cw = mean(fly(kk).Stim.(param)(cond).data(strt:stp,:));
                filter = fly(kk).Control.mean_speed(cond).data <= min_speed;
                cw(filter) = NaN;
            ccw = mean(fly(kk).Stim.(param)(cond+7).data(strt:stp,:));
                filter = fly(kk).Control.mean_speed(cond+7).data <= min_speed;
                ccw(filter) = NaN;
            stim(cond).(param).data = [stim(cond).(param).data, cw.*adj, ccw];

            % Control
            cntrlcond = 15;
            cw = mean(fly(kk).Stim.(param)(cntrlcond).data(strt:stp,:));
                filter = fly(kk).Control.mean_speed(cntrlcond).data <= min_speed;
                cw(filter) = NaN;
            ccw = mean(fly(kk).Stim.(param)(cntrlcond+7).data(strt:stp,:));
                filter = fly(kk).Control.mean_speed(cntrlcond+7).data <= min_speed;
                ccw(filter) = NaN;
            ctrl(cond).(param).data = [ctrl(cond).(param).data, cw.*adj, ccw];
        end  

        % Stats:
        stim(cond).(param).avg = nanmean(stim(cond).(param).data); 
        ctrl(cond).(param).avg = nanmean(ctrl(cond).(param).data); 
        stim(cond).(param).err = sem(stim(cond).(param).data,2);
        ctrl(cond).(param).err = sem(ctrl(cond).(param).data,2);

        % Restructure the data:
        s.stim.avg(idx-1) = stim(cond).(param).avg;
        s.ctrl.avg(idx-1) = ctrl(cond).(param).avg;
        s.stim.err(idx-1) = stim(cond).(param).err;
        s.ctrl.err(idx-1) = ctrl(cond).(param).err;
      end

      % PLOT THE DATA:
        figure(fig); subplot(2,2,pp+2); hold all
        %stim data
        fill_data = error_fill(x_data, s.stim.avg, s.stim.err);
            h = fill(fill_data.X, fill_data.Y, p.color_idx(1,:), 'EdgeColor','none');
            set(h, 'facealpha', 0.2)
        plot(x_data, s.stim.avg, 'LineStyle', p.lstyle, 'color', p.color_idx(2,:),...
             'LineWidth', p.thin);
        scatter(x_data, s.stim.avg, p.dotwide, p.color_idx(2,:),...
             'MarkerFaceColor', p.fill);

        %control data        
        fill_data = error_fill(x_data, s.ctrl.avg, s.ctrl.err);
            h = fill(fill_data.X, fill_data.Y, Color('Grey'), 'EdgeColor','none');
            set(h, 'facealpha', 0.2)
        plot(x_data, s.ctrl.avg, 'LineStyle', p.lstyle, 'color', 'k',...
             'LineWidth', p.thin);
        scatter(x_data, s.ctrl.avg, p.dotwide, 'k',...
             'MarkerFaceColor', p.fill2);
        title('TURNING')
        ylabel(param)
        fig = getaxes(fig, 15);
    end
end

save_figure(fig, [figures_dir, structure_name ' Tuning Curve light-length']); %

%% Time-course velocity and rotvel figures:
min_speed = 0.0;
x_data = -2:1/parameters.fps:2;
p.color_idx = [Color('black'); Color('DarkGreen', 'LimeGreen', 6)];
p.color_idx = [p.color_idx; p.color_idx];
p.lstyle = '-';
p.thin = 1;

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
                    filter = fly(kk).Control.mean_speed(cond).data <= min_speed;
                    cw(:,filter) = NaN;
                ccw = (fly(kk).(Type{tt}).(param)(cond+7).data);
                    filter = fly(kk).Control.mean_speed(cond+7).data <= min_speed;
                    ccw(:,filter) = NaN;
                S(idx).(Type{tt}).(param).data = [S(idx).(Type{tt}).(param).data, mean(cw, 2).*adj, mean(ccw,2)];
 
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

%
% Plot the data: 
fig = getfig('Time Course');
for pp = 1:2
    switch pp % Speed or Rotational Velocity
        case 1 
            param = 'speed'; adj2=10;
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
     title({structure_name; 'WALKING'})
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

save_figure(fig, [figures_dir, structure_name, ' speed and rotvel condition overlays']);

%% OLD FIGURES BELOW 

%% Figures of indivudual CW & CWW velocity and rotvelocity
%create label for combining CW and CCW data:
combined_label = create_combined_label(parameters);

% ROTATIONAL VELOCITY:
%cond placed as 2nd arg of MIDSTIM_(rot)velocity_fig returns cond & cond+7(CW&CCW)
%cond must be between 0-7;
[~, rotvelocity, rotvelocity_combined] = MIDSTIM_rotvelocity_fig(fly);

% SPEED | Velocity:
[~, velocity, velocity_combined] = MIDSTIM_velocity_fig(fly);

% % Save figures for all conditions with CW % CCW data visualized together
% fig_action = questdlg('Desired action for the rotational velocity figures? (''Esc'' for nothing)',...
%         'Rotation Velocity', 'View Only', 'Save Individual', 'Save and Append', 'View Only');
% if length(fig_action)>2 %anything but 'esc'
%     %generate the figures  
%     for cond = 1:7  
%         figure_name(cond) = {[cond_figures_dir filename(1:end-4) ' rotational velocity combined cond ' num2str(cond) '.pdf']};
%         fig(cond) = MIDSTIM_rotvel_and_speed_fig(fly, velocity, rotvelocity, cond);
%     end
%     %determine saving procedure
%   switch fig_action
%     case 'View Only'  
%        uiwait(fig(1)) %wait to proceed until the last figure is closed
%     case 'Save Individual'
%       for cond = 1:7
%         export_fig(fig(cond), figure_name{cond},'-pdf','-nocrop', '-r300' , '-painters', '-rgb');
%         close all
%       end
%     case 'Save and Append'
%       for cond = 7:-1:1
%         export_fig(fig(cond), figure_name{cond},'-pdf','-nocrop', '-r300' , '-painters', '-rgb');
%         close(fig(cond))
%       end
%       append_pdfs([figures_dir filename(1:end-4) ' rotational velocity combined.pdf'], figure_name{:})
%       close all
%       fprintf('\n Rotational velocity and speed figures appended!\n')
%   end
% end; clear figure_name fig_action

% Visualize speed|rotational velocity averages without slow flies
% Remove Slow Flies & show filtered vs unfiltered
% Fig   
[activity_level, fig] = Activity_Level_Index(fly, 0.3, 0.5); %Min activity percent, min speed
save_figure(fig, [figures_dir filename(1:end-4) ' filtered vs unfiltered  speed']);

% Fig: Speed of all conditions combined
fig = MIDSTIM_filtered_fig(fly, activity_level, velocity_combined);
save_figure(fig, [figures_dir filename(1:end-4) ' filtered speed all traces']);
 
% Filtered by Speed Struct with avg and err:
velocity_filtered = MIDSTIM_filter_combined_data(velocity_combined, activity_level);
rotvelocity_filtered = MIDSTIM_filter_combined_data(rotvelocity_combined, activity_level);

% %% Fourier Transform of Conditions Roll data: [mostly ignore this]
% % Fig
% cond = 5;
% fig = MIDSTIM_make_fourier_transform(fly, cond);
% save_figure(fig, [figures_dir filename(1:end-4) ' Roll Data Fourier Transform Cond ' num2str(cond)])
% 
% %% Analyze slope of effect and recovery after stimulus 
num.peak_select_buffer = 9;
% % ***unfiltered data***
% % rotvelocity -- all
% [rotvelocity_analysis_unfiltered, fig] = MIDSTIM_analysis_effect_n_recovery(rotvelocity_combined, fly, num.peak_select_buffer);
% save_figure(fig, [figures_dir filename(1:end-4) ' unfiltered rot velocity peak effect selection'])
% % speed -- all
% [velocity_analysis_unfiltered, fig] = MIDSTIM_analysis_effect_n_recovery(velocity_combined, fly, num.peak_select_buffer);
% save_figure(fig, [figures_dir filename(1:end-4) ' unfiltered velocity peak effect selection'])
% 
% % ***filtered data***
% %rotvelocity -- filtered
% [rotvelocity_analysis, fig] = MIDSTIM_analysis_effect_n_recovery(rotvelocity_filtered, fly, num.peak_select_buffer);
% save_figure(fig, [figures_dir filename(1:end-4) ' rot velocity peak effect selection'])
% %speed -- filtered
% [velocity_analysis, fig] = MIDSTIM_analysis_effect_n_recovery(velocity_filtered, fly, num.peak_select_buffer);
% save_figure(fig, [figures_dir filename(1:end-4) ' velocity peak effect selection'])
% 
% % Basic Trial Tuning Curves
% fig = MIDSTIM_basic_tuning_curves(rotvelocity_analysis, velocity_analysis, parameters);
% save_figure(fig, [figures_dir filename(1:end-4) ' basic tuning curves'])

% Single Trial Analysis - velocity
% Generate individual response (speed|velocity) analyses for each trial within each condition
for cond = 1:num.conds/2
   AD_velocity_analysis(cond) = MIDSTIM_trialsbytrial_analysis_effect_n_recovery...
                                (velocity_filtered, parameters, num.peak_select_buffer, cond);
end

% % Group the data into N bins based on their control steady state average
% num.bins = 6;
% [AD_velocity_single_trial_bins] = MIDSTIM_parameter_filter(AD_velocity_analysis, parameters, num.bins);
% 
% % Plot the single trial analysis tuning curves based on steady state
% % condition (3 figs: effect change, rate of response, rate of recovery)
% % replace the error bars with shaded regions
% fig = MIDSTIM_speed_based_tuning_curves_FIG(AD_velocity_analysis, AD_velocity_single_trial_bins, parameters);
% 
% Look at the area under the curves as a response analysis -- chunk out for
% various periods of time

% Overlay the average fly response to the various light lengths -- is the
% start the same? Do the reponses have similar trends? What about for speed
% based groups?

  
% Single Trial Analysis - rotational velocity
% Generate individual response (speed|velocity) analyses for each trial within each condition
for cond = 1:num.conds/2
   AD_rotvel_analysis(cond) = MIDSTIM_trialsbytrial_analysis_effect_n_recovery...
                             (rotvelocity_filtered, parameters, num.peak_select_buffer, cond);
end

% % Group the data into N bins based on their control steady state average
% num.bins = 6;
% [AD_rotvel_single_trial_bins] = MIDSTIM_parameter_filter(AD_rotvel_analysis, parameters, num.bins);
% 
% % Plot the single trial analysis tuning curves based on steady state
% % condition (3 figs: effect change, rate of response, rate of recovery)
% % replace the error bars with shaded regions
% 
% %use speed filter rather than rotvelcity for splitting up the data
% fig = MIDSTIM_speed_based_tuning_curves_FIG(AD_rotvel_analysis, AD_velocity_single_trial_bins, parameters);

% Overlay the average fly response to the various light lengths -- is the
% start the same? Do the reponses have similar trends? What about for speed
% based groups? 

% Rotational velocity aligned
fig = MIDSTIM_convergent_parameters_graphs(AD_rotvel_analysis, parameters);
save_figure(fig(1), [figures_dir filename(1:end-4) ' ' fig(1).Name ' light activation aligned']);
save_figure(fig(2), [figures_dir filename(1:end-4) ' ' fig(2).Name ' light activation aligned']);

% Speed aligned
fig = MIDSTIM_convergent_parameters_graphs(AD_velocity_analysis, parameters);
save_figure(fig(1), [figures_dir filename(1:end-4) ' ' fig(1).Name ' light activation aligned']);
save_figure(fig(2), [figures_dir filename(1:end-4) ' ' fig(2).Name ' light activation aligned']);

%% Manual Selection of Effects:
% ** NEED TO COMBINE THE TWO DIRECTIONS HERE ***
% Find the fly average and then the grouped average
for tt = 1:2
  for cond = 1:num.conds
    % Add average of each fly into the 'Grouped Data' structures
    for kk = 1:num.fly
       temp.avg_speed = mean(fly(kk).(Type{tt}).speed(cond).data,2);
       GD_velocity(cond).(Type{tt})(:,kk) = temp.avg_speed; 
       temp.avg_rot = mean(fly(kk).(Type{tt}).rotvelocity(cond).data,2);
       GD_rotvelocity(cond).(Type{tt})(:,kk) = temp.avg_rot; 
    end
    % Average across all flies
    temp.velmean = mean(GD_velocity(cond).(Type{tt}),2);
    GD_velocity(cond).([Type{tt} '_avg']) = temp.velmean;
    temp.rotmean = mean(GD_rotvelocity(cond).(Type{tt}),2);
    GD_rotvelocity(cond).([Type{tt} '_avg']) = temp.rotmean;
    % Error across all flies
    temp.velerr = std(GD_velocity(cond).(Type{tt}),0,2)/sqrt(num.fly);
    GD_velocity(cond).([Type{tt} '_err']) = temp.velerr;
    temp.roterr = std(GD_rotvelocity(cond).(Type{tt}),0,2)/sqrt(num.fly);
    GD_rotvelocity(cond).([Type{tt} '_err']) = temp.roterr;
  end
end


% figure! SPEED FOR ALL STIMULI
fig = figure; set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900])
    for cond = 1:num.conds
        subplot(4,7,cond)
        switch cond
            case num2cell(1:7)
                COND = 1;
            case num2cell(8:14)
                COND = 8;
            case num2cell(15:21)
                COND = 15;
            case num2cell(22:28)
                COND = 22;
        end
            temp.x = -num.stim_length:num.stim_length;
            temp.y = [GD_velocity(COND).Control_avg(1:end-1); GD_velocity(COND).Stim_avg];
            temp.err = [GD_velocity(COND).Control_err(1:end-1); GD_velocity(COND).Stim_err];
            hold all
            fill_data = error_fill(temp.x, temp.y, temp.err);
                h = fill(fill_data.X, fill_data.Y, get_color('grey'), 'EdgeColor','none');
                set(h, 'facealpha', 0.2)
            hh = plot(temp.x, temp.y, 'k', 'LineWidth', 2);
%          cond = 3;
            temp.x = -num.stim_length:num.stim_length;
            temp.y = [GD_velocity(cond).Control_avg(1:end-1); GD_velocity(cond).Stim_avg];
            temp.err = [GD_velocity(cond).Control_err(1:end-1); GD_velocity(cond).Stim_err];
            fill_data = error_fill(temp.x, temp.y, temp.err);
                h = fill(fill_data.X, fill_data.Y, get_color('lightgreen'), 'EdgeColor','none');
                set(h, 'facealpha', 0.2)
            hh = plot(temp.x, temp.y, 'color', Color('ForestGreen'), 'LineWidth', 2);

        vline([0, parameters.conds_matrix(cond).opto*num.fps], 'k-')
    end
       
subplot(4,7,1)
    ylabel({'Walking CW'; 'Speed'})
    title({structure_name; 'SPEED'})
subplot(4,7,8)
    ylabel({'Walking CCW'; 'Speed'})
subplot(4,7,15)
    ylabel({'Turning CW'; 'Speed'})
subplot(4,7,22)
    ylabel({'Turning CCW'; 'Speed'})

save_figure(fig, [figures_dir structure_name ' Speed All Conditions']);

% figure! SPEED FOR ALL STIMULI
fig = figure; set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900])
    for cond = 1:num.conds
        subplot(4,7,cond)
        switch cond
            case num2cell(1:7)
                COND = 1;
            case num2cell(8:14)
                COND = 8;
            case num2cell(15:21)
                COND = 15;
            case num2cell(22:28)
                COND = 22;
        end
            temp.x = -num.stim_length:num.stim_length;
            temp.y = [GD_rotvelocity(COND).Control_avg(1:end-1); GD_rotvelocity(COND).Stim_avg];
            temp.err = [GD_rotvelocity(COND).Control_err(1:end-1); GD_rotvelocity(COND).Stim_err];
            hold all
            fill_data = error_fill(temp.x, temp.y, temp.err);
            h = fill(fill_data.X, fill_data.Y, get_color('grey'), 'EdgeColor','none');
                set(h, 'facealpha', 0.2)
            hh = plot(temp.x, temp.y, 'k', 'LineWidth', 2);
%          cond = 3;
            temp.x = -num.stim_length:num.stim_length;
            temp.y = [GD_rotvelocity(cond).Control_avg(1:end-1); GD_rotvelocity(cond).Stim_avg];
            temp.err = [GD_rotvelocity(cond).Control_err(1:end-1); GD_rotvelocity(cond).Stim_err];
            fill_data = error_fill(temp.x, temp.y, temp.err);
            h = fill(fill_data.X, fill_data.Y, get_color('lightgreen'), 'EdgeColor','none');
                set(h, 'facealpha', 0.2)
            hh = plot(temp.x, temp.y, 'color', Color('ForestGreen'), 'LineWidth', 2);

        vline([0, parameters.conds_matrix(cond).opto*num.fps], 'k-')
    end
subplot(4,7,1)
    ylabel({'Walking CW'; 'Rot Vel'})
    title({structure_name; 'Rotational Velocity'})
subplot(4,7,8)
    ylabel({'Walking CCW'; 'Rot Vel'})
subplot(4,7,15)
    ylabel({'Turning CW'; 'Rot Vel'})
subplot(4,7,22)
    ylabel({'Turning CCW'; 'Rot Vel'})

save_figure(fig, [figures_dir, structure_name, ' Rot Velocity All Conditions']);
    
    
%% WT or no laser: RUN ALL CONDS TOGETHER
clear rotvelocity speed

for kk = 1:num.fly
    for tt = 1:2
        for cond = 1:num.conds
            raw = nanmean(fly(kk).(Type{tt}).speed(cond).data,2);
            speed.(Type{tt}).avg(cond).data(:,kk) = raw;
%             speed.(Type{tt}).avg(cond).err(:,kk) = raw;
        end
    end
end

for tt = 1:2
    clear raw
    idx = 0;
    for cond = 1:14
        idx = idx+1;
        raw = speed.(Type{tt}).avg(cond).data;
        speed.(Type{tt}).str(:,idx) = nanmean(raw,2);
    end
    speed.(Type{tt}).strerr = nanstd(speed.(Type{tt}).str,0,2)/sqrt(12);
    speed.(Type{tt}).stravg = nanmean(speed.(Type{tt}).str,2);
    
    clear raw
    idx = 0;
    for cond = 15:28
        idx = idx+1;
        raw = speed.(Type{tt}).avg(cond).data;
        speed.(Type{tt}).rot(:,idx) = nanmean(raw,2);
    end
    speed.(Type{tt}).roterr = nanstd(speed.(Type{tt}).rot,0,2)/sqrt(12);
    speed.(Type{tt}).rotavg = nanmean(speed.(Type{tt}).rot,2);
end
    
% Rotational Velocity

for kk = 1:num.fly
    for tt = 1:2
        for cond = 1:num.conds
            raw = nanmean(fly(kk).(Type{tt}).rotvelocity(cond).data,2);
            rotvelocity.(Type{tt}).avg(cond).data(:,kk) = raw;
        end
    end
end

for tt = 1:2
    clear raw
    idx = 0;
    for cond = 1:14
        idx = idx+1;
        raw = rotvelocity.(Type{tt}).avg(cond).data;
        
        if cond >=5 && cond <=14
            rotvelocity.(Type{tt}).str(:,idx) = nanmean(raw,2);
        else
            rotvelocity.(Type{tt}).str(:,idx) = nanmean(raw,2)*-1;
        end
       
    end
    rotvelocity.(Type{tt}).strerr = nanstd(rotvelocity.(Type{tt}).str,0,2)/sqrt(12);
    rotvelocity.(Type{tt}).stravg = nanmean(rotvelocity.(Type{tt}).str,2);
    
    clear raw
    idx = 0;
    for cond = 15:28
        idx = idx+1;
        raw = rotvelocity.(Type{tt}).avg(cond).data;
        
        if cond >=22 && cond <=28
            rotvelocity.(Type{tt}).rot(:,idx) = nanmean(raw,2);
        else
            rotvelocity.(Type{tt}).rot(:,idx) = nanmean(raw,2)*-1;
        end
    end
    rotvelocity.(Type{tt}).roterr = nanstd(rotvelocity.(Type{tt}).rot,0,2)/sqrt(12);
    rotvelocity.(Type{tt}).rotavg = nanmean(rotvelocity.(Type{tt}).rot,2);
end
       

fig = figure; set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900])
subplot(2,2,1) %str stim
hold all
% stim with light
    temp.x = (-60:60);
    temp.Y1 = speed.Control.stravg;
    temp.Y2 = speed.Stim.stravg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [speed.Control.strerr(1:end-1); speed.Stim.strerr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)
  
vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Speed (cm/s)')
title('Speed during Walking')

subplot(2,2,2) %str stim
hold all
% stim with light
    temp.Y1 = speed.Control.rotavg;
    temp.Y2 = speed.Stim.rotavg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [speed.Control.roterr(1:end-1); speed.Stim.roterr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)

vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Speed (cm/s)')
title('Speed during Turning')

% save_figure(fig, [figures_dir, filename(1:end-4), ' Walking and Turning Nof3'])

% fig = figure; set(fig, 'color', 'w', 'pos', [50, 50, 1450, 900])
subplot(2,2,3) %str stim
hold all
% stim with light
    temp.x = (-60:60);
    temp.Y1 = rotvelocity.Control.stravg;
    temp.Y2 = rotvelocity.Stim.stravg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [rotvelocity.Control.strerr(1:end-1); rotvelocity.Stim.strerr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)
  
vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Rot Vel (deg/s)')
title('Rotational Velocity during Walking')

subplot(2,2,4) %str stim
hold all
% stim with light
    temp.Y1 = rotvelocity.Control.rotavg;
    temp.Y2 = rotvelocity.Stim.rotavg;
    temp.Y = [temp.Y1(1:end-1); temp.Y2];
    temp.err = [rotvelocity.Control.roterr(1:end-1); rotvelocity.Stim.roterr];
    fill_data = error_fill(temp.x, temp.Y, temp.err);
    h = fill(fill_data.X, fill_data.Y, get_color('black'), 'EdgeColor','none');
        set(h, 'facealpha', 0.2)
    plot(temp.x, temp.Y, 'color', 'k', 'linewidth', 3)

vline([0,0.09*num.fps], 'g-')
xlabel('Time (fps)')
ylabel('Rot velocity (deg/s)')
title({'Rotational Velocity during Turning'; [fly(1).param.cross ', N = ' num2str(num.fly)]})

% save_figure(fig, [figures_dir, filename(1:end-4), ' Walking and Turning Nof3'])

save_figure(fig, [figures_dir, structure_name, ' Walking and Turning'])


toc


