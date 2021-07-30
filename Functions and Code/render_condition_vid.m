
function render_condition_vid(fly, cond, rep)
% render_conditon_vid(fly, cond, rep)
% Create a video that plots the speed, rot velocity, flatpath, and 
% arena pattern during a specific condition for both stim and control
% Inputs:
% 'fly' [structure with all_data]
% 'cond' [condition number to plot]
% 'rep' [rep number to plot]
% Outputs:
% video saved to the fly's analysis folder
%     
% ES Dickinson, University of Washington, Dec 2018


%load variables
[num, ~, ~, labels] = NUM(fly.parameters);
Type = {'Stim', 'Control'};

% variables:
x = (-num.OL:1/num.fps:num.OL);
num.frames = length(x);
% ----------ARENA INFORMATION------------%
% set arena information
for param = 1:2 %stim|control   
	arena.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.arena_X);
	laser.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.laser_trig);
	basler.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.camA_trig);
end
arena.all = [arena.Control(1:end-1); arena.Stim];
laser.all = [laser.Control(1:end-1); laser.Stim];
basler.all = [basler.Control(1:end-1); basler.Stim];
% ----------SPEED & ROTVELOCITY INFORMATION------------%
for param = 1:2 %stim|control   
    speed.(Type{param}) = fly.(Type{param}).speed(cond).data(:,rep);
    rotvelocity.(Type{param}) = fly.(Type{param}).rotvelocity(cond).data(:,rep);
end
speed.all = [speed.Control(1:end-1); speed.Stim];
speed_max = max(speed.all)+1;
rotvelocity.all = [rotvelocity.Control(1:end-1); rotvelocity.Stim];
rotvelocity.max = max(rotvelocity.all)+1;
rotvelocity.min = min(rotvelocity.all)-25;
% ----------ROTATION VECTORS INFORMATION------------%
for param = 1:2 %stim|control   
  rot_X.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.X);
  rot_Y.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.Y);
  rot_Z.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.Z);
end
rot_X.all = rad2deg([rot_X.Control(1:end-1); rot_X.Stim])*num.fps;
rot_Y.all = rad2deg([rot_Y.Control(1:end-1); rot_Y.Stim])*num.fps;
rot_Z.all = rad2deg([rot_Z.Control(1:end-1); rot_Z.Stim])*num.fps;
rot_X.max = max([rot_X.all; rot_Y.all; rot_Z.all]) + 50;
rot_X.min = min([rot_X.all; rot_Y.all; rot_Z.all]) - 50;
% ----------TRAJECTORY INFORMATION------------%
for param = 1:2 %stim|control
% Position data     
  int_X.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.int_x)*num.ball_radius;
  int_Y.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.int_y)*num.ball_radius;
  %total distance information
  int_X.distance.(Type{param}) = abs(diff(int_X.(Type{param})));
  int_Y.distance.(Type{param}) = abs(diff(int_Y.(Type{param})));
  total_distance.(Type{param}) = sum(sqrt(int_X.distance.(Type{param}).^2 + int_Y.distance.(Type{param}).^2));
end
int_X.all = [int_X.Control(1:end-1); int_X.Stim];
int_Y.all = [int_Y.Control(1:end-1); int_Y.Stim];
% Laser data
opto_off = round(fly.parameters.conds_matrix(cond).opto*num.fps); %frame where laser turned off
int_X.laser = fly.Stim.data(cond,rep).raw(1:opto_off,labels.int_x)*num.ball_radius;
int_Y.laser = fly.Stim.data(cond,rep).raw(1:opto_off,labels.int_y)*num.ball_radius;
int_X.max = max(int_X.all)+0.3;
int_X.min = min(int_X.all)-0.3;
int_Y.max = max(int_Y.all)+0.3;
int_Y.min = min(int_Y.all)-0.3;


% check/create folder for condition videos in data analysis within raw data
analysis_directory = [fly.parameters.save_location_matlab 'Analysis\'];
if ~isfolder(analysis_directory)
    mkdir(analysis_directory)
    fprintf('\n Created ''Analysis'' folder \n')
end
video_name = [fly.parameters.matlab_data_file ' Cond Plot Vid R' num2str(rep) 'C' num2str(cond)];

%% MAKE VIDEO: 
fig = figure; set(fig, 'pos',[5 10 1900 980], 'color','w');
hold on
v = VideoWriter([analysis_directory video_name '.avi'], 'Uncompressed AVI');
v.FrameRate = 10;
open(v);
for datapoint = 1:num.frames
set(fig, 'color','w');
%-------------------------------------------------------------------------%
%                           Figure Information
%-------------------------------------------------------------------------%
% SPEED SUBPLOT
subplot(2,3,1)
hold all
%plot arena info
xlim([-num.OL, num.OL])
    yyaxis right
    ylim([-20,10])
    plot(x(1:datapoint), arena.all(1:datapoint), 'b:');
    plot(x(1:datapoint), laser.all(1:datapoint), 'g:');
    plot(x(1:datapoint), basler.all(1:datapoint), 'r:');
    ylabel('Arena Parameters')
%plot Fictrac data
    yyaxis left
    ylim([0, speed_max])
    plot(x(1:datapoint), speed.all(1:datapoint), 'k', 'LineWidth', width)
    % add lines for laser
    vline(0, 'g-')
    vline(fly.parameters.conds_matrix(cond).opto, 'g-')
    % add lines camera times
    vline((0-fly.parameters.basler_delay), 'r-')
    vline((fly.parameters.OL_time-fly.parameters.basler_delay), 'r-')
    hline(fly.parameters.min_speed, 'k:')
    %labels
    title({['Speed | ' fly.parameters.condition_label{cond}],...
          ['rep ' num2str(rep) ' | cond ' num2str(cond)]})
    xlabel('Time (sec)')
    ylabel('Speed (cm/sec)')
hold off

%-------------------------------------------------------------------------%

% ROTATIONAL VELOCITY SUBPLOT
subplot(2,3,4)
hold all
%plot arena info
xlim([-num.OL, num.OL])
    yyaxis right  
    ylim([-20,10])
    plot(x(1:datapoint), arena.all(1:datapoint), 'b:');
    plot(x(1:datapoint), laser.all(1:datapoint), 'g:');
    plot(x(1:datapoint), basler.all(1:datapoint), 'r:');
    ylabel('Arena Parameters')
%plot Fictrac data
    yyaxis left
    ylim([rotvelocity.min, rotvelocity.max])
    plot(x(1:datapoint), rotvelocity.all(1:datapoint), 'k', 'LineWidth', width)
% add lines for laser
    vline(0, 'g-')
    vline(fly.parameters.conds_matrix(cond).opto, 'g-')
    % add lines camera times
    vline((0-fly.parameters.basler_delay), 'r-')
    vline((fly.parameters.OL_time-fly.parameters.basler_delay), 'r-')
    hline(0, 'k:')
    %labels
    title({['Rot Velocity | ' fly.parameters.condition_label{cond}],...
          ['rep ' num2str(rep) ' | cond ' num2str(cond)]})
    xlabel('Time (sec)')
    ylabel('Rotational Velocity (deg/sec)')
hold off

%-------------------------------------------------------------------------%

% DELTA ROTATION VECTORS SUBPLOT
subplot(2,3,2)
hold all
    ylim([rot_X.min, rot_X.max])
    xlim([-num.OL, num.OL])
    % plot data
    plot(x(1:datapoint), rot_X.all(1:datapoint), 'r', 'LineWidth', width);
    plot(x(1:datapoint), rot_Y.all(1:datapoint), 'b', 'LineWidth', width);
    plot(x(1:datapoint), rot_Z.all(1:datapoint), 'g', 'LineWidth', width);
    % add lines for laser
    vline(0, 'g')
    vline(fly.parameters.conds_matrix(cond).opto, 'g')
    % add lines camera times
    vline((0-fly.parameters.basler_delay), 'r')
    vline((fly.parameters.OL_time-fly.parameters.basler_delay), 'r')
    legend({'roll | tilt', 'pitch | forward', 'yaw | heading'}, 'Location', 'southeast')
    legend('boxoff')
    ylabel('Change in roll, pitch, & yaw (deg/sec)')
    xlabel('Time (sec)')
    title('Delta rotation vector')
hold off

%-------------------------------------------------------------------------%

% TRAJECTORY SUBPLOT
subplot(2,3,5)
hold all
    xlim([int_X.min, int_X.max])
    ylim([int_Y.min, int_Y.max])
    % plot raw data
    switch datapoint
        case num2cell(1:num.stim_length) %Control period; plot only control
            scatter(int_X.Control(1:datapoint), int_Y.Control(1:datapoint),...
                    20, Color('grey'), 'filled') 
        case num2cell(num.stim_length:num.stim_length+opto_off) %Stim period; plot control&stim%laser
            if opto_off == 0 %for the control laser condition
                scatter(int_X.Control, int_Y.Control, 20, Color('grey'), 'filled') 
                scatter(int_X.Stim(1), int_Y.Stim(1), 20, Color('blue'), 'filled') 
            else %conditions with laser on
                scatter(int_X.Control, int_Y.Control, 20, Color('grey'), 'filled') 
                scatter(int_X.Stim(1:datapoint-num.stim_length), int_Y.Stim(1:datapoint-num.stim_length),...
                    20, Color('blue'), 'filled') 
                scatter(int_X.Stim(1:datapoint-num.stim_length), int_Y.Stim(1:datapoint-num.stim_length),...
                    28, 'g')
            end 
        case num2cell(num.stim_length+opto_off+1:num.frames) %Post-stim data; plot control
            scatter(int_X.Control, int_Y.Control, 20, Color('grey'), 'filled') 
            scatter(int_X.Stim(1:datapoint-num.stim_length), int_Y.Stim(1:datapoint-num.stim_length),...
                    20, Color('blue'), 'filled') 
            scatter(int_X.laser, int_Y.laser, 28, 'g')     
    end
    scatter(int_X.all(datapoint), int_Y.all(datapoint), 10, 'k', 'filled')
    legend({'Stim', 'Control', 'laser'}, 'Location', 'southeast')
    legend('boxoff')
    xlabel('x postion (cm)')
    ylabel('y position (cm)') 
    title({['Stim distance: ' num2str(total_distance.Stim) ' cm']...
           ['Control distance: ' num2str(total_distance.Control) ' cm']})
hold off

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Collect image info for movie file
f = getframe(fig);
writeVideo(v, f)
clf('reset')
clear y r s z a
end

close(v)
close all
fprintf(['\n Saved ' video_name '\n'])
end
