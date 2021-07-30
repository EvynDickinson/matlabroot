
% function combined_vid_n_plots(fly, cond, rep)
% 
% filename = ['G:\My Drive\Data\FicTrac Raw Data\11.16.18\Fly 1_0\pose-2d-filtered\11162018_fly1_0 R1C1 Cam-A str-cw-0 sec.h5'];
% info = hdf5info(filename);
% name = info.GroupHierarchy.Groups.Datasets.Name;
% data = h5read(filename,name);


%load variables
[num, ~, ~, labels] = NUM(fly.parameters);
Type = {'Stim', 'Control'};  
% variables:
x = (-num.OL:1/num.fps:num.OL);
num.frames = length(x);
width = 2.5;
n = 10; %upsampling interp value
% warning('off', 'MSID')

% % Raw Videos not processed videos:
% % get grayscale values for the specific fly walking video
% video_path = [fly.parameters.save_location_matlab 'Raw Video\'];
% video_name = [fly.parameters.matlab_data_file ' R' num2str(rep) 'C' num2str(cond) ' Cam-A '...
%               fly.parameters.conds_matrix(cond).label '.avi'];  
% [~, pret_images_A] = convert_vid_to_pixel_values([video_path, video_name], fly.parameters);
% video_name = [fly.parameters.matlab_data_file ' R' num2str(rep) 'C' num2str(cond) ' Cam-C '...
%               fly.parameters.conds_matrix(cond).label '.avi'];  
% [~, pret_images_C] = convert_vid_to_pixel_values([video_path, video_name], fly.parameters);

% % Processed Video Information: 
% % videoC = 3D model of fly legs
video_name = [fly.parameters.matlab_data_file ' R' num2str(rep) 'C' num2str(cond) ' Cam-C '...
              fly.parameters.conds_matrix(cond).label '.avi'];  %
video_path = [fly.parameters.vid_folder 'videos-labeled\'];
[~, pret_images_C] = convert_vid_to_pixel_values([video_path, video_name], fly.parameters, 1);

% get grayscale|color values for the specific fly walking video
video_path = fly.parameters.vid_3d;
% % videoA = labeled fly vid
video_name = [fly.parameters.matlab_data_file ' R' num2str(rep) 'C' num2str(cond) '  '...
              fly.parameters.conds_matrix(cond).label '.avi'];  %
[~, pret_images_A] = convert_vid_to_pixel_values([video_path, video_name], fly.parameters, 1);
         


%% GATHER INFORMATION
% ----------ARENA INFORMATION------------%
% set arena information
for param = 1:2 %stim|control   
    arena.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.arena_X);
    laser.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.laser_trig);
    basler.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.camA_trig);
end
arena.all = interp([arena.Control(1:end-1); arena.Stim], n);
laser.all = interp([laser.Control(1:end-1); laser.Stim], n);
basler.all = interp([basler.Control(1:end-1); basler.Stim], n);
% ----------SPEED & ROTVELOCITY INFORMATION------------%
for param = 1:2 %stim|control   
    speed.(Type{param}) = fly.(Type{param}).speed(cond).data(:,rep);
    rotvelocity.(Type{param}) = fly.(Type{param}).rotvelocity(cond).data(:,rep);
end
speed.all = interp([speed.Control(1:end-1); speed.Stim], n);
speed_max = max(speed.all)+1;
rotvelocity.all = interp([rotvelocity.Control(1:end-1); rotvelocity.Stim], n);
rotvelocity.max = max(rotvelocity.all)+1;
rotvelocity.min = min(rotvelocity.all)-25;
% ----------ROTATION VECTORS INFORMATION------------%
for param = 1:2 %stim|control   
  rot_X.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.X);
  rot_Y.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.Y);
  rot_Z.(Type{param}) = fly.(Type{param}).data(cond,rep).raw(:,labels.Z);
end
rot_X.all = interp(rad2deg([rot_X.Control(1:end-1); rot_X.Stim])*num.fps, n);
rot_Y.all = interp(rad2deg([rot_Y.Control(1:end-1); rot_Y.Stim])*num.fps, n);
rot_Z.all = interp(rad2deg([rot_Z.Control(1:end-1); rot_Z.Stim])*num.fps, n);
rot_X.max = max([rot_X.all; rot_Y.all; rot_Z.all]) + 50;
rot_X.min = min([rot_X.all; rot_Y.all; rot_Z.all]) - 50;
% ----------TRAJECTORY INFORMATION------------%
for param = 1:2 %stim|control
% Position data     
  int_X.(Type{param}) = interp(fly.(Type{param}).data(cond,rep).raw(:,labels.int_x)*num.ball_radius, n);
  int_Y.(Type{param}) = interp(fly.(Type{param}).data(cond,rep).raw(:,labels.int_y)*num.ball_radius, n);
  %total distance information
  int_X.distance.(Type{param}) = abs(diff(int_X.(Type{param})));
  int_Y.distance.(Type{param}) = abs(diff(int_Y.(Type{param})));
  total_distance.(Type{param}) = sum(sqrt(int_X.distance.(Type{param}).^2 + int_Y.distance.(Type{param}).^2));
end
int_X.all = ([int_X.Control(1:end-n); int_X.Stim]);
int_Y.all = ([int_Y.Control(1:end-n); int_Y.Stim]);
% Laser data
opto_off = round(fly.parameters.conds_matrix(cond).opto*num.fps); %frame where laser turned off

int_X.max = max(int_X.all)+0.3;
int_X.min = min(int_X.all)-0.3;
int_Y.max = max(int_Y.all)+0.3;
int_Y.min = min(int_Y.all)-0.3;

% ------------lEG ANGLES-------------%
%load the leg angles:
excelname = [fly.parameters.vid_dir, 'summaries\angles.csv'];
% [~, ~, excelfile] = xlsread(excelname); %OR
load('summaries') %loads excelfile from previously stored version
% Excel column for the headers:
angle_names = {'L1_CF', 'L1_FTi', 'L1_TiTa'};
% angle_names = {'L1_FTi', 'L2_FTi', 'L3_FTi'};
num.angles = length(angle_names);
Excel.headers = excelfile(1,:); %xltitles
Excel.date = find(strcmpi('folder_1',Excel.headers) == 1);
Excel.flynum = find(strcmpi('folder_2',Excel.headers) == 1);
Excel.filename = find(strcmpi('filename',Excel.headers) == 1);
row_ind = strcmpi(excelfile(:,Excel.filename), video_name(1:end-4)); %row locations index
if sum(row_ind) == 0
    error('Warning: No data found in the excel file')
end
% Extract the data: 
for ii = 1:num.angles
    Excel.(angle_names{ii}) = find(strcmpi((angle_names{ii}),Excel.headers) == 1);
    angles.(angle_names{ii}) = cell2mat(excelfile(row_ind, Excel.(angle_names{ii})));
    a = angle_names{ii};
    [b, c] = strtok( a, '_' );
    c = strtok( c, '_' );
    
    angles.names{ii} = [b, '-', c]; clear a b c
end

% 
% figure;
% hold all
% strt_buffer = (fly.parameters.OL_time-fly.parameters.basler_delay)*num.fps*n;
% end_buffer = (fly.parameters.OL_time-fly.parameters.basler_length+fly.parameters.basler_delay)*(num.fps*n);
% for ii = 1:num.angles
%     angles.([angle_names{ii} '_Y']) = [nan(strt_buffer,1); angles.(angle_names{ii}); nan(end_buffer,1)];
%     angles.size = length(angles.([angle_names{ii} '_Y']));
%     intrvl = abs((x(end)-x(1))/angles.size);
%     angles.X = x(1):intrvl:x(end);
%     plot(angles.X(1:angles.size), angles.([angle_names{ii} '_Y']), 'LineWidth', 3)
% end
% xlim([-2, 2])
% xlabel('Time (sec)')
% ylabel('Joint Angle (deg)')
% legend(angles.names, 'Location', 'NorthWest')
% vline(0, 'k')
% vline(fly.parameters.conds_matrix(cond).opto)

% **** WORKING HERE******** %



%%
% fprintf(['\n Cond: ' num2str(cond)])
tic

x = interp(-num.OL:1/num.fps:num.OL, n);  
num.frames = length((-num.OL:1/num.fps:num.OL));
analysis_directory = filepath;
opto_on = (2*num.fps*n)-1;
off_opto = (2*num.fps*n)+(fly.parameters.conds_matrix(cond).opto*n*num.fps);


% MAKE VIDEO

fig = figure; set(fig, 'pos',[15 50 1900 940], 'color','w'); 
% set(0,'DefaultFigureVisible','off'); % Don't display figure
hold on
% v = VideoWriter([analysis_directory video_name '.avi'], 'Uncompressed AVI');
v = VideoWriter([analysis_directory video_name '.avi'], 'Uncompressed AVI');

v.FrameRate = 100;
open(v);
for datapoint = 1:1200
set(fig, 'color','w');
%-------------------------------------------------------------------------%
%                           Figure Information
%-------------------------------------------------------------------------%
% 
if datapoint >= opto_on && datapoint <= off_opto 
    light_on = true;
else
    light_on = false; 
end
%-------------------------------------------------------------------------%
% SPEED SUBPLOT
subplot(2,3,1)
hold all
%plot arena info
xlim([-num.OL, num.OL])
    yyaxis right
    ylim([-20,10])
    plot(x(1:datapoint), (arena.all(1:datapoint)), 'b-');
    plot(x(1:datapoint), (laser.all(1:datapoint)), 'g-');
    plot(x(1:datapoint), (basler.all(1:datapoint)), 'r-');
    ylabel('Arena Parameters')
%plot Fictrac data
    yyaxis left
    ylim([0, speed_max])
    plot(x(1:datapoint), (speed.all(1:datapoint)), 'k', 'LineWidth', width)
    % add lines for laser
    vline(0, 'g-')
    vline(fly.parameters.conds_matrix(cond).opto, 'g-')
    % add lines camera times
    vline((0-fly.parameters.basler_delay), 'r-')
    vline((fly.parameters.OL_time-fly.parameters.basler_delay), 'r-')
    hline(fly.parameters.min_speed, 'k:')
    %labels
%     title({['Speed | ' fly.parameters.condition_label{cond}],...
%           ['rep ' num2str(rep) ' | cond ' num2str(cond)]})
%     title({['Speed | ' fly.parameters.condition_label{cond}]})
    title('Speed')
    xlabel('Time (sec)')
    ylabel('Speed (cm/sec)')
hold off

%-------------------------------------------------------------------------%

% ROTATIONAL VELOCITY SUBPLOT
subplot(2,3,4)
hold all
xlim([-num.OL, num.OL])
%plot Fictrac data
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
    title({'Rot Velocity'}) % ['rep ' num2str(rep) ' | cond ' num2str(cond)]
    xlabel('Time (sec)')
    ylabel('Rotational Velocity (deg/sec)')
hold off

%-------------------------------------------------------------------------%
% LEG ANGLES SUBPLOT
subplot(2,3,5)
    color_matrix = {'green', 'purple', 'pink'};
    hold all
    strt_buffer = (fly.parameters.OL_time-fly.parameters.basler_delay)*num.fps*n;
    end_buffer = (fly.parameters.OL_time-fly.parameters.basler_length+fly.parameters.basler_delay)*(num.fps*n);
    %find angle max and min
    for ii = 1:3
        angles.allmax(ii) = nanmax(angles.(angle_names{ii}));
        angles.allmin(ii) = nanmin(angles.(angle_names{ii}));
    end
    angles.max = max(angles.allmax)+40;
    angles.min = min(angles.allmin);
    
    for ii = 1:2
        angles.([angle_names{ii} '_Y']) = [nan(strt_buffer,1); angles.(angle_names{ii}); nan(end_buffer,1)];
        angles.size = length(angles.([angle_names{ii} '_Y']));
        intrvl = abs((x(end)-x(1))/angles.size);
        angles.X = x(1):intrvl:x(end);
        h = plot(angles.X(1:datapoint), angles.([angle_names{ii} '_Y'])(1:datapoint));
        set(h, 'LineWidth', 1.5, 'color', get_color(color_matrix{ii}))
    end
    ylim([0, angles.max])
    xlim([-0.75, 2])
    xlabel('Time (sec)')
    ylabel('Joint Angle (deg)')
    legend(angles.names(1:2), 'Location', 'NorthEast')
    title('Joint Angles')
    vline(0, 'g-')
    % add lines camera times
    vline((0-fly.parameters.basler_delay), 'r-')
    vline((fly.parameters.OL_time-fly.parameters.basler_delay), 'r-')
    vline(fly.parameters.conds_matrix(cond).opto, 'g-')
%-------------------------------------------------------------------------%

% TRAJECTORY SUBPLOT
subplot(2,3,2)
hold all
    xlim([int_X.min, int_X.max])
    ylim([int_Y.min, int_Y.max])
    % plot raw data    
    control_start = 1;
    control_end = num.stim_length*n;
    %case 2 timing points:
    stim_start = num.stim_length*n;
    stim_end = (num.stim_length+opto_off)*n;
    %case 3 timing points:
    post_start = (num.stim_length+opto_off)*n;
    post_end = num.frames*n;
    switch datapoint
        case num2cell(control_start:control_end) %Control period; plot only control
            scatter(int_X.Control(1:datapoint), int_Y.Control(1:datapoint),...
                    20, Color('grey'), 'filled') 
        case num2cell(stim_start:stim_end) %Stim period; plot control&stim%laser
            if opto_off == 0 %for the control laser condition
                scatter(int_X.Control, int_Y.Control, 20, Color('grey'), 'filled') 
                scatter(int_X.Stim(1), int_Y.Stim(1), 20, Color('blue'), 'filled') 
            else %conditions with laser on
                scatter(int_X.Control, int_Y.Control, 20, Color('grey'), 'filled') 
                scatter(int_X.Stim(1:datapoint-(num.stim_length*n)), int_Y.Stim(1:datapoint-(num.stim_length*n)),...
                    20, Color('blue'), 'filled') 
                scatter(int_X.Stim(1:datapoint-(num.stim_length*n)), int_Y.Stim(1:datapoint-(num.stim_length*n)),...
                    35, 'g')
            end 
        case num2cell(post_start:post_end) %Post-stim data; plot control
            scatter(int_X.Control, int_Y.Control, 20, Color('grey'), 'filled') 
            scatter(int_X.Stim(1:datapoint-(num.stim_length*n)), int_Y.Stim(1:datapoint-(num.stim_length*n)),...
                    20, Color('blue'), 'filled') 
                if sum(opto_off) > 0
                    scatter(int_X.Stim(1:opto_off*n), int_Y.Stim(1:opto_off*n),...
                    28, 'g')
%                     scatter(int_X.laser, int_Y.laser, 28, 'g')     
                end
    end
    scatter(int_X.all(datapoint), int_Y.all(datapoint), 20, 'r', 'filled')
    title({'2D Fly Path'; [fly.parameters.condition_label{cond} ' activation']})
    
%         legend({'Stim', 'Control', 'laser'}, 'Location', 'southeast')
%         legend('boxoff')
%     xlabel('x postion (cm)')
%     ylabel('y position (cm)') 
%     title({['Stim distance: ' num2str(total_distance.Stim) ' cm']...
%            ['Control distance: ' num2str(total_distance.Control) ' cm']})
hold off
%-------------------------------------------------------------------------%

% LABELED VIDEO
pos = [0.64 0.53 0.4 0.4];
subplot('Position', pos) 
% subplot(2,3,3)
imshow(flip(pret_images_A(datapoint).data,2), 'border', 'tight')
title('Labeled Joints')
%-------------------------------------------------------------------------%

% 3D LEG CONSTRUCTION
pos = [0.64 0.06 0.4 0.4];
subplot('Position', pos) 
hold all
% subplot(2,3,6);
imshow(pret_images_C(datapoint).data, 'border', 'tight')
if light_on == true
    scatter(25, 25, 500, 'g', 'filled')
end
title('3D Leg Reconstruction')



%
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

toc

% end

% fig = figure; set(fig, 'pos',[5 10 1900 980], 'color','w');
% set(0,'DefaultFigureVisible','on'); %Display figure
% imshow(pret_images_A(600).data, 'border', 'tight')




