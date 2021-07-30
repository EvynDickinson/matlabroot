
function [num, Fictrac, Matlab, labels] = NUM(parameters)
% [num, Fictrac, Matlab, labels] = NUM(parameters)
% Extract values into a structure 'num' from the parameters
% of a given experiment; loads some label variables too
% Fictrac columns, etc.
% Input: 
% 'parameters' [information on reps and conditions, etc.]
% Ouput: 
% 'num' [structure with fields of different values]
%        e.g. num.conds || num.reps
% 'Fictrac' [structure with column numbers for fictrac channels]
% 'Matlab'  [structure with column numbers/IDs
%            for matlab analogue channels]
% 'labels' ['all_data' combined column labels (for when 
%            matlab and fictrac are combined together)]
%           
%
% ES Dickinson, University of Washington, Dec 2018

% Update Variables on timing:
num.reps = parameters.num_reps;                         %repetitions
num.conds = parameters.num_conds;                       %stimulus conditions
num.OL = parameters.OL_time;                            %length of stimulus condition
num.control_time = parameters.control_time;             %length of control period
num.CL = parameters.CL_time;                            %length of closed loop
num.stim_length = parameters.fps * parameters.OL_time;             %frames per stimulus
num.control_length = parameters.fps * parameters.control_time;     %frames per CL-time
num.fps = parameters.fps;
num.ball_radius = 0.478; 
num.rad2deg = (parameters.fps*180/pi);
load('Fictrac_column_labels', 'Fictrac');
Matlab.digital_counter = 1;
Matlab.fictrac_trig = 2;
Matlab.arena_X = 3;
Matlab.arena_Y = 4;
Matlab.cond_sig = 5;
Matlab.unknown = 6;
Matlab.laser_trig = 7;
Matlab.camA_trig = 8;
Matlab.camB_trig = 9;
Matlab.camC_trig = 10;

%Fictrac labels
labels.frame_counter = 1;          %frame count of the recording
labels.X = 6;                      %X rotation of the ball=right to left
labels.Y = 7;                      %Y rotation of the ball=down side to up-side
labels.Z = 8;                      %Z rotation of the ball=forward and backward
labels.int_x = 15;                 %integrated x position of the fly (incorporates animal heading)
labels.int_y = 16;                 %integrated y position of the fly (incorporates animal heading)
labels.heading = 17;               %overall direction that the fly is heading
labels.inst_dir = 18;              %direction of the fly in that frame (in radians)
labels.speed = 19;                 %the fly's walking speed
labels.diffheading = 24;
%Matlab analogue signals
labels.digital_counter = 25;
labels.fictrac_trig = 26;
labels.arena_X = 27;
labels.arena_Y = 28;
labels.cond_sig = 29;
labels.ignore = 30;
labels.laser_trig = 31;
labels.camA_trig = 32;
labels.camB_trig = 33;
labels.camC_trig = 34;
labels.camD_trig = 35;
labels.camE_trig = 36;
labels.camF_trig = 37;

% 
% figure; hold all
% % plot(data(1,:), 'b'); %digital counter
% % plot(data(2,:), 'g'); %fly camera trigger
% % plot(data(3,:), 'r'); %LED x pos
% % plot(data(4,:), 'c'); %LED y pos
% % plot(data(5,:), 'k'); %condition signal
% % plot(data(6,:), 'm'); %FicTrac Heading
% plot(data(7,:), 'b'); %LED trigger  
% plot(data(8,:), 'r'); %Basler camera trigger
% plot(data(9,:), 'y'); %Basler camera trigger
% plot(data(10,:), 'c'); %Basler camera trigger
% ylabel('volts')
% xlabel('time')
% % legend(parameters.chan_labels);
% % vline(peaks, 'k')


% 
% min_speed = 0.5;                         %minimum speed level acceptable for 'active'
% fps = 30;                                %recorded frames per second
% ball_radius =                       %9.55mm diameter-->0.478cm radius [this varies by ball]
% stim_length = fps * parameters.OL_time;             %frames per stimulus
% control_length = fps * parameters.control_time;     %frames per CL-time
% cond_unit = stim_length + control_length;%total frames per CL and stim together (condition unit)
% Radtodeg =                  %conversion for radian to deg/sec + timing info


end