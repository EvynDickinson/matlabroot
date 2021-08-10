%   load('Current_data_file')
%   MIDSTIM_Find_conditions_and_alignment(fly) 
%   start_gui

% Export_Cross_Ns

% Next Step:
% MIDSTIM_Step_2_sort_and_condense_data
imaqreset
closepreview  
clear all; close all; clc 

param = load_parameters;
param.fly_num = '4_1';
param.LED_intensity = [0, 7, 7, 7, 7, 7, 7];
% param.LED_intensity = [0, 4, 4, 4, 4, 4, 4];
% param.LED_intensity = [0, 3, 4, 5, 6, 7, 8];  % Intensity of the laser
param.light_length =  [0.720];  % LED trigger time for optigentic flies

param.start_pause =  30;
param.num_reps =      1; 

%%  Set up experiment param and conditions:
%Initialize DAQ & zero condition signal
Panel_com('set_ao', [3 0]); 

% Create directories and blank data files
[blank_data, param] = create_dat_dir_n_file(param);

% Generate the conditions matrix
param.conds_matrix = generate_conditions_matrix(param);
fprintf('\n Loading arena configuration \n ... \n ... \n ...\n')

%% Create Sessions  
% LASER|Basler session
s_vid_light = daq.createSession('ni');
s_vid_light.Rate = param.basler_samplerate;
% add analog output channels for LASER|Basler
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao0', 'Voltage'); %LASER output
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao1', 'Voltage'); %Basler output

% Create sessions
% FlyCap camera trigger (tracking camera):
s = daq.createSession('ni');
s.IsContinuous = true;
param.cam_trig_fs = 4*param.fps;
param.acq_in_fs = 1000; %flycap sampling rate
s.Rate = param.acq_in_fs;

% Analog-In channels
aa = addAnalogInputChannel(s,'Dev1', (0:param.num_inchans-1), 'Voltage');

% Channel Outputs:
% counter for triggering FlyCap
ch = addCounterOutputChannel(s,'Dev1', 0, 'PulseGeneration'); %Flycap camera trigger
ch.Frequency = param.fps;

% input channels    
for i = 1:length(aa)
    aa(i).TerminalConfig =  'SingleEnded'; % 'Differential'; 
end
aa(5).TerminalConfig = 'Differential'; % change the condition signal to differential
aa(6).TerminalConfig = 'Differential'; % change the laser input signal to differential
aa(7).TerminalConfig = 'Differential'; % change the side basler Signal to differential
aa(8).TerminalConfig = 'Differential'; % change the top basler Signal to differential


% Setup source and side_vid objects for Basler Camera
for i = 1:param.num_cams
%     [param.([Alphabet(i) '_vid']), param.([Alphabet(i) '_Basler_src'])] = ...
%         initiate_camera(param.(['Cam' Alphabet(i)]));
    param.([Alphabet(i) '_vid'])= ...
        initiate_camera(param.(['Cam' Alphabet(i)]));
end

% Set up LED and Basler data queues -- ONE/STIM TYPE
[param.conds_matrix, param.LED] = create_laser_signal(param.conds_matrix, param);    

%% Plot data as it is acquired and log data to file
fig = figure; set(fig, 'pos', [10, 480, 600, 500], 'color', 'w')
lh = addlistener(s,'DataAvailable', @plotData); 
lh2 = addlistener(s,'DataAvailable',@(src, event)logData(src, event, blank_data));
micropause = 0.001;
% Send one frame to pt.grey cam to start Fictrac loading process
Panel_com('set_ao', [3 0]); 
fprintf('\n Sent one frame to FicTrac... \n');
aoMax = 32767/10*9;

startBackground(s); 
s.stop 

Panel_com('set_ao', [3 aoMax]);  %to check timing of signal out
pause(17) %wait for FicTrac to catch up

% Start the daq acquisition
Panel_com('set_ao', [3 0]);
Panel_com('all_off'); 
fprintf('Acquisition starting...\n');
startBackground(s); 

pause(param.start_pause-30); % Pause before starting stimulus procedure 
beep
pause(30);
param.target_loop_time = ((param.OL_time*2)+param.control_time)*3.5;

% Start arena experiment
pause(micropause)
Panel_com('set_ao', [3 0]);
pause(micropause)
ind = 0;
temp.vid_strt = [param.Basler_folder_name '\' param.matlab_data_file  ' R'];
for ii = 1:param.num_reps
   %set condition randomization for this rep
    conds = randperm(param.num_conds);
    param.conditions_rand(:,ii) = conds;
% ROTATE THROUGH ALL CONDITIONS:
  for jj = 1:length(conds) 
      tick_total = tic;
      tic
      ind = ind+1;
    flushdata(param.B_vid); flushdata(param.C_vid); flushdata(param.A_vid) 
    CC = conds(jj);                       % unrandomizes condition - actual condition 'number'   
    pat = param.conds_matrix(CC).pat;     % pattern ID for this condition
    Dir = param.conds_matrix(CC).dir;     % direction for the pattern 
    
    %queue the basler and LED:
    queueOutputData(s_vid_light, param.conds_matrix(CC).LED)
    
    % Prep the Basler cameras
    temp.vid_end  = [param.conds_matrix(CC).label '.avi'];
    for icam = 1:param.num_cams
        video_name = [temp.vid_strt num2str(ii) 'C' num2str(CC) ' Cam-' Alphabet(icam) ' ' temp.vid_end];
        diskLogger = VideoWriter(video_name, 'Grayscale AVI');
        diskLogger.FrameRate = param.Basler_fps; 
        param.([Alphabet(icam) '_vid']).DiskLogger = diskLogger;  
    end

    %print rep and cond and opto info
    fprintf(['\nrep = ' num2str(ii) '|num = ' num2str(jj) ':\n']);
    fprintf(['Stim ' num2str(CC) ' | ' param.conds_matrix(CC).label  '\n']);
  %-------------CONTROL---------------%
  
  % Stationary Stripe:          
  pause(micropause)
    Panel_com('set_pattern_id', param.stripe_Pat);   
    Panel_com('set_mode', [ 4 0 ]);                    
    Panel_com('set_position', [48 1]);  
    Panel_com('set_posfunc_id', [1 0]);
  pause(micropause)  
    pause(param.control_time)
    Panel_com('stop');
  %-------------STIM trial------------%
    % Cond signal:
%     Panel_com('set_ao', [3 32767/10*param.conds_matrix(CC).cond_sig]);
 pause(micropause)
  %prep basler cam
    start(param.A_vid); start(param.B_vid); start(param.C_vid);start(param.D_vid); start(param.E_vid); start(param.F_vid);
 pause(micropause)
    Panel_com('set_ao', [3 aoMax]);  
    % Set the arena stim conditions 
    Panel_com('set_pattern_id', param.conds_matrix(CC).arenapat);   
    Panel_com('set_mode', param.conds_matrix(CC).mode);                    
    Panel_com('set_position', param.conds_matrix(CC).pos);  
    Panel_com('set_posfunc_id', param.conds_matrix(CC).posfunc);
    Panel_com('send_gain_bias', param.conds_matrix(CC).gain_bias);
  pause(micropause)
    Panel_com('start'); % start the pattern   
    pause(param.OL_time-param.basler_delay);  
    startBackground(s_vid_light);           
    pause(param.OL_time)
    Panel_com('stop');
  pause(micropause)
    Panel_com('set_ao', [3 0]); % set cond signal back to zero
    %--------- extra write to disk time as needed by the program ---------%
    % Stationary Stripe:
    Panel_com('set_pattern_id', param.stripe_Pat);   
    Panel_com('set_mode', [ 4 0 ]);                    
    Panel_com('set_position', [48 1]);  
    Panel_com('set_posfunc_id', [1 0]);
  pause(micropause)
    % stops the Basler camera recording post stim
    stop(param.A_vid); stop(param.B_vid); stop(param.C_vid);stop(param.D_vid); stop(param.E_vid); stop(param.F_vid);
    actual_loop_time = toc;
    param.writetodiskpause = param.target_loop_time-actual_loop_time;
    pause(param.writetodiskpause)
    toc(tick_total)
  end %cond loop
end %rep loop

% Stop Everything
% End with Stationary Stripe:
Panel_com('set_pattern_id', param.stripe_Pat);   
Panel_com('set_mode', [ 4 0 ]);                    
Panel_com('set_position', [48 1]);  
Panel_com('set_posfunc_id', [1 0]); 
Panel_com('stop');
pause(param.end_pause);
s.stop;
delete(lh);  delete(lh2);
fclose(blank_data);
release(s)
clear s s1 

% PLOT the data
temp.fileroot = 'C:\matlabroot\data\';
d = dir(temp.fileroot);
fid2 = fopen([temp.fileroot d(length(d)).name]);
[data,count] = fread(fid2,[param.num_inchans+1,inf],'double');
fclose(fid2);

% SAVE DATA
fly.data = data;
fly.param = param; 

save('Current_data_file')  
close(fig)

%% Throw data to the analysis script:
MIDSTIM_Find_conditions_and_alignment(fly) 

%% Open camera previews to realign the next fly
answer = questdlg('Open camera gui ?', 'Notes', 'Yes', 'No', 'Yes');
switch answer
    case 'Yes'
        start_gui
    case 'No'
end

%% 

% 
% figure;hold all;
% plot(data(1,:), 'b'); %digital counter
% % plot(data(2,:), 'g'); %fly camera trigger
% % plot(data(3,:), 'r'); %LED x pos
% % plot(data(4,:), 'c'); %LED y pos
% plot(data(5,:), 'k'); %condition signal
% plot(data(6,:), 'm'); %FicTrac Heading ??????
% plot(data(7,:), 'r'); %LED trigger  
% plot(data(8,:), 'color', get_color('orange')); %Basler camera trigger
% plot(data(9,:), 'r'); %Basler camera trigger
% plot(data(10,:), 'b'); %Basler camera trigger 
% plot(data(11,:), 'r'); %Basler camera trigger
% plot(data(12,:), 'k'); %Basler camera trigger
% plot(data(13,:), 'r'); %Basler camera trigger
% 
% ylabel('volts')
% xlabel('time') 
% legend(param.chan_labels);
% % vline(peaks, 'k')

