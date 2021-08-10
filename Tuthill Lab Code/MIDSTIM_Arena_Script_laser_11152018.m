
closepreview
clear all; close all; clc

% Genetic Cross Information:
Cross = {'BDP-gal4xUAS-csChrimson', 'BDP-gal4xUAS-gtACR1', 'IAV-gal4xUAS-csChrimson',...        % 1:3
         'IAV-gal4xUAS-gtACR1', 'OK371-gal4xUAS-gtACR1', 'OK371-gal4xUAS-ReachR-YFP',...        % 4:6
         'Berlin-WT', 'JR204-gal4xUAS-csChrimson', 'JR204-gal4xUAS-gtACR1',...                  % 7:9
         'JR279-gal4xUAS-csChrimson', 'JR279-gal4xUAS-gtACR1', 'JR185-gal4xUAS-csChrimson',...  % 10:12
         'JR185-gal4xUAS-gtACR1', '35C09-gal4xUAS-csChrimson', '35C09-gal4xUAS-gtACR1'...       % 13:15
         'pJFRC7-22A08-gal4xUAS-ChR88', '23C02-LexAxUAS-gtACR1', 'NorpA',...                    % 16:18 
         '81A06-gal4xcsChrimson', '81A06-gal4xcsChrimson'};  % 19:21 
     
parameters.cross = Cross{20};  


% initialize arena and experiment parameters    
FormatOut = 'mmddyyyy';
Today = datestr(datetime,FormatOut);

fly_num = '1_0';
parameters.num_reps =      1;  
parameters.vid_on =        1;                  % 1 for on, 0 for off
parameters.light_length =  [0, 0.03, 0.06, 0.09, 0.18, 0.36, 0.720];  % LED trigger time for optigentic flies
parameters.LED_light =  'green';

%start|end timing:
parameters.start_pause =   1;   % length of time to run tracking without stimulus at experiment start
parameters.end_pause =     5;     % pause after all stimulus

%light and video info
parameters.LED_intensity = 5;     % Intensity of the LED light
parameters.basler_length = 2;     % Basler trigger length in seconds
parameters.basler_delay =  0.5;   % How long before laser to start the cameras
parameters.Basler_fps =    300;   % Desired fps from the basler cam
parameters.fps =           30;    % Pt. Grey camera fps

%arena stimiulus
parameters.num.pats =  2;         % number of patterns
parameters.dir =       [1 0];     % 'direction' of rotation for pattern
parameters.speeds =    [68 -68];  % stimulus group D
parameters.width =     4;         % width of stripe in the arena (by pixel)
parameters.stripe_Pat =    2;     % wiggle pattern
parameters.grating_Pat =   3;     % dark stripe, bright background-rotating pattern

%timing info
parameters.CL_time =            3;
parameters.OL_time =            2; %time before and time after the laser onset
parameters.control_time =       3;  
parameters.writetodiskpause =   5;

parameters.num_conds = parameters.num.pats*length(parameters.light_length)*length(parameters.dir);
Panel_com('set_ao', [3 0]);   % zero condition signal and initialize daq
parameters.num_inchans = 9; 
parameters.chan_labels = {'trigger', 'x pos', 'y pos', 'cond sig', 'FicTrac heading',...
                          'LED trigger', 'Side Camera', 'Front Camera', 'Back Camera'};

% Humidity and temp:
parameters.temp = inputdlg('Temperature?'); 
parameters.humidity = inputdlg('Humidity?');
                      
%create folder for data to go into
directory = 'E:\FicTrac Raw Data\';
% directory = 'C:\matlabroot\FicTrac Raw Data\';
matlab_data_file = [Today '_fly' fly_num];
folder_date = Num2Month(date);    % today's date
if isfolder([directory folder_date]) == 1
    save_location_matlab = [directory folder_date];
else
mkdir(directory, folder_date);
    save_location_matlab = [directory folder_date];
end

if ~isdir(save_location_matlab)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', Basler_folder_name);
    uiwait(warndlg(errorMessage));
    return;
end

fprintf('\n Parameters generated \n')
fprintf('\n Loading arena configuration \n ... \n ... \n ...\n')

% Create conditions matrix|labels
Pattern = {'str', 'rot'};
Dir = {'cw', 'ccw'};
for ii = 1:length(parameters.light_length)
    Opto{ii} = [num2str(parameters.light_length(ii)) ' sec'];
end

ind = 1;
for kk = 1:parameters.num.pats
    for jj = 1:length(parameters.dir)
        for ll = 1:length(parameters.light_length)
            conds_matrix(ind).pat = kk;  
            conds_matrix(ind).dir = parameters.dir(jj);
            conds_matrix(ind).opto = parameters.light_length(ll);
            conds_matrix(ind).cond_sig = 10/(parameters.num_conds+1)*ind;
            condition_label(ind) = {[Pattern{kk} '-' Dir{jj} '-' Opto{ll}]};
            conds_matrix(ind).label = condition_label{ind};
            ind = ind+1;
        end
        % Add the parameters for the arena configurations
    end
end  
% SET ARENA PAT TO CORRESPOND TO MATRIX
for ind = 1:parameters.num_conds
%add in the arena parameters
    switch conds_matrix(ind).pat %pattern choice
        case 1 % Small wiggle: CW|CCW
                conds_matrix(ind).arenapat = 2;
                conds_matrix(ind).mode = [4 0];
                conds_matrix(ind).pos = [48 1];
                conds_matrix(ind).gain_bias = [0 0 0 0];
             if conds_matrix(ind).dir == 1 
                conds_matrix(ind).posfunc = [1 3];
             else   
                conds_matrix(ind).posfunc = [1 4];
             end

        case 2 % Rotation: CW|CCW  
                conds_matrix(ind).arenapat = 4;
                conds_matrix(ind).mode = [0 0];
                conds_matrix(ind).pos = [48 1];
             if conds_matrix(ind).dir == 1 
                conds_matrix(ind).gain_bias = [parameters.speeds(1) 0 0 0];
             else   
                conds_matrix(ind).gain_bias = [parameters.speeds(2) 0 0 0];
             end 
                conds_matrix(ind).posfunc = [1 0];

    end %pattern selection switch board
end    

% Specify the data files
% matlab file saving
filename = 'test_'; now = datetime('now','TimeZone','local');
formatOut = 'yyyymmddHHMMSS'; fileroot = 'C:\matlabroot\data\';
fid1 = fopen([fileroot filename datestr(now,formatOut)],'w'); %opens this data file to write into it

% create folder for VIDEO data:
if isfolder([save_location_matlab '\Video']) == 0
    mkdir(save_location_matlab, '\Video');
end 
% Create folder for videos for each fly 
switch isfolder([save_location_matlab '\Video\Fly ' fly_num])
    case 0 % no folder
        mkdir([save_location_matlab '\Video\'], ['Fly ' fly_num]);
        Basler_folder_name = [save_location_matlab '\Video\Fly ' fly_num];
    case 1 % folder exists
        Basler_folder_name = [save_location_matlab '\Video\Fly ' fly_num];
end
% error for inexistant video folder
if ~isdir(Basler_folder_name)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', Basler_folder_name);
    uiwait(warndlg(errorMessage));
    return;
end

% LASER|Basler session
s_vid_light = daq.createSession('ni');
s_vid_light.Rate = 10000;
% add analog output channels for LASER|Basler
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao0', 'Voltage'); %LASER output
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao1', 'Voltage'); %Basler output

% Create sessions
% FlyCap camera trigger (tracking camera):
s = daq.createSession('ni');
s.IsContinuous = true;
parameters.cam_trig_fs = 4*parameters.fps;
parameters.acq_in_fs = 1000; %flycap sampling rate
s.Rate = parameters.acq_in_fs;

% Analog-In channels
aa = addAnalogInputChannel(s,'Dev1', [0:parameters.num_inchans-1], 'Voltage');

% Channel Outputs:
% counter for triggering FlyCap
ch = addCounterOutputChannel(s,'Dev1', 0, 'PulseGeneration'); %Flycap camera trigger
ch.Frequency = parameters.fps;

% input channels
for i = 1:length(aa)
    aa(i).TerminalConfig = 'SingleEnded';
end
aa(5).TerminalConfig = 'Differential'; % change the condition signal to differential
aa(6).TerminalConfig = 'Differential'; % change the laser input signal to differential
aa(7).TerminalConfig = 'Differential'; % change the side basler Signal to differential
aa(8).TerminalConfig = 'Differential'; % change the top basler Signal to differential

% Setup source and side_vid objects for Basler Camera
% ------ side left camera ------:
A_vid = videoinput('gentl', 2, 'Mono8');
triggerconfig(A_vid, 'hardware');                 % set trigger to come from hardware i.e. line in
A_vid.LoggingMode = 'disk&memory';
A_vid.FramesPerTrigger = parameters.Basler_fps*parameters.basler_length;
% side_vid.ROIPosition = [0 0 832 632];              % full-frame video acq.
A_vid.ROIPosition = [80 175 700 450];            % fly-only video acq.
A_Basler_src = getselectedsource(A_vid);
A_Basler_src.TriggerMode = 'Off';
% Set up for input wires from Basler Cam
A_Basler_src.LineSelector = 'Line4';             % brings up settings for line4
A_Basler_src.LineMode = 'output';                % should be 'output'; Basler cam info
A_Basler_src.LineInverter = 'False';             % should be 'False'
A_Basler_src.LineSource = 'ExposureActive';      % send out signal of when the exposure was active on basler
A_Basler_src.LineSelector = 'Line3';             % brings up settings for line3
A_Basler_src.LineMode = 'input';                 % should be 'output'; sends trig info to cam
A_Basler_src.TriggerSelector = 'FrameStart';     % start frame with trigger
A_Basler_src.TriggerMode = 'Off';
A_Basler_src.LineSelector = 'Line3';             % brings up settings for line3
A_Basler_src.TriggerActivation = 'RisingEdge';
A_Basler_src.TriggerMode = 'On';
A_Basler_src.GainAuto = 'Off';
A_Basler_src.ExposureTime = 2000;                % Exposure setting for Basler
A_Basler_src.Gain = 11.993167134727036;          % Lighting gain
A_Basler_src.Gamma = 0.6999969482421875;         % White enhancement on video
fprintf('\n Side camera configuration completed \n')

% ------back left camera ------:
B_vid = videoinput('gentl', 3, 'Mono8');
triggerconfig(B_vid, 'hardware');                 % set trigger to come from hardware i.e. line in
B_vid.LoggingMode = 'disk&memory';
B_vid.FramesPerTrigger = parameters.Basler_fps*parameters.basler_length;
% back_vid.ROIPosition = [0 0 832 632];              % full-frame video acq.
B_vid.ROIPosition = [55 160 550 360];              % fly-only video acq.
B_Basler_src = getselectedsource(B_vid);
B_Basler_src.TriggerMode = 'Off';
% Set up for input wires from Basler Cam
B_Basler_src.LineSelector = 'Line4';             % brings up settings for line4
B_Basler_src.LineMode = 'output';                % should be 'output'; Basler cam info
B_Basler_src.LineInverter = 'False';             % should be 'False'
B_Basler_src.LineSource = 'ExposureActive';      % send out signal of when the exposure was active on basler
B_Basler_src.LineSelector = 'Line3';             % brings up settings for line3
B_Basler_src.LineMode = 'input';                 % should be 'output'; sends trig info to cam
B_Basler_src.TriggerSelector = 'FrameStart';     % start frame with trigger
B_Basler_src.TriggerMode = 'Off';
B_Basler_src.LineSelector = 'Line3';             % brings up settings for line3
B_Basler_src.TriggerActivation = 'RisingEdge';
B_Basler_src.TriggerMode = 'On';
B_Basler_src.GainAuto = 'Off';
B_Basler_src.ExposureTime = 2000;                % Exposure setting for Basler
B_Basler_src.Gain = 9.3050319678579498;          % Lighting gain
B_Basler_src.Gamma = 0.7413787841796875;         % White enhancement on video
fprintf('\n Top camera configuration completed \n')

% ------ front left camera ------:
C_vid = videoinput('gentl', 1, 'Mono8');
triggerconfig(C_vid, 'hardware');                % set trigger to come from hardware i.e. line in
C_vid.LoggingMode = 'disk&memory';
C_vid.FramesPerTrigger = parameters.Basler_fps*parameters.basler_length;
C_vid.ROIPosition = [180 250 500 350];           % fly-only video acq.
% front_vid.ROIPosition = [0 0 832 632];             % full-frame video acq.
C_Basler_src = getselectedsource(C_vid);
C_Basler_src.TriggerMode = 'Off';
% Set up for input wires from Basler Cam
C_Basler_src.LineSelector = 'Line4';             % brings up settings for line4
C_Basler_src.LineMode = 'output';                % should be 'output'; Basler cam info
C_Basler_src.LineInverter = 'False';             % should be 'False'
C_Basler_src.LineSource = 'ExposureActive';      % send out signal of when the exposure was active on basler
C_Basler_src.LineSelector = 'Line3';               % brings up settings for line3
C_Basler_src.LineMode = 'input';                 % should be 'output'; sends trig info to cam
C_Basler_src.TriggerSelector = 'FrameStart';     % start frame with trigger
C_Basler_src.TriggerMode = 'Off';
C_Basler_src.LineSelector = 'Line3';             % brings up settings for line3
C_Basler_src.TriggerActivation = 'RisingEdge';
C_Basler_src.TriggerMode = 'On';
C_Basler_src.GainAuto = 'Off';
C_Basler_src.ExposureTime = 2000;                % Exposure setting for Basler
C_Basler_src.Gain = 6.7256621521589084;          % Lighting gain
C_Basler_src.Gamma = 0.6999969482421875;         % White enhancement on video
fprintf('\n Front camera configuration completed \n')

% Set up LED and Basler data queues -- ONE/STIM TYPE

signal_length = s_vid_light.Rate*parameters.basler_length;

for jj = 2:length(parameters.light_length)
    
    % ---- LED trigger data ---- %
    %ratio: ON:OFF
    LED(jj).frequency = 1200; %hz of LED frequency
    LED(jj).ratio_num = 0.5; % 1:num ratio of light off
    LED(jj).rate = round(s_vid_light.Rate/LED(jj).frequency); %adjustment ratio to get LED signal speed correct 
    LED(jj).DutyCycle = 1/(1+LED(jj).rate);
    LED(jj).light_length = parameters.light_length(jj);
    
    % creating the LED trigger signal
    LED(jj).pulse = [(ones(LED(jj).rate,1)*parameters.LED_intensity); zeros(round(LED(jj).ratio_num*LED(jj).rate),1)];
    LED(jj).pulse_length = length(LED(jj).pulse);
    LED(jj).pulse_num = round(s_vid_light.Rate/LED(jj).pulse_length)*parameters.light_length(jj); %should equal the desired light Hz  
    LED(jj).sig = LED(jj).pulse;
    %concatenate the individual pulses to reach the total length of light exposure
    for ii = 1:(LED(jj).pulse_num-1)
        LED(jj).sig = [LED(jj).sig; LED(jj).pulse];
    end  
    

    LED(jj).sig(end-(2*LED(jj).pulse_length-1):end) = 0; %two units of off at the end of the signal
    LED(jj).pre_sig = zeros(s_vid_light.Rate*parameters.basler_delay,1); %actual signal
    LED(jj).post_sig = parameters.basler_length-parameters.basler_delay-parameters.light_length(jj); %time (sec) of post LED w/basler on

    LED(jj).ON_outsig = [LED(jj).pre_sig; LED(jj).sig; zeros(round(s_vid_light.Rate*LED(jj).post_sig), 1)]; %; zeros((10*parameters.light_length(jj)),1)
    xxlength = signal_length-size(LED(jj).ON_outsig,1);
    LED(jj).ON_outsig = [LED(jj).ON_outsig; zeros(xxlength,1)];
    
    a = size(LED(jj).ON_outsig);
    LED(jj).OFF_outsig = zeros(a(1), 1);
    
end

LED(1).light_length = parameters.light_length(1);
LED(1).ON_outsig = LED(2).OFF_outsig;
LED(1).OFF_outsig = LED(2).OFF_outsig;

% ---- Basler trigger data ---- %
basler_volts = 9;
basler_outsig = zeros(a(1), 1);
basler_rate = round(s_vid_light.Rate/parameters.Basler_fps);
basler_outsig(1:basler_rate:end) = basler_volts;
clear a xxlength

%error checking analogue output signal
if size(LED(3).ON_outsig) == size(basler_outsig)
    fprintf('\n Laser signal and camera signal aligned \n')
else 
    warndlg('Mismatched laser and camera signal')
    return
end

% ---- COMBINED DATA QUEUE ---- %
for jj = 1:parameters.num_conds
    switch conds_matrix(jj).opto
        case parameters.light_length(1)
            conds_matrix(jj).LED = [LED(1).ON_outsig, basler_outsig];
        case parameters.light_length(2)
            conds_matrix(jj).LED = [LED(2).ON_outsig, basler_outsig]; 
        case parameters.light_length(3)
            conds_matrix(jj).LED = [LED(3).ON_outsig, basler_outsig]; 
        case parameters.light_length(4)
            conds_matrix(jj).LED = [LED(4).ON_outsig, basler_outsig]; 
        case parameters.light_length(5)
            conds_matrix(jj).LED = [LED(5).ON_outsig, basler_outsig]; 
        case parameters.light_length(6)
            conds_matrix(jj).LED = [LED(6).ON_outsig, basler_outsig]; 
        case parameters.light_length(7)
            conds_matrix(jj).LED = [LED(7).ON_outsig, basler_outsig];
    end
end
clear ind
%% Plot data as it is acquired and log data to file
figure(1);clf;
lh = addlistener(s,'DataAvailable', @plotData); 
lh2 = addlistener(s,'DataAvailable',@(src, event)logData(src, event, fid1));

% Send one frame to pt.grey cam to start Fictrac loading process
Panel_com('set_ao', [3 0]); 
fprintf('\n Sent one frame to FicTrac... \n');
                parameters.CLOCK(1,:) = clock;
startBackground(s); 
% pause(5); % max num frames it will send is 11.
s.stop 
  Panel_com('set_ao', [3 32767/10*9]);  %to check timing of signal out
pause(17) %wait for FicTrac to catch up
                parameters.CLOCK(3,:) = clock;
% Start the daq acquisition
Panel_com('set_ao', [3 0]);
Panel_com('all_off');
fprintf('Acquisition starting...\n');
startBackground(s); 
                parameters.CLOCK(4,:) = clock;
pause(parameters.start_pause); % Pause before starting stimulus procedure 
                parameters.CLOCK(5,:) = clock;
                
parameters.target_loop_time = ((parameters.OL_time*2)+parameters.control_time)*2.5;
% Start arena experiment
pause(0.001)
Panel_com('set_ao', [3 0]);
pause(0.001)
ind = 0;
for ii = 1:parameters.num_reps
    %set condition randomization for this rep
    conds = randperm(parameters.num_conds);
    parameters.conditions_rand(:,ii) = conds;
% ROTATE THROUGH ALL CONDITIONS:
  for jj = 1:length(conds) 
      tick_total = tic;
     tic
      ind = ind+1;
    flushdata(B_vid); flushdata(C_vid); flushdata(A_vid) 
    CC = conds(jj);                 % unrandomizes condition - actual condition 'number'   
    pat = conds_matrix(CC).pat;     % pattern ID for this condition
    Dir = conds_matrix(CC).dir;     % direction for the pattern 
    
    %queue the basler and LED:
    queueOutputData(s_vid_light, conds_matrix(CC).LED)
    
    % Prep the Basler cameras
    %A cam
    diskLogger = VideoWriter([Basler_folder_name '\' matlab_data_file  ' R' num2str(ii) 'C' num2str(CC) ' Cam-A ' conds_matrix(CC).label '.avi'],...
    'Grayscale AVI');
    diskLogger.FrameRate = parameters.Basler_fps; 
    A_vid.DiskLogger = diskLogger;  
    %B cam
    diskLogger = VideoWriter([Basler_folder_name '\' matlab_data_file ' R' num2str(ii) 'C' num2str(CC) ' Cam-B ' conds_matrix(CC).label '.avi'],...
    'Grayscale AVI');
    diskLogger.FrameRate = parameters.Basler_fps; 
    B_vid.DiskLogger = diskLogger;  
    %C cam
    diskLogger = VideoWriter([Basler_folder_name '\' matlab_data_file ' R' num2str(ii) 'C' num2str(CC) ' Cam-C ' conds_matrix(CC).label '.avi'],...
    'Grayscale AVI');
    diskLogger.FrameRate = parameters.Basler_fps; 
    C_vid.DiskLogger = diskLogger;
    
    %print rep and cond and opto info
    fprintf(['\nrep = ' num2str(ii) '|num = ' num2str(jj) ':\n']);
    fprintf(['Stim ' num2str(CC) ' | ' conds_matrix(CC).label  '\n']);
  %-------------CONTROL---------------%
  
   format shortg; CLOCK.control(ind,:) = clock;
  % Stationary Stripe:          
  pause(0.001)
    Panel_com('set_pattern_id', parameters.stripe_Pat);   
    Panel_com('set_mode', [ 4 0 ]);                    
    Panel_com('set_position', [48 1]);  
    Panel_com('set_posfunc_id', [1 0]);
  pause(0.001)
    
    pause(parameters.control_time)
    Panel_com('stop');
  %-------------STIM trial------------%
    % Cond signal:
%     Panel_com('set_ao', [3 32767/10*conds_matrix(CC).cond_sig]);
 pause(0.001)
  %prep basler cam
  
    start(A_vid);
    start(B_vid);
    start(C_vid);
 pause(0.001)
    Panel_com('set_ao', [3 32767/10*9]);  
    % Set the arena stim conditions 
    Panel_com('set_pattern_id', conds_matrix(CC).arenapat);   
    Panel_com('set_mode', conds_matrix(CC).mode);                    
    Panel_com('set_position', conds_matrix(CC).pos);  
    Panel_com('set_posfunc_id', conds_matrix(CC).posfunc);
    Panel_com('send_gain_bias', conds_matrix(CC).gain_bias);
  pause(0.001)
                format shortg; CLOCK.prestim(ind,:) = clock;
    Panel_com('start'); % start the pattern
    
    pause(parameters.OL_time-parameters.basler_delay);
   
    startBackground(s_vid_light); 
                format shortg; CLOCK.stim(ind,:) = clock;
    pause(parameters.OL_time)
                format shortg; CLOCK.post(ind,:) = clock;
    Panel_com('stop');
  pause(0.001)
    Panel_com('set_ao', [3 0]); % set cond signal back to zero
    %---- extra write to disk time as needed by the program---%%
    % Stationary Stripe:
    Panel_com('set_pattern_id', parameters.stripe_Pat);   
    Panel_com('set_mode', [ 4 0 ]);                    
    Panel_com('set_position', [48 1]);  
    Panel_com('set_posfunc_id', [1 0]);
  pause(0.001)
    % stops the Basler camera recording post stim
    stop(A_vid); 
    stop(B_vid);
    stop(C_vid);
    actual_loop_time = toc;
    parameters.writetodiskpause = parameters.target_loop_time-actual_loop_time;
    pause(parameters.writetodiskpause)
    toc(tick_total)
  end %cond loop
end %rep loop

% Stop Everything
% End with Stationary Stripe:
Panel_com('set_pattern_id', parameters.stripe_Pat);   
Panel_com('set_mode', [ 4 0 ]);                    
Panel_com('set_position', [48 1]);  
Panel_com('set_posfunc_id', [1 0]); 
Panel_com('stop');
pause(parameters.end_pause);
s.stop;
delete(lh);  delete(lh2);
fclose(fid1);
release(s)
clear s s1 

% PLOT the data
d = dir(fileroot);
fid2 = fopen([fileroot d(length(d)).name]);
[data,count] = fread(fid2,[parameters.num_inchans+1,inf],'double');
fclose(fid2);

% SAVE DATA
parameters.conds_matrix = conds_matrix;
parameters.A_Basler_src = A_Basler_src;
parameters.B_Basler_src = B_Basler_src;
parameters.C_Basler_src = C_Basler_src;
parameters.A_vid = A_vid;
parameters.B_vid = B_vid;
parameters.C_vid = C_vid;
parameters.LED = LED;
parameters.conditions_timestamps = CLOCK;
parameters.condition_label = condition_label;
parameters.directory = directory;
parameters.fly_num = fly_num;
parameters.folder_date = folder_date;
parameters.matlab_data_file = matlab_data_file;
parameters.save_location_matlab = save_location_matlab;

fly(1).data = data;
fly(1).parameters = parameters;

save('Current_data_file')
% save('TESTING_FILE', 'fly')

%% Throw data to the analysis script:
MIDSTIM_Step_2_find_conditions_and_alignment(fly)

% figure(2);hold all;
% % plot(data(1,:), 'b'); %digital counter
% plot(data(2,:), 'g'); %fly camera trigger
% plot(data(3,:), 'r'); %LED x pos
% plot(data(4,:), 'c'); %LED y pos
% plot(data(5,:), 'k'); %condition signal
% plot(data(6,:), 'm'); %FicTrac Heading
% plot(data(7,:), 'b'); %LED trigger  
% plot(data(8,:), 'color', color.orange); %Basler camera trigger
% plot(data(9,:), 'y'); %Basler camera trigger
% ylabel('volts')
% xlabel('time')
% legend(parameters.chan_labels);
% % vline(peaks, 'k')

% figure;
% hold all
% plot(data(8,:), 'r--'); %Basler camera trigger
% vline(peaks, 'k')

%% Open camera previews to realign the next fly
answer = questdlg('Open camera previews ?', 'Notes', 'Yes', 'No', 'Yes');
switch answer
    case 'Yes'
        A_Basler_src.TriggerMode = 'Off';
        B_Basler_src.TriggerMode = 'Off';   
        C_Basler_src.TriggerMode = 'Off';

        % Video sizes for POSITIONING the fly
        C_vid.ROIPosition = [320 285 500 350]; %gent 1
        B_vid.ROIPosition = [220 180 550 360]; %gent 3
        A_vid.ROIPosition = [320 250 350 300]; %gent 3

        preview(A_vid)
        preview(B_vid)
        preview(C_vid)
    case 'No'
end


