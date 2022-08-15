


% RECORDS TRIAL VIDEOS FOR TESTING CAMERA FUNCTION INTO 'BASLER TRIG' FOLDER
clear all; close all; clc
closepreview     
save_vid = 1;


%% Setup source and side_vid objects for Basler Camera
% PARAMETERS:
parameters.Basler_fps = 10;
parameters.light_length = 0.5;
parameters.basler_delay = .5; %sec of basler running before LED comes on|arena stimulus  
parameters.LED_intensity = 4.5;
parameters.basler_length = 2;

% parameters.basler_length = parameters.light_length+2;

clc
imaqreset

% ------ side left camera ------:  
side_vid = videoinput('gentl', 2, 'Mono8');
triggerconfig(side_vid, 'hardware');                 % set trigger to come from hardware i.e. line in
side_vid.LoggingMode = 'disk&memory';
side_vid.FramesPerTrigger = parameters.Basler_fps*parameters.basler_length;
% side_vid.ROIPosition = [0 0 832 632];              % full-frame video acq.
side_vid.ROIPosition = [80 175 700 450];             % fly-only video acq.
side_Basler_src = getselectedsource(side_vid);
side_Basler_src.TriggerMode = 'Off';
% Set up for input wires from Basler Cam  
side_Basler_src.LineSelector = 'Line4';             % brings up settings for line4
side_Basler_src.LineMode = 'output';                % should be 'output'; Basler cam info
side_Basler_src.LineInverter = 'False';             % should be 'False'
side_Basler_src.LineSource = 'ExposureActive';      % send out signal of when the exposure was active on basler
side_Basler_src.LineSelector = 'Line3';               % brings up settings for line3
side_Basler_src.LineMode = 'input';                 % should be 'output'; sends trig info to cam
side_Basler_src.TriggerSelector = 'FrameStart';     % start frame with trigger
side_Basler_src.TriggerMode = 'Off';
side_Basler_src.LineSelector = 'Line3';             % brings up settings for line3
side_Basler_src.TriggerActivation = 'RisingEdge';
side_Basler_src.TriggerMode = 'On';
side_Basler_src.GainAuto = 'Off';
side_Basler_src.ExposureTime = 2000;                % Exposure setting for Basler
side_Basler_src.Gain = 11.993167134727036;          % Lighting gain
side_Basler_src.Gamma = 0.6999969482421875;         % White enhancement on video
fprintf('\n Side camera configuration completed \n')

% ------ back left camera ------:
back_vid = videoinput('gentl', 3, 'Mono8');
triggerconfig(back_vid, 'hardware');                 % set trigger to come from hardware i.e. line in
back_vid.LoggingMode = 'disk&memory';
back_vid.FramesPerTrigger = parameters.Basler_fps*parameters.basler_length;
% back_vid.ROIPosition = [0 0 832 632];              % full-frame video acq.
back_vid.ROIPosition = [55 160 550 360];             % fly-only video acq.
back_Basler_src = getselectedsource(back_vid);
back_Basler_src.TriggerMode = 'Off';
% Set up for input wires from Basler Cam
back_Basler_src.LineSelector = 'Line4';             % brings up settings for line4
back_Basler_src.LineMode = 'output';                % should be 'output'; Basler cam info
back_Basler_src.LineInverter = 'False';             % should be 'False'
back_Basler_src.LineSource = 'ExposureActive';      % send out signal of when the exposure was active on basler
back_Basler_src.LineSelector = 'Line3';             % brings up settings for line3
back_Basler_src.LineMode = 'input';                 % should be 'output'; sends trig info to cam
back_Basler_src.TriggerSelector = 'FrameStart';     % start frame with trigger
back_Basler_src.TriggerMode = 'Off';
back_Basler_src.LineSelector = 'Line3';             % brings up settings for line3
back_Basler_src.TriggerActivation = 'RisingEdge';
back_Basler_src.TriggerMode = 'On';
back_Basler_src.GainAuto = 'Off';
back_Basler_src.ExposureTime = 2000;                % Exposure setting for Basler
back_Basler_src.Gain = 9.3050319678579498;          % Lighting gain
back_Basler_src.Gamma = 0.7413787841796875;         % White enhancement on video
fprintf('\n Top camera configuration completed \n')

% ------ front left camera ------:
front_vid = videoinput('gentl', 1, 'Mono8');
triggerconfig(front_vid, 'hardware');                 % set trigger to come from hardware i.e. line in
front_vid.LoggingMode = 'disk&memory';
front_vid.FramesPerTrigger = parameters.Basler_fps*parameters.basler_length;
front_vid.ROIPosition = [180 250 500 350];              %front leg only  video acq.
% front_vid.ROIPosition = [0 0 832 632];             % full-frame video acq.
front_Basler_src = getselectedsource(front_vid);
front_Basler_src.TriggerMode = 'Off';
% Set up for input wires from Basler Cam
front_Basler_src.LineSelector = 'Line4';             % brings up settings for line4
front_Basler_src.LineMode = 'output';                % should be 'output'; Basler cam info
front_Basler_src.LineInverter = 'False';             % should be 'False'
front_Basler_src.LineSource = 'ExposureActive';      % send out signal of when the exposure was active on basler
front_Basler_src.LineSelector = 'Line3';               % brings up settings for line3
front_Basler_src.LineMode = 'input';                 % should be 'output'; sends trig info to cam
front_Basler_src.TriggerSelector = 'FrameStart';     % start frame with trigger
front_Basler_src.TriggerMode = 'Off';
front_Basler_src.LineSelector = 'Line3';             % brings up settings for line3
front_Basler_src.TriggerActivation = 'RisingEdge';
front_Basler_src.TriggerMode = 'On';
front_Basler_src.GainAuto = 'Off';
front_Basler_src.ExposureTime = 2000;                      % Exposure setting for Basler
front_Basler_src.Gain = 6.7256621521589084;          % Lighting gain
front_Basler_src.Gamma = 0.6999969482421875;         % White enhancement on video
fprintf('\n Front camera configuration completed \n')

%% SET UP DAQ SESSIONS FOR VIDEOS AND LASER
% LED|Basler session
s_vid_light = daq.createSession('ni');     
s_vid_light.Rate = 10000;

% add analog output channels for LED|Basler
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao0', 'Voltage'); %LASER output
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao1', 'Voltage'); %Basler outputs

%% CREATE ANALOGUE SIGNAL FOR LASER AND CAMERAS

for vv = 7:max_intensity
parameters.LED_intensity = vv;

jj = 1;
% ---- LASER trigger data ---- %
%ratio: ON:OFF
LED(jj).frequency = 1200; %hz of LED frequency
LED(jj).ratio_num = 0.5; % 1:num ratio of light on:off
LED(jj).rate = round(s_vid_light.Rate/ LED(jj).frequency); %adjustment ratio to get LED signal speed correct

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
LED(jj).pre_sig = zeros(s_vid_light.Rate*parameters.basler_delay,1);
LED(jj).post_sig = parameters.basler_length-parameters.basler_delay-parameters.light_length(jj); %sec post LED w/basler on

LED(jj).ON_outsig = [zeros(s_vid_light.Rate*parameters.basler_delay,1); LED(jj).sig; zeros(s_vid_light.Rate*LED(jj).post_sig, 1); zeros((10*parameters.light_length),1)];
a = size(LED(jj).ON_outsig);
LED(jj).OFF_outsig = zeros(a(1), 1);

% ---- Basler trigger data ---- %
basler_volts = 9;
basler_outsig = zeros(a(1), 1);
basler_rate = round(s_vid_light.Rate/parameters.Basler_fps);
basler_outsig(1:basler_rate:end) = basler_volts;
outsig = [LED(jj).ON_outsig, basler_outsig];

figure;
hold all
plot(basler_outsig)
plot(LED(jj).ON_outsig)

% ERROR CHECKING/TROUBLE SHOOTING
% size(LED(jj).ON_outsig)
% size(basler_outsig)
% clear a
% figure; hold all
% plot(basler_outsig)
% plot(LED(jj).ON_outsig)

%% RECORD AND WRITE VIDEO
for pp = 1:1
    Basler_folder_name = 'E:\Basler Trig\';
% Basler_folder_name = 'C:\matlabroot\basler_trig\';

duty_cycle = (1/(1+LED(jj).ratio_num));

% PREP CAMERAS
%Side cam
diskLogger = VideoWriter([Basler_folder_name date ' ' fly_cross ' ' num2str(vv) 'V ' filler ' C.avi'],...
'Grayscale AVI');
diskLogger.FrameRate = parameters.Basler_fps; 
% diskLogger.Quality = 85; 
side_vid.DiskLogger = diskLogger;  
%Back cam
diskLogger = VideoWriter([Basler_folder_name date ' ' fly_cross ' ' num2str(vv) 'V ' filler ' B.avi'],...
'Grayscale AVI');
% diskLogger = VideoWriter([Basler_folder_name date ' ' num2str(pp) ' back.avi'],...
% 'Grayscale AVI');
diskLogger.FrameRate = parameters.Basler_fps; 
% diskLogger.Quality = 85; 
back_vid.DiskLogger = diskLogger;  
%Front side cam
diskLogger = VideoWriter([Basler_folder_name date ' ' fly_cross ' ' num2str(vv) 'V ' filler ' A.avi'],...
'Grayscale AVI');
diskLogger.FrameRate = parameters.Basler_fps; 
% diskLogger.Quality = 85; 
front_vid.DiskLogger = diskLogger;  
%  'VideoFormat', 'Grayscale', 'Motion JPEG AVI'

% START RECORDING
fprintf('\n STARTING VIDEOS NOW \n')
queueOutputData(s_vid_light,  outsig) 
if save_vid == 1
    start(side_vid);
    start(back_vid);
    start(front_vid);
end

startBackground(s_vid_light); 
pause(parameters.basler_length *1.5)
fprintf('\n Videos finished, writing to disk now... \n')

if save_vid == 1
    stop(side_vid); 
    stop(back_vid);
    stop(front_vid);
end
s_vid_light.stop
fprintf(' \n done \n')  
end    
end    
%% OPEN PREVIEW FOR EACH CAMERA
side_Basler_src.TriggerMode = 'Off';
back_Basler_src.TriggerMode = 'Off';   
front_Basler_src.TriggerMode = 'Off';

% Video sizes for POSITIONING the fly
front_vid.ROIPosition = [320 285 500 350]; %gent 1
back_vid.ROIPosition = [220 180 550 360]; %gent 3
side_vid.ROIPosition = [320 250 350 300]; %gent 3

preview(side_vid)
preview(back_vid)
preview(front_vid)

% close all     











