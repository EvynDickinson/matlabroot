
% RECORDS TRIAL VIDEOS FOR TESTING CAMERA FUNCTION INTO 'BASLER TRIG' FOLDER
% Next Step:
% MIDSTIM_Step_1_Arena_Script_jan2019

clear all; close all; clc
closepreview; imaqreset    

to_do = {'film', 'calibrate', 'position'};
record_type = questdlg('What type of recording do you want?', 'record type', to_do{1}, to_do{2}, to_do{3}, to_do{1});

filler = ' fly 2 '; %if you need to add anything else into the name of the vid files
Basler_folder_name = 'E:\Basler Trig\';

%% Setup source and side_vid objects for Basler Camera
% PARAMETERS:
cam_spec.Basler_fps = 300;
cam_spec.light_length = 0.5; 
cam_spec.basler_length = cam_spec.light_length+2.5;
cam_spec.basler_delay = 0.5;    %sec of basler running before LED comes on|arena stimulus  
cam_spec.start_intensity = 1;
idx = 1;  
for ii = cam_spec.start_intensity:10 
    temp.a{idx} = num2str(ii);idx = idx+1;
end

switch record_type
    case 'film'
        fly_cross = select_cross;
        intensities = listdlg('ListString', temp.a, 'SelectionMode', 'Multiple',...
                        'ListSize', [80, 150], 'PromptString', 'Laser intensity?');        
    case 'calibrate'
        fly_cross = 'Calibration Test';
        cam_spec.start_intensity = 1;
        intensities = 1;
        filler = '';
end

% Camera ROIs for positioning: 
cam_spec.CamA.ROI_positioning = [320 250 350 300]; 
cam_spec.CamB.ROI_positioning = [220 180 550 360]; 
cam_spec.CamC.ROI_positioning = [320 285 500 350]; 
% Camera ROIs for filming: 
cam_spec.CamA.ROI_film = [80 175 700 450];
cam_spec.CamB.ROI_film = [55 160 550 360];
cam_spec.CamC.ROI_film = [180 250 500 350];
% Camera ROIs for calibration: 
cam_spec.CamA.ROI_FULL = [0 0 832 632];
cam_spec.CamB.ROI_FULL = [0 0 832 632];  
cam_spec.CamC.ROI_FULL = [0 0 832 632];

% Setup source and side_vid objects for Basler Camera
FramesPerTrigger = cam_spec.Basler_fps*cam_spec.basler_length;
% Cam Specs A:
cam_spec.CamA.num = 2;
cam_spec.CamA.ROI = [80 175 700 450];
cam_spec.CamA.Gain = 11.993167134727036;
cam_spec.CamA.Gamma = 0.6999969482421875;
cam_spec.CamA.ExposureTime = 2000;
cam_spec.CamA.FramesPerTrigger = FramesPerTrigger;

% Cam Specs B:
cam_spec.CamB.num = 3;
cam_spec.CamB.ROI = [55 160 550 360];
cam_spec.CamB.Gain = 9.3050319678579498;
cam_spec.CamB.Gamma = 0.7413787841796875;
cam_spec.CamB.ExposureTime = 2000;
cam_spec.CamB.FramesPerTrigger = FramesPerTrigger;

% Cam Specs C:
cam_spec.CamC.num = 1;
cam_spec.CamC.ROI = [180 250 500 350];  
cam_spec.CamC.Gain = 6.7256621521589084;  
cam_spec.CamC.Gamma = 0.6999969482421875;
cam_spec.CamC.ExposureTime = 2000;
cam_spec.CamC.FramesPerTrigger = FramesPerTrigger;

switch record_type
    case to_do{1}
        save_vid = 1;
        cam_spec.CamA.ROI = cam_spec.CamA.ROI_film;
        cam_spec.CamB.ROI = cam_spec.CamB.ROI_film;
        cam_spec.CamC.ROI = cam_spec.CamC.ROI_film;
    case to_do{2}
        save_vid = 1;
        cam_spec.CamA.ROI = cam_spec.CamA.ROI_FULL;
        cam_spec.CamB.ROI = cam_spec.CamB.ROI_FULL;
        cam_spec.CamC.ROI = cam_spec.CamC.ROI_FULL;
        cam_spec.basler_length = 120;
        cam_spec.Basler_fps = 30;
    case to_do{3}
        save_vid = 0;
        cam_spec.CamA.ROI = cam_spec.CamA.ROI_positioning;
        cam_spec.CamB.ROI = cam_spec.CamB.ROI_positioning;
        cam_spec.CamC.ROI = cam_spec.CamC.ROI_positioning;    
end

% Grab Cameras and Initiate Set up:
[A_vid, A_Basler_src] = initiate_camera(cam_spec.CamA);
[B_vid, B_Basler_src] = initiate_camera(cam_spec.CamB);
[C_vid, C_Basler_src] = initiate_camera(cam_spec.CamC);

if save_vid == 0
    A_Basler_src.TriggerMode = 'Off'; preview(A_vid)
    B_Basler_src.TriggerMode = 'Off'; preview(B_vid)
    C_Basler_src.TriggerMode = 'Off'; preview(C_vid)
    return
end

%% SET UP DAQ SESSIONS FOR VIDEOS AND LASER
% LED|Basler session
s_vid_light = daq.createSession('ni');     
s_vid_light.Rate = 10000;

% add analog output channels for LED|Basler
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao0', 'Voltage'); %LASER output
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao1', 'Voltage'); %Basler outputs

%% CREATE ANALOGUE SIGNAL FOR LASER AND CAMERAS
for vv = 1:length(intensities)
    cam_spec.LED_intensity = intensities(vv);
    jj = 1;
    % ---- LASER trigger data ---- %
    %ratio: ON:OFF
    LED(jj).frequency = 1200; %hz of LED frequency
    LED(jj).ratio_num = 0.5; % 1:num ratio of light on:off
    LED(jj).rate = round(s_vid_light.Rate/ LED(jj).frequency); %adjustment ratio to get LED signal speed correct

    % creating the LED trigger signal
    LED(jj).pulse = [(ones(LED(jj).rate,1)*cam_spec.LED_intensity); zeros(round(LED(jj).ratio_num*LED(jj).rate),1)];
    LED(jj).pulse_length = length(LED(jj).pulse);
    LED(jj).pulse_num = round(s_vid_light.Rate/LED(jj).pulse_length)*cam_spec.light_length(jj); %should equal the desired light Hz  
    LED(jj).sig = LED(jj).pulse;
    %concatenate the individual pulses to reach the total length of light exposure
    for ii = 1:(LED(jj).pulse_num-1)
        LED(jj).sig = [LED(jj).sig; LED(jj).pulse];
    end

    LED(jj).sig(end-(2*LED(jj).pulse_length-1):end) = 0; %two units of off at the end of the signal
    LED(jj).pre_sig = zeros(s_vid_light.Rate*cam_spec.basler_delay,1);
    LED(jj).post_sig = cam_spec.basler_length-cam_spec.basler_delay-cam_spec.light_length(jj); %sec post LED w/basler on

    temp.A = zeros(s_vid_light.Rate*cam_spec.basler_delay,1);
    temp.B = LED(jj).sig;
    temp.C = zeros(s_vid_light.Rate*LED(jj).post_sig, 1);
    temp.D = zeros((10*cam_spec.light_length),1);

    LED(jj).ON_outsig = [temp.A; temp.B; temp.C; temp.D];
    a = size(LED(jj).ON_outsig);
    LED(jj).OFF_outsig = zeros(a(1), 1);

    % ---- Basler trigger data ---- %
    basler_volts = 9;
    basler_outsig = zeros(a(1), 1);
    basler_rate = round(s_vid_light.Rate/cam_spec.Basler_fps);
    basler_outsig(1:basler_rate:end) = basler_volts;
    outsig = [LED(jj).ON_outsig, basler_outsig];

    figure;
    hold all
    plot(basler_outsig)
    plot(LED(jj).ON_outsig)

    % RECORD AND WRITE VIDEO
    duty_cycle = (1/(1+LED(jj).ratio_num));

    % PREP CAMERAS
    volts = num2str(cam_spec.LED_intensity);
    %Side cam
    diskLogger = VideoWriter([Basler_folder_name date ' ' fly_cross ' ' volts 'V ' filler ' C.avi'],...
    'Grayscale AVI');
    diskLogger.FrameRate = cam_spec.Basler_fps; 
    % diskLogger.Quality = 85; 
    A_vid.DiskLogger = diskLogger;  
    %Back cam
    diskLogger = VideoWriter([Basler_folder_name date ' ' fly_cross ' ' volts 'V ' filler ' B.avi'],...
    'Grayscale AVI');
    % diskLogger = VideoWriter([Basler_folder_name date ' ' num2str(pp) ' back.avi'],...
    % 'Grayscale AVI');
    diskLogger.FrameRate = cam_spec.Basler_fps; 
    % diskLogger.Quality = 85; 
    B_vid.DiskLogger = diskLogger;  
    %Front side cam
    diskLogger = VideoWriter([Basler_folder_name date ' ' fly_cross ' ' volts 'V ' filler ' A.avi'],...
    'Grayscale AVI');
    diskLogger.FrameRate = cam_spec.Basler_fps; 
    % diskLogger.Quality = 85; 
    C_vid.DiskLogger = diskLogger;  
    %  'VideoFormat', 'Grayscale', 'Motion JPEG AVI'

    
    
    % START RECORDING
    fprintf('\n STARTING VIDEOS NOW \n')
    queueOutputData(s_vid_light,  outsig) 
    
    
    
    tic
%     if save_vid == 1
        start(A_vid); start(B_vid); start(C_vid);
%     end
    startBackground(s_vid_light); 
    
    while (C_vid.FramesAcquired ~= C_vid.DiskLoggerFrameCount)
        pause(.1)
    end
    
%     pause(cam_spec.basler_length * 1.2)
%     fprintf('\n Videos finished, writing to disk now... \n')
%     if save_vid == 1
        stop(A_vid); stop(B_vid); stop(C_vid);
%     end
    s_vid_light.stop  
    
    toc
    
    
    fprintf([' \n Done w/ intensity: '  num2str(cam_spec.LED_intensity) 'V\n'])  
end    
 fprintf('DONE \n')
 
 close all
%% OPEN PREVIEW FOR EACH CAMERA
% 
% % Video sizes for POSITIONING the fly
% C_vid.ROIPosition = [320 285 500 350]; %gent 1
% B_vid.ROIPosition = [220 180 550 360]; %gent 3
% A_vid.ROIPosition = [320 250 350 300]; %gent 3
% 
% close all     
% 
% [cameras.A]
% offset = [224, 120, 832, 704]
% 
% [cameras.B]
% offset = [224, 230, 800, 700]
% 
% [cameras.C]
% offset = [160, 250, 864, 700]
%  s = struct()
% s.cameras.A.offset = [224, 120, 832, 704];
% s.cameras.B.offset = [224, 120, 832, 704];
% s.cameras.C.offset = [224, 120, 832, 704];
% 
% addpath C:/Users/Tuthill/Pierre/builds/matlab-toml
% toml.write('config.toml', s)



% 
% hAvid = preview(A_vid);
% set(hAvid.Parent, 'position',  [0.5    0.5  704.0000  450.0000])
% 
% objProperties.Parent.Position
% 
% 
% [0.85,0.47,0.36,0.43]
% 
%  preview(A_vid)
%     B_Basler_src.TriggerMode = 'Off'; preview(B_vid)
%     C_Basler_src.TriggerMode = 'Off'; preview(C_vid)


