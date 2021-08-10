
% RECORDS TRIAL VIDEOS FOR TESTING CAMERA FUNCTION INTO 'BASLER TRIG' FOLDER
% Next Step:
% MIDSTIM_Step_1_Arena_Script_jan2019

clear all; close all; clc  
closepreview; imaqreset  

to_do = {'film', 'calibrate', 'position'};
record_type = questdlg('What type of recording do you want?', 'record type', to_do{1}, to_do{2}, to_do{3}, to_do{1});
filler = ' lens '; %if you need to add anything else into the name of the vid files
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
                        'ListSize', [80, 300], 'PromptString', 'Laser intensity?');  
        if isempty(intensities)
            intensities = str2double(cell2mat((inputdlg('Laser Intensity?'))));
        end
    case 'calibrate'
        fly_cross = 'Calibration Test';
        cam_spec.start_intensity = 1;
        intensities = 1;
        filler = ' TEST 3';
end
num.cams = 6;

CamList = {'CamF', 'CamB', 'CamC', 'CamA', 'CamE', 'CamD'}; %order of cameras

% Setup source and side_vid objects for Basler Camera
FramesPerTrigger = cam_spec.Basler_fps*cam_spec.basler_length;
cam_spec = load_cam_parameters(cam_spec, num, FramesPerTrigger);

% Live Feed Adjustments:  
% E_vid.ROIPosition = [ 200 260 200 136];
% F_Basler_src.Gain = 6.7256621521589084;
% F_Basler_src.Gamma = 0.899993896484375;

switch record_type
    case to_do{1}
        save_vid = 1; % FILM
        for icam = 1:num.cams
            cam_spec.(['Cam' Alphabet(icam)]).ROI = cam_spec.(['Cam' Alphabet(icam)]).ROI_film;
        end
    case to_do{2} % CALIBRATE
        save_vid = 1;
        cam_spec.basler_length = 60;
        cam_spec.Basler_fps = 30;
        for icam = 1:num.cams
            cam_spec.(['Cam' Alphabet(icam)]).FramesPerTrigger = cam_spec.Basler_fps*cam_spec.basler_length;
            cam_spec.(['Cam' Alphabet(icam)]).ROI = cam_spec.(['Cam' Alphabet(icam)]).ROI_FULL;
        end
    case to_do{3} % POSITION
        save_vid = 0;
        for icam = 1:num.cams
            cam_spec.(['Cam' Alphabet(icam)]).ROI = cam_spec.(['Cam' Alphabet(icam)]).ROI_positioning; %ROI_positioning;
        end    
end


% Grab Cameras and Initiate Set up:
[A_vid, A_Basler_src] = initiate_camera(cam_spec.CamA);
[B_vid, B_Basler_src] = initiate_camera(cam_spec.CamB);
[C_vid, C_Basler_src] = initiate_camera(cam_spec.CamC);
[D_vid, D_Basler_src] = initiate_camera(cam_spec.CamD);
[E_vid, E_Basler_src] = initiate_camera(cam_spec.CamE);
[F_vid, F_Basler_src] = initiate_camera(cam_spec.CamF);

if save_vid == 0
    A_Basler_src.TriggerMode = 'Off'; preview(A_vid)
    B_Basler_src.TriggerMode = 'Off'; preview(B_vid)
    C_Basler_src.TriggerMode = 'Off'; preview(C_vid)
    D_Basler_src.TriggerMode = 'Off'; preview(D_vid)
    E_Basler_src.TriggerMode = 'Off'; preview(E_vid)
    F_Basler_src.TriggerMode = 'Off'; preview(F_vid)
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
    
    tic
    cam_spec.LED_intensity = 4;% intensities(vv);
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

%     figure;
%     hold all
%     plot(basler_outsig)
%     plot(LED(jj).ON_outsig)

    % RECORD AND WRITE VIDEO
    duty_cycle = (1/(1+LED(jj).ratio_num));

    % PREP CAMERAS
    volts = num2str(cam_spec.LED_intensity);
    
    % SAVING STRUCTURES FOR THE CAMERAS
    
    for icam = 1:num.cams
        diskLogger = VideoWriter([Basler_folder_name date ' ' fly_cross ' ' volts 'V ' filler ' Cam-' Alphabet(icam) '.avi'],...
        'Grayscale AVI');
        diskLogger.FrameRate = cam_spec.Basler_fps; 
        switch icam
            case 1; A_vid.DiskLogger = diskLogger;    
            case 2; B_vid.DiskLogger = diskLogger; 
            case 3; C_vid.DiskLogger = diskLogger; 
            case 4; D_vid.DiskLogger = diskLogger;  
            case 5; E_vid.DiskLogger = diskLogger; 
            case 6; F_vid.DiskLogger = diskLogger; 
        end
    end
    
    
    % START RECORDING
    fprintf('\n STARTING VIDEOS NOW \n')
    queueOutputData(s_vid_light,  outsig) 
    if save_vid == 1
        start(A_vid); 
        start(B_vid); start(C_vid);
        start(D_vid); start(E_vid); start(F_vid);
    end
    startBackground(s_vid_light) 
    
    if strcmpi(record_type, to_do{2})
        preview(B_vid);preview(D_vid);preview(F_vid)
    end
    
    pause(cam_spec.basler_length * 1.0)
    fprintf('\n Videos finished, writing to disk now... \n')
    if save_vid == 1
        stop(A_vid); 
        stop(B_vid); stop(C_vid);
        stop(D_vid); stop(E_vid); stop(F_vid);
    end
    s_vid_light.stop  
    fprintf([' \n Done w/ intensity: '  num2str(cam_spec.LED_intensity) 'V\n'])  
    toc
end    
 fprintf('DONE \n')
 
 close all
 
% size(outsig)
% 1200006/s_vid_light.Rate



