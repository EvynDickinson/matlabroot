function varargout = start_gui(varargin)
% START_GUI MATLAB code for start_gui.fig
%      START_GUI, by itself, creates a new START_GUI or raises the existing
%      singleton*.
%
%      H = START_GUI returns the handle to a new START_GUI or the handle to
%      the existing singleton*.
%
%      START_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in START_GUI.M with the given input arguments.
%
%      START_GUI('Property','Value',...) creates a new START_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before start_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to start_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help start_gui

% Last Modified by GUIDE v2.5 19-May-2020 13:04:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @start_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @start_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% --- Executes just before start_gui is made visible.
function start_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to start_gui (see VARARGIN)

% handles.fighandles = load();
handles.figure1.Position = [256 3.6923 124.4000 42.3077];

% fill in the data to the fly cross list:
file_id = fopen('Cross.txt');
A = fscanf(file_id, '%s');
Cross = strsplit(A,','); 
set(handles.listbox1,'String', Cross);

% Choose default command line output for start_gui
handles.output = hObject;

% set all cameras to ON
SelectAllCamsButton_Callback(hObject, eventdata, handles)


% handles.posfig1 = handles.figure1.Position;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes start_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = start_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

% contents = cellstr(get(hObject,'String')) %returns listbox1 contents as cell array
% contents{get(hObject,'Value')} % returns selected item from listbox1


% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in FILM BUTTON.
function filmbutton_Callback(hObject, eventdata, handles)

pushbutton7_Callback(hObject, eventdata, handles) % close all figs and reset cams

% Select the saving pathway:
Basler_folder_name = 'E:\Basler Trig\';
Cross = get(handles.listbox1,'String');
idx = get(handles.listbox1,'Value');
fly_cross = Cross{idx};
handles.fly_cross = fly_cross;

% get the laser intensity   %this could be seriously cleaned up....
laser_idx = get(handles.listbox5, 'Value');
all_intensities = {'1' '1.2' '1.4' '1.6' '1.8' '2.0' '2.5' '3.0' '3.5' '4.0' '4.5' '5.0' '6.0' '7.0' '8.0' '9.0'};
ALL_intensities = [1 1.2 1.4 1.6 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0 6.0 7.0 8.0 9.0];
intensity = all_intensities{laser_idx(1)};

% get light length:
lengths = [0.1 0.2 0.3 0.5 0.72 1 1.5 2 2.5 3];
light_idx = get(handles.listbox4, 'Value');
light_length = lengths(light_idx);
set(handles.text6, 'String', light_length)

% PARAMETERS: 
cam_spec.Basler_fps = 300;
cam_spec.light_length = light_length; 
cam_spec.basler_length = cam_spec.light_length+2.5;
cam_spec.basler_delay = 0.5;    %sec of basler running before light comes on  

% sort and find the cameras that are being selected to run:
cam = getCamList(hObject, eventdata, handles); %find toggled cameras
num.cam = sum(cam); %total number of selected cameras
%       end tag = [400, 395, 083, 483, 635, 652];
%fullOrderNames = {'A', 'B', 'C', 'D', 'E', 'F'};
fullOrderNums = [4, 2, 3, 6, 5, 1]; % order that cameras are read in by their serial number (small2large)
loc = logical(cam); % true or false binary 
fullOrderNums(~loc) = NaN; %give unselected cameras a NaN value

idx = sort(fullOrderNums, 'ascend'); %sort the existing camera numbers by size
num.camNums = NaN(1,6); %blank set the size of the num of cameras (total, not selected)
for ii = 1:num.cam %for each selected camera:
   a = fullOrderNums==(idx(ii)); %find the orginal number in the list of all cameras
   num.camNums(a) = ii; %assign that new input number to the camera 
end

pushbutton7_Callback(hObject, eventdata, handles)

% Setup source and side_vid objects for Basler Camera
FramesPerTrigger = cam_spec.Basler_fps*cam_spec.basler_length;
cam_spec = load_cam_parameters(cam_spec, num, FramesPerTrigger);
for icam = 1:6
    cam_spec.(['Cam' Alphabet(icam)]).ROI = cam_spec.(['Cam' Alphabet(icam)]).ROI_film;
end

% Grab Cameras and Initiate Set up:
if cam(1) == true
    A_vid = initiate_camera(cam_spec.CamA);
end
if cam(2) == true
    B_vid = initiate_camera(cam_spec.CamB);
end
if cam(3) == true
    C_vid = initiate_camera(cam_spec.CamC);
end
if cam(4) == true
    D_vid = initiate_camera(cam_spec.CamD);
end
if cam(5) == true
    E_vid = initiate_camera(cam_spec.CamE);
end
if cam(6) == true
    F_vid = initiate_camera(cam_spec.CamF);
end

% LED|Basler session
s_vid_light = daq.createSession('ni');     
s_vid_light.Rate = 10000;

% add analog output channels for LED|Basler
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao0', 'Voltage'); %LASER output
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao1', 'Voltage'); %Basler outputs

% get the file name filler:
file_tag = get(handles.edit1, 'String');

for vv = 1:length(laser_idx)
    tic
    cam_spec.LED_intensity = ALL_intensities(laser_idx(vv));
    
    % set the pathway into the file name spot
    Basler_folder_root = [Basler_folder_name date ' ' fly_cross ' ' num2str(cam_spec.LED_intensity) 'V ' file_tag];
    set(handles.text6, 'String', Basler_folder_root);
    handles.Basler_folder_root = Basler_folder_root;
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

    % PREP CAMERAS
    volts = num2str(cam_spec.LED_intensity);
    % SAVING STRUCTURES FOR THE CAMERAS
    
    for icam = 1:6 %do for all cameras, but later only denoted the selected ones
        if cam(icam) == 1
            diskLogger = VideoWriter([Basler_folder_root ' Cam-' Alphabet(icam) '.avi'],'Grayscale AVI');
            diskLogger.FrameRate = cam_spec.Basler_fps; 
            switch icam
                case 1; A_vid.DiskLogger = diskLogger; % start(A_vid);   
                case 2; B_vid.DiskLogger = diskLogger; % start(B_vid);
                case 3; C_vid.DiskLogger = diskLogger; % start(C_vid);
                case 4; D_vid.DiskLogger = diskLogger; % start(D_vid);
                case 5; E_vid.DiskLogger = diskLogger; % start(E_vid);
                case 6; F_vid.DiskLogger = diskLogger; % start(F_vid);
            end
        end
    end
    
    for icam = 1:6 %do for all cameras, but later only denoted the selected ones
        if cam(icam) == 1
            switch icam
                case 1; start(A_vid);   
                case 2; start(B_vid);
                case 3; start(C_vid);
                case 4; start(D_vid);
                case 5; start(E_vid);
                case 6; start(F_vid);
            end
        end
    end
    
    % START RECORDING
    fprintf('\n STARTING VIDEOS NOW \n')
    queueOutputData(s_vid_light,  outsig) 
    % start the camera recordings:
    startBackground(s_vid_light) 
    pause(cam_spec.basler_length * 1.0)
    fprintf('\n Videos finished, writing to disk now... \n')
    % save the recordings:
    for icam = 1:6 %do for all cameras, but later only denoted the selected ones
        if cam(icam) == 1
            switch icam
                case 1; stop(A_vid);   
                case 2; stop(B_vid);
                case 3; stop(C_vid);
                case 4; stop(D_vid);
                case 5; stop(E_vid);
                case 6; stop(F_vid);
            end
        end
    end
%     stop(A_vid); stop(B_vid); stop(C_vid); stop(D_vid); stop(E_vid); stop(F_vid);
    s_vid_light.stop  
    fprintf([' \n Done w/ intensity: '  num2str(cam_spec.LED_intensity) 'V\n'])  
    toc   
end
fprintf('DONE \n')

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton2.% CALIBRATE OPTION
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pushbutton7_Callback(hObject, eventdata, handles)

Basler_folder_name = 'E:\Basler Trig\';
fly_cross = 'Calibration Test';
filler = get(handles.edit1, 'String');
Basler_folder_root = [Basler_folder_name date ' ' fly_cross ' ' filler];
set(handles.text6, 'String', Basler_folder_root)
cam_spec.Basler_fps = 30;
cam_spec.light_length = 0.5; 
cam_spec.basler_length = 60;    %record for 60 seconds
cam_spec.basler_delay = 0.5;    %sec of basler running before LED comes on|arena stimulus  
cam_spec.start_intensity = 1;
% Setup source and side_vid objects for Basler Camera
num.cams = 6;
CamList = {'CamF', 'CamB', 'CamC', 'CamA', 'CamE', 'CamD'}; %order of cameras
FramesPerTrigger = cam_spec.Basler_fps*cam_spec.basler_length;
cam_spec = load_cam_parameters(cam_spec, num, FramesPerTrigger);
for icam = 1:num.cams
    cam_spec.(['Cam' Alphabet(icam)]).FramesPerTrigger = cam_spec.Basler_fps*cam_spec.basler_length;
    cam_spec.(['Cam' Alphabet(icam)]).ROI = cam_spec.(['Cam' Alphabet(icam)]).ROI_FULL;
end
% Grab Cameras and Initiate Set up:
[A_vid, ~] = initiate_camera(cam_spec.CamA);
[B_vid, ~] = initiate_camera(cam_spec.CamB);
[C_vid, ~] = initiate_camera(cam_spec.CamC);
[D_vid, ~] = initiate_camera(cam_spec.CamD);
[E_vid, ~] = initiate_camera(cam_spec.CamE);
[F_vid, ~] = initiate_camera(cam_spec.CamF);
% LED|Basler session
s_vid_light = daq.createSession('ni');     
s_vid_light.Rate = 10000;

% add analog output channels for LED|Basler
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao0', 'Voltage'); %LASER output
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao1', 'Voltage'); %Basler outputs

tic
    cam_spec.LED_intensity = 0.5;
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

    % PREP CAMERAS
    for icam = 1:num.cams
        diskLogger = VideoWriter([Basler_folder_root ' Cam-' Alphabet(icam) '.avi'],'Grayscale AVI');
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
    fprintf('\n STARTING VIDEOS NOW \n'); 
    load splat; sound(y, Fs*(2/3))
    queueOutputData(s_vid_light,  outsig) 
    start(A_vid); start(B_vid); start(C_vid); start(D_vid); start(E_vid); start(F_vid);
    startBackground(s_vid_light) 
    % preview cameras while recording
%     preview(D_vid)
    preview(A_vid)
%         preview(B_vid);;preview(F_vid); ;preview(C_vid);preview(E_vid)
    pause(cam_spec.basler_length * 1.0)
    fprintf('\n Videos finished, writing to disk now... \n')
    stop(A_vid); stop(B_vid); stop(C_vid); stop(D_vid); stop(E_vid); stop(F_vid);
    sound(y, Fs*(2/3))
    s_vid_light.stop  
    toc
 fprintf('DONE \n')
closepreview


% --- Executes on button press in pushbutton3 / POSITIONING PREVIEWS.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Setup source and side_vid objects for Basler Camera

%find which cameras are being used:
cam = getCamList(hObject, eventdata, handles);
%        end tag = [400, 395, 083, 483, 635, 652];
% fullOrderNames = {'A', 'B', 'C', 'D', 'E', 'F'};
fullOrderNums = [4, 2, 3, 6, 5, 1];
loc = logical(cam);
fullOrderNums(~loc) = NaN;
[~,num.camNums] = sort(fullOrderNums);

 
pushbutton7_Callback(hObject, eventdata, handles)

cam_spec.Basler_fps = 300;
cam_spec.light_length = 0.5; 
cam_spec.basler_length = cam_spec.light_length+2.5;
cam_spec.basler_delay = 0.5;    %sec of basler running before light comes on  
num.cams = sum(cam);


FramesPerTrigger = cam_spec.Basler_fps*cam_spec.basler_length;
cam_spec = load_cam_parameters(cam_spec, num, FramesPerTrigger);
for icam = 1:6
    cam_spec.(['Cam' Alphabet(icam)]).ROI = cam_spec.(['Cam' Alphabet(icam)]).ROI_positioning; %ROI_positioning;
end 

% adjust camera gentl-# for the number of cameras running:



% Grab Cameras and Initiate Set up:
if cam(1) == true
    [A_vid, A_Basler_src] = initiate_camera(cam_spec.CamA);
    A_Basler_src.TriggerMode = 'Off'; preview(A_vid)
end
if cam(2) == true
    [B_vid, B_Basler_src] = initiate_camera(cam_spec.CamB);
    B_Basler_src.TriggerMode = 'Off'; preview(B_vid)
end
if cam(3) == true
    [C_vid, C_Basler_src] = initiate_camera(cam_spec.CamC);
    C_Basler_src.TriggerMode = 'Off'; preview(C_vid)
end
if cam(4) == true
    [D_vid, D_Basler_src] = initiate_camera(cam_spec.CamD);
    D_Basler_src.TriggerMode = 'Off'; preview(D_vid)
end
if cam(5) == true
    [E_vid, E_Basler_src] = initiate_camera(cam_spec.CamE);
    E_Basler_src.TriggerMode = 'Off'; preview(E_vid)
end
if cam(6) == true
    [F_vid, F_Basler_src] = initiate_camera(cam_spec.CamF);
    F_Basler_src.TriggerMode = 'Off'; preview(F_vid)
end


% --- Executes on button press in pushbutton6 CAMERA PREVIEW.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pushbutton7_Callback(hObject, eventdata, handles)
 warning('off', 'imaq:gentl:hardwareTriggerTriggerModeOff');
 warning('off', 'imaq:gentl:adaptorSetROIModified');
% find which cameras are being used:
cam = getCamList(hObject, eventdata, handles);
num.cam = sum(cam);
%       end tag = [400, 395, 083, 483, 635, 652];
%fullOrderNames = {'A', 'B', 'C', 'D', 'E', 'F'};
fullOrderNums = [4, 2, 3, 6, 5, 1];
loc = logical(cam);
fullOrderNums(~loc) = NaN; %only use selected cameras

idx = sort(fullOrderNums, 'ascend');
num.camNums = NaN(1,6);
for ii = 1:sum(cam)
   a = fullOrderNums==(idx(ii));
   num.camNums(a) = ii;
end

pushbutton7_Callback(hObject, eventdata, handles)

cam_spec.Basler_fps = 300;
cam_spec.light_length = 0.5; 
cam_spec.basler_length = cam_spec.light_length+2.5;
cam_spec.basler_delay = 0.5;    %sec of basler running before light comes on  
% num.cams = 6;
FramesPerTrigger = cam_spec.Basler_fps*cam_spec.basler_length;
cam_spec = load_cam_parameters(cam_spec, num, FramesPerTrigger);
for icam = 1:6
    cam_spec.(['Cam' Alphabet(icam)]).ROI = cam_spec.(['Cam' Alphabet(icam)]).ROI_film; %ROI_positioning;
end 

z = 1.0; %zoom ratio for previews

% Grab Cameras and Initiate Set up:
if cam(1) == true
    [A_vid, A_Basler_src] = initiate_camera(cam_spec.CamA);
    A_Basler_src.TriggerMode = 'Off'; 
    AFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam A 400 gentl-4');
    n = cam_spec.CamA.ROI_film(4);
    m = cam_spec.CamA.ROI_film(3);
    hImage = image(zeros(n, m, 1));
    AFig.Position = [800 570 m*z n*z];
    handles.Cams.AFig = AFig; %save the figure handle
    preview(A_vid, hImage);
end
if cam(2) == true
    [B_vid, B_Basler_src] = initiate_camera(cam_spec.CamB);
    B_Basler_src.TriggerMode = 'Off'; 
    BFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam B 395 gentl-2');
    n = cam_spec.CamB.ROI_film(4);
    m = cam_spec.CamB.ROI_film(3);
    hImage = image(zeros(n, m, 1));
    BFig.Position = [1338 633 m n];
    preview(B_vid, hImage);
    handles.Cams.BFig = BFig; %save the figure handle
end
if cam(3) == true
    [C_vid, C_Basler_src] = initiate_camera(cam_spec.CamC);
    C_Basler_src.TriggerMode = 'Off'; 
    CFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam C 083 gentl-3');
    n = cam_spec.CamC.ROI_film(4);
    m = cam_spec.CamC.ROI_film(3);
    hImage = image(zeros(n, m, 1));
    CFig.Position = [30 500 m n];
    preview(C_vid, hImage);
    handles.Cams.CFig = CFig; %save the figure handle
end
if cam(4) == true
    [D_vid, D_Basler_src] = initiate_camera(cam_spec.CamD);
    D_Basler_src.TriggerMode = 'Off';
    DFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam D 483 gentl-6');
    n = cam_spec.CamD.ROI_film(4);
    m = cam_spec.CamD.ROI_film(3);
    hImage = image(zeros(n, m, 1));
    DFig.Position = [910 98 m n];
    preview(D_vid, hImage);
    handles.Cams.DFig = DFig; %save the figure handle
end
if cam(5) == true
    [E_vid, E_Basler_src] = initiate_camera(cam_spec.CamE);
    E_Basler_src.TriggerMode = 'Off'; 
    EFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam E 635 gentl-5');
    n = cam_spec.CamE.ROI_film(4);
    m = cam_spec.CamE.ROI_film(3);
    hImage = image(zeros(n, m, 1));
    EFig.Position = [395 90 m n];
    preview(E_vid, hImage);
    handles.Cams.EFig = EFig; %save the figure handle
end
if cam(6) == true
    [F_vid, F_Basler_src] = initiate_camera(cam_spec.CamF);
    F_Basler_src.TriggerMode = 'Off'; 
    FFig = figure('Toolbar','none','Menubar', 'none',...
       'NumberTitle','Off','Name','Cam F 652 gentl-1');
    n = cam_spec.CamF.ROI_film(4);
    m = cam_spec.CamF.ROI_film(3);        
    hImage = image(zeros(n, m, 1));
    FFig.Position = [10 150 m  n];
    preview(F_vid, hImage);
    handles.Cams.FFig = FFig; %save the figure handle
end

% Update handles structure
guidata(hObject, handles);





% --- Executes on button press in pushbutton4. LASER ON
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
daqreset
s_vid_light = daq.createSession('ni');
s_vid_light.Rate = 100;

% get the laser intensity
laser_idx = get(handles.listbox5, 'Value');
ALL_intensities = [1 1.2 1.4 1.6 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0 6.0 7.0 8.0 9.0];
intensity = ALL_intensities(laser_idx(1));
% intensity = 1.4 ; 
light_length = 0.5;
on = intensity*ones(s_vid_light.Rate*light_length,1); % ends on
addAnalogOutputChannel(s_vid_light,'Dev1', 'ao0', 'Voltage');
queueOutputData(s_vid_light, on) 
startBackground(s_vid_light); 

fprintf(['\n Laser ON : ' num2str(intensity) 'V \n']) 
set(handles.text6, 'String', 'LASER ON')

% --- Executes on button press in pushbutton5 LASER OFF.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
daqreset
s_vid_light = daq.createSession('ni');
s_vid_light.Rate = 100;

addAnalogOutputChannel(s_vid_light,'Dev1', 'ao0', 'Voltage');
queueOutputData(s_vid_light, [0 0 0 0 0]') 
startBackground(s_vid_light);    
fprintf(' Laser OFF \n') 
set(handles.text6, 'String', 'LASER OFF')


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

 
% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% disp(handles.figure1.Position)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton7. CAMERA RESET
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton8_Callback(hObject, eventdata, handles)
imaqreset
 warning('off', 'imaq:gentl:hardwareTriggerTriggerModeOff');
 warning('off', 'imaq:gentl:adaptorSetROIModified');

% --- Executes on button press in pushbutton8. CLOSE PREVIEWS
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.figure1, 'HandleVisibility', 'off');
close all;
set(handles.figure1, 'HandleVisibility', 'on');


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for icam = 1:6
    Cams.([Alphabet(icam) '_Cam']) = handles.Cams.([Alphabet(icam) 'Fig']).Position;
end
save('Film Previews Config', 'Cams')
disp(handles.Cams)

% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox5


% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in PControlButton.
function PControlButton_Callback(hObject, eventdata, handles)
% hObject    handle to PControlButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PControl


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CamAbutton.
function CamAbutton_Callback(hObject, eventdata, handles)
% hObject    handle to CamAbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CamAbutton


% --- Executes on button press in CamBbutton.
function CamBbutton_Callback(hObject, eventdata, handles)
% hObject    handle to CamBbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CamBbutton


% --- Executes on button press in CamCbutton.
function CamCbutton_Callback(hObject, eventdata, handles)
% hObject    handle to CamCbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CamCbutton


% --- Executes on button press in CamDbutton.
function CamDbutton_Callback(hObject, eventdata, handles)
% hObject    handle to CamDbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CamDbutton


% --- Executes on button press in CamEbutton.
function CamEbutton_Callback(hObject, eventdata, handles)
% hObject    handle to CamEbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CamEbutton


% --- Executes on button press in CamFbutton.
function CamFbutton_Callback(hObject, eventdata, handles)
% hObject    handle to CamFbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CamFbutton


% --- Executes on button press in SelectAllCamsButton.
function SelectAllCamsButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectAllCamsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.CamAbutton, 'Value', 1)
set(handles.CamBbutton, 'Value', 1)
set(handles.CamCbutton, 'Value', 1)
set(handles.CamDbutton, 'Value', 1)
set(handles.CamEbutton, 'Value', 1)
set(handles.CamFbutton, 'Value', 1)

% --- Executes on button press in SelectNOcameras.
function SelectNOcameras_Callback(hObject, eventdata, handles)
% hObject    handle to SelectNOcameras (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.CamAbutton, 'Value', 0)
set(handles.CamBbutton, 'Value', 0)
set(handles.CamCbutton, 'Value', 0)
set(handles.CamDbutton, 'Value', 0)
set(handles.CamEbutton, 'Value', 0)
set(handles.CamFbutton, 'Value', 0)

function cam = getCamList(hObject, ~, handles)
camNames = {'CamAbutton', 'CamBbutton', 'CamCbutton',...
            'CamDbutton', 'CamEbutton', 'CamFbutton'};
for icam = 1:length(camNames)
   cam(icam) = get(handles.(camNames{icam}), 'Value') ;
end


% --- Executes on button press in pushbutton15. CLC clear command window
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
