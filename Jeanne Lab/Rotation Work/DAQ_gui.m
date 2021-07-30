function varargout = DAQ_gui(varargin)
% DAQ_GUI MATLAB code for DAQ_gui.fig
%      DAQ_GUI, by itself, creates a new DAQ_GUI or raises the existing
%      singleton*.
%
%      H = DAQ_GUI returns the handle to a new DAQ_GUI or the handle to
%      the existing singleton*.
%
%      DAQ_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DAQ_GUI.M with the given input arguments.
%
%      DAQ_GUI('Property','Value',...) creates a new DAQ_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DAQ_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DAQ_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DAQ_gui

% Last Modified by GUIDE v2.5 10-Sep-2020 08:44:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DAQ_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @DAQ_gui_OutputFcn, ...
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


% --- Executes just before DAQ_gui is made visible.
function DAQ_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DAQ_gui (see VARARGIN)
% set gui position
% handles.figure1.Position = [256 3.6923 124.4000 42.3077];
handles.volts = 5;
handles.fs = 10E3;

% Load DAQ session:
s = daq.createSession('ni');
s.addAnalogOutputChannel('Dev1',[0:3] , 'Voltage'); %LED,ODOR,SHOCK,CAM
s.Rate = handles.fs;
disp('DAQ connected')
handles.s = s;

% Choose default command line output for DAQ_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DAQ_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DAQ_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in COMBO_start.
function COMBO_start_Callback(hObject, eventdata, handles)
% hObject    handle to COMBO_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

volts = handles.volts;
combo = false(1,3);
stim_length = nan(1,3);
fs = handles.fs;

% print selection:
disp('Selected stimulus:')
% get status of each toggle:
if get(handles.LED_combo, 'value')==true
    stim_length(1) = str2double(get(handles.LED_length, 'String'));
    combo(1) = true;
    ROI(1) = round(fs*stim_length(1));
    disp('LED')
end
if get(handles.ODOR_combo, 'value')==true
    stim_length(2) = str2double(get(handles.ODOR_length, 'String'));
    combo(2) = true;
    ROI(2) = round(fs*stim_length(2));
    disp('Odor')
end   
if get(handles.SHOCK_combo, 'value')==true
    stim_length(3) = str2double(get(handles.SHOCK_length, 'String'));
    combo(3) = true;
    ROI(3) = round(fs*stim_length(3));
    disp('Shock')
end

% Set the analogue output vectors
total_l = max(stim_length);
vector_l = zeros(round(handles.fs*total_l),1);

all = repmat(vector_l,1,4);

for ii = 1:length(combo)
    if combo(ii)==true
        all(1:ROI(ii),ii) = volts;
    end
end

all(end,:) = zeros(4,1);
% figure; plot(all); ylim([-1,6])
queueOutputData(handles.s, all) 
startBackground(handles.s) 
disp('Group stimulus start')
    



% WORKING HERE 9.10.2020

% ODOR_length = str2double(get(handles.ODOR_length, 'String'));
% vector_l = round(ODOR_length*handles.fs);
% volts = handles.volts;
% outsig = volts*ones(vector_l,1);
% outsig(end,1) = 0;
% 
% % channels for LED, odor, shock, cam:
% dummy = zeros(vector_l,1);
% data = [dummy, outsig, dummy, dummy];
% 
% %load data to DAQ:
% queueOutputData(handles.s, data) 
% startBackground(handles.s) 
% disp('Odor start')





% --- Executes on button press in LED_combo.
function LED_combo_Callback(hObject, eventdata, handles)
% hObject    handle to LED_combo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LED_combo


% --- Executes on button press in ODOR_combo.
function ODOR_combo_Callback(hObject, eventdata, handles)
% hObject    handle to ODOR_combo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ODOR_combo


% --- Executes on button press in SHOCK_combo.
function SHOCK_combo_Callback(hObject, eventdata, handles)
% hObject    handle to SHOCK_combo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SHOCK_combo



function SHOCK_length_Callback(hObject, eventdata, handles)
% hObject    handle to SHOCK_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SHOCK_length as text
%        str2double(get(hObject,'String')) returns contents of SHOCK_length as a double


% --- Executes during object creation, after setting all properties.
function SHOCK_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SHOCK_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SHOCK_start.
function SHOCK_start_Callback(hObject, eventdata, handles)
% hObject    handle to SHOCK_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SHOCK_length = str2double(get(handles.SHOCK_length, 'String'));
vector_l = round(SHOCK_length*handles.fs);
volts = handles.volts;
outsig = volts*ones(vector_l,1);
outsig(end,1) = 0;

% channels for LED, odor, shock, cam:
dummy = zeros(vector_l,1);
data = [dummy, dummy, outsig, dummy];

%load data to DAQ:
queueOutputData(handles.s, data) 
startBackground(handles.s) 
disp('Shock start')

% --- Executes on button press in ODOR_start.
function ODOR_start_Callback(hObject, eventdata, handles)
% hObject    handle to ODOR_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ODOR_length = str2double(get(handles.ODOR_length, 'String'));
vector_l = round(ODOR_length*handles.fs);
volts = handles.volts;
outsig = volts*ones(vector_l,1);
outsig(end,1) = 0;

% channels for LED, odor, shock, cam:
dummy = zeros(vector_l,1);
data = [dummy, outsig, dummy, dummy];

%load data to DAQ:
queueOutputData(handles.s, data) 
startBackground(handles.s) 
disp('Odor start')


function ODOR_length_Callback(hObject, eventdata, handles)
% hObject    handle to ODOR_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ODOR_length as text
%        str2double(get(hObject,'String')) returns contents of ODOR_length as a double


% --- Executes during object creation, after setting all properties.
function ODOR_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ODOR_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LED_length_Callback(hObject, eventdata, handles)
% hObject    handle to LED_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LED_length as text
%        str2double(get(hObject,'String')) returns contents of LED_length as a double


% --- Executes during object creation, after setting all properties.
function LED_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LED_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LED_on.
function LED_on_Callback(hObject, eventdata, handles)
% hObject    handle to LED_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vector_l = round(0.02*handles.fs);

LED_volt = handles.volts;
outsig = LED_volt*ones(vector_l,1);
% outsig(end,1) = 0;

% dummy channels for odor, shock, cam:
dummy = zeros(vector_l,1);
data = [outsig, dummy, dummy, dummy];

%load data to LED:
queueOutputData(handles.s, data) 
startBackground(handles.s) 
disp('LED on')


% --- Executes on button press in LED_off.
function LED_off_Callback(hObject, eventdata, handles)
% hObject    handle to LED_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vector_l = round(0.01*handles.fs);

LED_volt = handles.volts;
outsig = LED_volt*zeros(vector_l,1);

% dummy channels for odor, shock, cam:
dummy = zeros(vector_l,1);
data = [outsig, dummy, dummy, dummy];

%load data to LED:
queueOutputData(handles.s, data) 
startBackground(handles.s) 
disp('LED off')


% --- Executes on button press in LED_start.
function LED_start_Callback(hObject, eventdata, handles)
% hObject    handle to LED_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LED_length = str2double(get(handles.LED_length, 'String'));

vector_l = round(LED_length*handles.fs);

LED_volt = handles.volts;
outsig = LED_volt*ones(vector_l,1);
outsig(end,1) = 0;

% dummy channels for odor, shock, cam:
dummy = zeros(vector_l,1);
data = [outsig, dummy, dummy, dummy];

%load data to LED:
queueOutputData(handles.s, data) 
startBackground(handles.s) 
disp('LED start')




















