function varargout = behavior_gui(varargin)
% BEHAVIOR_GUI MATLAB code for behavior_gui.fig
%      BEHAVIOR_GUI, by itself, creates a new BEHAVIOR_GUI or raises the existing
%      singleton*.
%  
%      H = BEHAVIOR_GUI returns the handle to a new BEHAVIOR_GUI or the handle to
%      the existing singleton*.
%
%      BEHAVIOR_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEHAVIOR_GUI.M with the given input arguments.
%
%      BEHAVIOR_GUI('Property','Value',...) creates a new BEHAVIOR_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before behavior_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to behavior_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help behavior_gui

% Last Modified by GUIDE v2.5 08-Apr-2020 11:26:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @behavior_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @behavior_gui_OutputFcn, ...
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


% --- Executes just before behavior_gui is made visible.
function behavior_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to behavior_gui (see VARARGIN)

% Choose default command line output for behavior_gui
handles.output = hObject;
handles.num.reps = 3; %Default values -- updates later automatically
handles.num.conds = 28;
warning('off')
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes behavior_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Executes on KEY CODING with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Character
    case 'a'
        disp('other')
        OtherButton_Callback(hObject, eventdata, handles)
    case 's'
        disp('grooming')
        GroomingButton_Callback(hObject, eventdata, handles)
    case 'd'
        disp('walking')
        WalkingButton_Callback(hObject, eventdata, handles)
    case 'c'
        disp('stationary')
        StationaryButton_Callback(hObject, eventdata, handles)
    case 'k'
        disp('swing')
        SwingButton_Callback(hObject, eventdata, handles)
    case 'm'
        disp('stance')
        StanceButton_Callback(hObject, eventdata, handles)
    case 'n'
        disp('Next Video')
        NextButton_Callback(hObject, eventdata, handles)
%     case 'ButtonDownFcn'
%         disp('Next Video')
%         NextButton_Callback(hObject, eventdata, handles)
    case 'v'
        disp('Previous Video') 
        PreviousButton_Callback(hObject, eventdata, handles)
    
    case 'A'
        set(handles.cameralistbox, 'Value', 1)
        figure(handles.figure1)
    case 'z'
        NextButton_Callback(hObject, eventdata, handles)
    case 'Z'
        TrackingFrameUpButton_Callback(hObject, eventdata, handles)
    case 'X'
        TrackingFrameDownButton_Callback(hObject, eventdata, handles)
    case 'B'
        set(handles.cameralistbox, 'Value', 2)  
        figure(handles.figure1)
    case 'C'
        set(handles.cameralistbox, 'Value', 3) 
        figure(handles.figure1)
    case 'D'
        set(handles.cameralistbox, 'Value', 4)
        figure(handles.figure1)
    case 'E'
        set(handles.cameralistbox, 'Value', 5)
        figure(handles.figure1)
    case 'F'
        set(handles.cameralistbox, 'Value', 6)
        figure(handles.figure1)
end
if strcmpi(eventdata.Key,'return')
        disp('Replay Video') 
        PlayButton_Callback(hObject, eventdata, handles)
end



% --- Outputs from this function are returned to the command line.
function varargout = behavior_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------- DONT EDIT NORTH OF HERE ------------------------




% --- Executes on button press in StationaryButton.
function StationaryButton_Callback(hObject, eventdata, handles)
% hObject    handle to StationaryButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.CurrentBehaviorTag, 'String', 'Stationary')

% --- Executes on button press in WalkingButton.
function WalkingButton_Callback(hObject, eventdata, handles)
% hObject    handle to WalkingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.CurrentBehaviorTag, 'String', 'Walking')


% --- Executes on button press in GroomingButton.
function GroomingButton_Callback(hObject, eventdata, handles)
% hObject    handle to GroomingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.CurrentBehaviorTag, 'String', 'Grooming')

% --- Executes on button press in OtherButton.
function OtherButton_Callback(hObject, eventdata, handles)
% hObject    handle to OtherButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.CurrentBehaviorTag, 'String', 'Other')


% --- Executes on button press in SelectStructButton.
function SelectStructButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectStructButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ------------- SELECT CROSS FROM LIST ----------------
handles.excelstartpoint = 470; %TODO
strtpoint = handles.excelstartpoint; %line in excel file to start on
[excelfile, Excel] = load_flysummary;
% Pull out the structure names from the file:
structure_names.excelfile = excelfile(strtpoint:end, Excel.structure);
ind = 1;
for ii = 1:length(structure_names.excelfile)
    if ischar(structure_names.excelfile{ii})
        structure_names.text{ind} = (structure_names.excelfile{ii});
        ind = ind+1;
    end
end
%find the unique structure names:
structure_names.unique = unique(structure_names.text);
indx = listdlg('ListString', structure_names.unique, 'SelectionMode', 'Single', 'ListSize', [300, 400]);

structure_name = structure_names.unique{indx};
handles.structure_name = structure_name;
% update structure name:
set(handles.StructureName, 'String', {' '; structure_name})

location = find(strcmpi(structure_name, excelfile(strtpoint:end,Excel.new_struct_name)))+strtpoint-1;
% save the fly information to the handle
handles.excelfile = excelfile(location,:);
bb = 1;
for tt = 1:length(location)
    a = cell2mat(excelfile(location(tt), Excel.structurenum));
    if ischar(a) %no value for the structure
    else %has a number in the structure
        flynum(bb) = a; %find total number of flies
        Location.rownum(bb) = location(tt); %get the rownumber for each fly
        bb = bb+1; %increase the index for the number of flies
    end
end

num.flies = max(flynum);
handles.num.flies = num.flies;
  
% Error Check
if sum(num.flies) == 0
    warndlg('Check name of flies in Excel File, no matching files found')
    return
end

% -------- load in parameter variables ----------
t_file_root = 'G:\My Drive\Evyn\Data\FicTrac Raw Data\'; %TODO

t_folder = excelfile{location(1),Excel.date};
t_flyNum = excelfile{location(1),Excel.flynum};
t_flyID = generate_flyID(t_folder, t_flyNum);
t_fullfile = [t_file_root, t_folder, '\Fly ', t_flyNum, '\FicTrac Data\', t_flyID '.mat'];
% e.g. G:\My Drive\Evyn\Data\FicTrac Raw Data\8.13.19\Fly 1_3\FicTrac Data\08132019_fly1_3.mat'
% load the matlab data file:
MTLB_data = load(t_fullfile);

% extract the parameter data:
param = MTLB_data.fly.param;
handles.num.conds = param.num_conds;
handles.num.reps = param.num_reps;
% one label per condition:
for ii = 1:param.num_conds
    labels{ii} = param.conds_matrix(ii).label;
end
handles.conds_labels = labels;
% -----------------------------------------------

% % set the folder directory that contains the video labels (e.g. behavior
% class and joint angles etc)
% folderpath = uigetdir('','Select Behavior Folder');
% handles.folderpath = folderpath;
folderpath = 'C:\matlabroot\behavior class'; %TODO


%load or create behavior structure
try load([folderpath, '\' structure_name ' behavior class'], 'group')
    fprintf('\n Loaded behavior record\n')
    handles.group = group;
catch
    fprintf('\n No existing behavior record found\n')
    group = struct;
    for kk = 1:num.flies
        for icond = 1:handles.num.conds
            for irep = 1:handles.num.reps
                group(kk).behavior{icond, irep} = '-';
                group(kk).phase{icond, irep} = '-';
            end
        end
    end
   handles.group = group; 
end

%load or create response structure
try load([folderpath, '\' structure_name ' response data'], 'ResponseData')
    fprintf(' Loaded response record\n')
    handles.response_data = ResponseData;
catch
    fprintf(' No existing responses record found\n')
    response_data = struct;
    for kk = 1:num.flies
        for icond = 1:handles.num.conds
            for irep = 1:handles.num.reps
                response_data(kk).behavior{icond, irep} = '-';
            end
        end
    end
   handles.response_data = response_data; 
end

%load or create Joint Angle structures
try load([folderpath, '\' structure_name ' joint angle data'], 'JA_Data')
    fprintf(' Loaded joint angles \n')
    handles.JA_Data = JA_Data;
catch
    fprintf(' No existing joint angle record found\n')
    JA_Data = struct;
    for kk = 1:num.flies
        for icond = 1:handles.num.conds
            for irep = 1:handles.num.reps
                JA_Data(kk).angles(icond, irep).data = zeros(4,2);
                JA_Data(kk).CoFe{icond, irep} = '-';
                JA_Data(kk).FeTi{icond, irep} = '-';
                % add structure for series joint tracking
                
                JA_Data(kk).FeTi{icond, irep} = '-';
            end
        end
    end
   handles.JA_Data = JA_Data; 
end

try load([folderpath, '\' structure_name ' tracking data'], 'tracking')
    fprintf(' Loaded tracking struct \n')
    handles.tracking = tracking;
catch
    fprintf(' No existing tracking record found \n')
end

handles.folderpath = folderpath;
set(handles.cameralistbox, 'Value', 3)
set(handles.fpslistbox, 'Value', 3)
guidata(hObject, handles);

% ------------- DISPLAY DATA FOR FLY #1 IN STRUCTURE ------------------
%SET COND=1 AND REP=1
ResetRepNCond(handles)
ResetFlyNum(handles)

UpdateLabels(handles);


function ResetFlyNum(handles)
set(handles.FlyNumEdit, 'String', '1')


function UpdateLabels(handles)
ShowStateLabels(handles)
% find and set fly ID
[~, cond, flynum] = getflyinfo(handles);

dateflynum = handles.excelfile{flynum,2};
datestring = dateconverter(handles.excelfile{flynum,1});
fly_ID = [datestring '_fly' dateflynum];
set(handles.FlyIDEdit, 'String', fly_ID)
% find and set the stimulus info
% cond_label = getcond(cond);
cond_label = handles.conds_labels{cond};
set(handles.StimulusEdit, 'String', cond_label)


function ShowStateLabels(handles)
[rep, cond, flynum] = getflyinfo(handles);
phase = handles.group(flynum).phase{cond, rep};
behavior = handles.group(flynum).behavior{cond, rep};
response = handles.response_data(flynum).behavior{cond, rep};
CoFe = handles.JA_Data(flynum).CoFe{cond, rep};
FeTi = handles.JA_Data(flynum).FeTi{cond, rep};
%display current behavior and phase info
set(handles.CurrentBehaviorTag, 'String', behavior)
set(handles.PhaseTag, 'String', phase)
set(handles.CurrentResponseLabel, 'String', response)
set(handles.CoFe_joint_edit, 'String', CoFe)
set(handles.FeTi_joint_edit, 'String', FeTi)


function handles = saveCurrentLabels(handles)
handles = SavePhase(handles);
handles = SaveBehavior(handles);
handles = SaveResponse(handles);



% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = saveCurrentLabels(handles);
% handles = SaveBehavior(handles);
% handles = SavePhase(handles);
group = handles.group;
save([handles.folderpath '\' handles.structure_name ' behavior class'], 'group')


function [rep, cond, flynum] = getflyinfo(handles)
rep = str2double(get(handles.RepEdit, 'String'));
cond = str2double(get(handles.ConditionEdit, 'String'));
flynum = str2double(get(handles.FlyNumEdit, 'String'));


function handles = SaveBehavior(handles)
[rep, cond, flynum] = getflyinfo(handles);
newbehavior = get(handles.CurrentBehaviorTag, 'String');
handles.group(flynum).behavior{cond, rep} = newbehavior;


function handles = SaveResponse(handles)
[rep, cond, flynum] = getflyinfo(handles);
newresponse = get(handles.CurrentResponseLabel, 'String');
handles.response_data(flynum).behavior{cond, rep} = newresponse;


function handles = SavePhase(handles)
[rep, cond, flynum] = getflyinfo(handles);
newphase = get(handles.PhaseTag, 'String');
handles.group(flynum).phase{cond, rep} = newphase;



function ResetRepNCond(handles)
handles = saveCurrentLabels(handles);
set(handles.RepEdit, 'String', '1')
set(handles.ConditionEdit, 'String', '1')


% --- Executes on button press in SwingButton.
function SwingButton_Callback(hObject, eventdata, handles)
% hObject    handle to SwingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.PhaseTag, 'String', 'Swing')

% --- Executes on button press in StanceButton.
function StanceButton_Callback(hObject, eventdata, handles)
% hObject    handle to StanceButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.PhaseTag, 'String', 'Stance')

% --- Executes on button press in PreviousButton.
function PreviousButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreviousButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try close(handles.fig); catch; end
handles = SaveBehavior(handles);
handles = SavePhase(handles);
handles = SaveResponse(handles);

guidata(hObject, handles);
%update the condition and rep numbers
[rep, cond, flynum] = getflyinfo(handles);
if cond == 1 && rep == 1
    if flynum > 1
        cond = handles.num.conds;
        rep = handles.num.reps;
        flynum = flynum-1;
    else 
        warndlg('No data before this')
        return
    end
elseif cond == 1 && rep > 1
    cond = handles.num.conds;
    rep = rep-1;
else
    cond = cond-1;
end

guidata(hObject, handles);
set(handles.RepEdit, 'String', rep)
set(handles.ConditionEdit, 'String', cond)
set(handles.FlyNumEdit, 'String', flynum)
UpdateLabels(handles)

switch handles.JA_radiobutton.Value
    case 0
        PlayButton_Callback(hObject, eventdata, handles)
    case 1
        Display_JA_image_button_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in PlayButton.
function PlayButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try close(handles.fig); catch; end
[rep, cond, flynum] = getflyinfo(handles);
fly_date = handles.excelfile{flynum, 1}; 
fly_num = handles.excelfile{flynum, 2};
fly_ID = get(handles.FlyIDEdit, 'String');  
camlist = get(handles.cameralistbox, 'String');
Cam = camlist{get(handles.cameralistbox, 'Value')};
fpslist = get(handles.fpslistbox, 'String');
fps = str2double(fpslist{get(handles.fpslistbox, 'Value')});
starttime = str2double(get(handles.StartTimeEdit, 'String'));
stoptime = str2double(get(handles.StopTimeEdit, 'String'));

% Play video
drive_root = get(handles.filerootedit, 'String');
vid_root = [drive_root char(fly_date) '\Fly ' char(fly_num) '\Raw Video\'];

vid_name = [vid_root, fly_ID, ' R' num2str(rep), 'C', num2str(cond), ' Cam-' Cam '*'];
FilePath = dir(vid_name);
if isempty(FilePath)
    warndlg('Check video file path')
    return
end
vid_name = FilePath.name;
FrameRate = fps;
handles.fig = vidplayback([vid_root, vid_name], FrameRate, starttime, stoptime);
% disp(vid_name) 
figure(handles.figure1)
guidata(hObject, handles);

% --- Executes on button press in NextButton.
function NextButton_Callback(hObject, eventdata, handles)
% hObject    handle to NextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try close(handles.fig); catch; end
handles = SaveBehavior(handles);
handles = SavePhase(handles);
handles = SaveResponse(handles);

guidata(hObject, handles);
%update the condition and rep numbers
[rep, cond, flynum] = getflyinfo(handles);
if cond == handles.num.conds && rep < handles.num.reps
    cond = 1;
    rep = rep+1;
elseif cond == handles.num.conds && rep == handles.num.reps
    if flynum < handles.num.flies
        flynum = flynum+1;
        cond = 1;
        rep = 1;
    elseif flynum == handles.num.flies
        warndlg('Finished Flies')
        return
    end
else 
    cond = cond+1;
end
set(handles.RepEdit, 'String', rep)
set(handles.ConditionEdit, 'String', cond)
set(handles.FlyNumEdit, 'String', flynum)
UpdateLabels(handles)

switch handles.JA_radiobutton.Value
    case 0 
       PlayButton_Callback(hObject, eventdata, handles)
    case 1
       Display_JA_image_button_Callback(hObject, eventdata, handles)
end
    


function ConditionEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ConditionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ConditionEdit as text
%        str2double(get(hObject,'String')) returns contents of ConditionEdit as a double


% --- Executes during object creation, after setting all properties.
function ConditionEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ConditionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RepEdit_Callback(hObject, eventdata, handles)
% hObject    handle to RepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RepEdit as text
%        str2double(get(hObject,'String')) returns contents of RepEdit as a double


% --- Executes during object creation, after setting all properties.
function RepEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RepEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FlyIDEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FlyIDEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FlyIDEdit as text
%        str2double(get(hObject,'String')) returns contents of FlyIDEdit as a double


% --- Executes during object creation, after setting all properties.
function FlyIDEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FlyIDEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StimulusEdit_Callback(hObject, eventdata, handles)
% hObject    handle to StimulusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StimulusEdit as text
%        str2double(get(hObject,'String')) returns contents of StimulusEdit as a double


% --- Executes during object creation, after setting all properties.
function StimulusEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StimulusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FlyNumEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FlyNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FlyNumEdit as text
%        str2double(get(hObject,'String')) returns contents of FlyNumEdit as a double


% --- Executes during object creation, after setting all properties.
function FlyNumEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FlyNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FlyNumIncreaseButton.
function FlyNumIncreaseButton_Callback(hObject, eventdata, handles)
% hObject    handle to FlyNumIncreaseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = saveCurrentLabels(handles);
oldflynum = str2double(get(handles.FlyNumEdit, 'String'));
switch oldflynum
    case num2cell(1:handles.num.flies-1) %All fly numbers except the last in group
        set(handles.FlyNumEdit, 'String', num2str(oldflynum+1))
        UpdateLabels(handles)
    case handles.num.flies %last fly in group==reset to 1
        set(handles.FlyNumEdit, 'String', '1')
%         ResetFlyNum(handles)
%         ResetRepNCond(handles)
        UpdateLabels(handles)
end
guidata(hObject, handles);

% --- Executes on button press in FlyNumDownButton.
function FlyNumDownButton_Callback(hObject, eventdata, handles)
% hObject    handle to FlyNumDownButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = saveCurrentLabels(handles);
oldflynum = str2double(get(handles.FlyNumEdit, 'String'));
switch oldflynum
    case num2cell(2:handles.num.flies) %All fly numbers except the first
        set(handles.FlyNumEdit, 'String', num2str(oldflynum-1))
        UpdateLabels(handles)
    case 1 %first fly in group==reset to biggest number
        set(handles.FlyNumEdit, 'String', num2str(handles.num.flies))
%         ResetFlyNum(handles)
%         ResetRepNCond(handles)
        UpdateLabels(handles)
end
guidata(hObject, handles);


% --- Executes on button press in RepDownButton.
function RepDownButton_Callback(hObject, eventdata, handles)
% hObject    handle to RepDownButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = saveCurrentLabels(handles);
oldnum = str2double(get(handles.RepEdit, 'String'));
switch oldnum
    case num2cell(2:handles.num.reps)
        set(handles.RepEdit, 'String', num2str(oldnum-1))
        UpdateLabels(handles)
    case 1
        set(handles.RepEdit, 'String', num2str(handles.num.reps))
%         ResetFlyNum(handles)
%         ResetRepNCond(handles)
        UpdateLabels(handles)
end
guidata(hObject, handles);

% --- Executes on button press in RepUpButton.
function RepUpButton_Callback(hObject, eventdata, handles)
% hObject    handle to RepUpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = saveCurrentLabels(handles);
oldnum = str2double(get(handles.RepEdit, 'String'));
switch oldnum
    case num2cell(1:handles.num.reps-1)
        set(handles.RepEdit, 'String', num2str(oldnum+1))
        UpdateLabels(handles)
    case handles.num.reps
        set(handles.RepEdit, 'String', '1')
%         ResetFlyNum(handles)
%         ResetRepNCond(handles)
        UpdateLabels(handles)
end
guidata(hObject, handles);

% --- Executes on button press in CondDownButton.
function CondDownButton_Callback(hObject, eventdata, handles)
% hObject    handle to CondDownButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = saveCurrentLabels(handles);
oldnum = str2double(get(handles.ConditionEdit, 'String'));
switch oldnum
    case num2cell(2:handles.num.conds)
        set(handles.ConditionEdit, 'String', num2str(oldnum-1))
        UpdateLabels(handles)
    case 1
        set(handles.ConditionEdit, 'String', num2str(handles.num.conds))
%         ResetFlyNum(handles)
%         ResetRepNCond(handles)
        UpdateLabels(handles)
end
guidata(hObject, handles);

% --- Executes on button press in CondUpButton.
function CondUpButton_Callback(hObject, eventdata, handles)
% hObject    handle to CondUpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = saveCurrentLabels(handles);
oldnum = str2double(get(handles.ConditionEdit, 'String'));
switch oldnum
    case num2cell(1:handles.num.conds-1)
        set(handles.ConditionEdit, 'String', num2str(oldnum+1))
        UpdateLabels(handles)
    case handles.num.conds
        set(handles.ConditionEdit, 'String', '1')
%         ResetFlyNum(handles)
%         ResetRepNCond(handles)
        UpdateLabels(handles)
end
guidata(hObject, handles);

% --- Executes on button press in FlySetButton.
function FlySetButton_Callback(hObject, eventdata, handles)
% hObject    handle to FlySetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateLabels(handles)
guidata(hObject, handles);
 


function CurrentBehaviorTag_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentBehaviorTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrentBehaviorTag as text
%        str2double(get(hObject,'String')) returns contents of CurrentBehaviorTag as a double


% --- Executes during object creation, after setting all properties.
function CurrentBehaviorTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentBehaviorTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PhaseTag_Callback(hObject, eventdata, handles)
% hObject    handle to PhaseTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PhaseTag as text
%        str2double(get(hObject,'String')) returns contents of PhaseTag as a double


% --- Executes during object creation, after setting all properties.
function PhaseTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PhaseTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cameralistbox.
function cameralistbox_Callback(hObject, eventdata, handles)
% hObject    handle to cameralistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cameralistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cameralistbox


% --- Executes during object creation, after setting all properties.
function cameralistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cameralistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fpslistbox.
function fpslistbox_Callback(hObject, eventdata, handles)
% hObject    handle to fpslistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fpslistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fpslistbox


% --- Executes during object creation, after setting all properties.
function fpslistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fpslistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StartTimeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to StartTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setTotalFrames(handles)


% Hints: get(hObject,'String') returns contents of StartTimeEdit as text
%        str2double(get(hObject,'String')) returns contents of StartTimeEdit as a double


% --- Executes during object creation, after setting all properties.
function StartTimeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StartTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function StopTimeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to StopTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StopTimeEdit as text
%        str2double(get(hObject,'String')) returns contents of StopTimeEdit as a double
setTotalFrames(handles)


% --- Executes during object creation, after setting all properties.
function StopTimeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StopTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CloseAllButton.
function CloseAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to CloseAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.figure1, 'HandleVisibility', 'off');
close all; clc
set(handles.figure1, 'HandleVisibility', 'on');



% --- Executes on button press in SaveResponsesButton.
function SaveResponsesButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveResponsesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = SaveResponse(handles);
ResponseData = handles.response_data;
save([handles.folderpath '\' handles.structure_name ' response data'], 'ResponseData')


function responseEdit1_Callback(hObject, eventdata, handles)
% hObject    handle to responseEdit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of responseEdit1 as text
%        str2double(get(hObject,'String')) returns contents of responseEdit1 as a double


% --- Executes during object creation, after setting all properties.
function responseEdit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to responseEdit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in responseButton1.
function responseButton1_Callback(hObject, eventdata, handles)
% hObject    handle to responseButton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
behavior_response = get(handles.responseEdit1, 'String');
set(handles.CurrentResponseLabel, 'String', behavior_response)


function responseEdit2_Callback(hObject, eventdata, handles)
% hObject    handle to responseEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of responseEdit2 as text
%        str2double(get(hObject,'String')) returns contents of responseEdit2 as a double


% --- Executes during object creation, after setting all properties.
function responseEdit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to responseEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in responseButton2.
function responseButton2_Callback(hObject, eventdata, handles)
% hObject    handle to responseButton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
behavior_response = get(handles.responseEdit2, 'String');
set(handles.CurrentResponseLabel, 'String', behavior_response)


function responseEdit3_Callback(hObject, eventdata, handles)
% hObject    handle to responseEdit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of responseEdit3 as text
%        str2double(get(hObject,'String')) returns contents of responseEdit3 as a double


% --- Executes during object creation, after setting all properties.
function responseEdit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to responseEdit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in responseButton3.
function responseButton3_Callback(hObject, eventdata, handles)
% hObject    handle to responseButton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
behavior_response = get(handles.responseEdit3, 'String');
set(handles.CurrentResponseLabel, 'String', behavior_response)


function responseEdit4_Callback(hObject, eventdata, handles)
% hObject    handle to responseEdit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of responseEdit4 as text
%        str2double(get(hObject,'String')) returns contents of responseEdit4 as a double


% --- Executes during object creation, after setting all properties.
function responseEdit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to responseEdit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in responseButton4.
function responseButton4_Callback(hObject, eventdata, handles)
% hObject    handle to responseButton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
behavior_response = get(handles.responseEdit4, 'String');
set(handles.CurrentResponseLabel, 'String', behavior_response)


function responseEdit5_Callback(hObject, eventdata, handles)
% hObject    handle to responseEdit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of responseEdit5 as text
%        str2double(get(hObject,'String')) returns contents of responseEdit5 as a double


% --- Executes during object creation, after setting all properties.
function responseEdit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to responseEdit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in responseButton5.
function responseButton5_Callback(hObject, eventdata, handles)
% hObject    handle to responseButton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
behavior_response = get(handles.responseEdit5, 'String');
set(handles.CurrentResponseLabel, 'String', behavior_response)


function responseEdit6_Callback(hObject, eventdata, handles)
% hObject    handle to responseEdit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of responseEdit6 as text
%        str2double(get(hObject,'String')) returns contents of responseEdit6 as a double


% --- Executes during object creation, after setting all properties.
function responseEdit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to responseEdit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in responseButton6.
function responseButton6_Callback(hObject, eventdata, handles)
% hObject    handle to responseButton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

behavior_response = get(handles.responseEdit6, 'String');
set(handles.CurrentResponseLabel, 'String', behavior_response)

% % get current ID info
% [rep, cond, flynum] = getflyinfo(handles);
% % pull the behavior tag
% behavior_response = get(handles.responseEdit6, 'String');
% % save the behavior tag into the structure
% handles.response_data(flynum).behavior{cond, rep} = behavior_response;

% guidata(hObject, handles);


% --- Executes on button press in PreloadResponseButton.
function PreloadResponseButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreloadResponseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load the labels from the text file:
file_id = fopen([handles.folderpath '\Behavior Response Labels.txt']);
A = fscanf(file_id, '%s');
Cross = strsplit(A,','); 

% Insert the labels into the Edits:
for ii = 1:6
   label = ['responseEdit' num2str(ii)];
   set(handles.(label), 'String', Cross{ii}) 
end

guidata(hObject, handles);


% --- Executes on button press in SaveResponseLabelsButton.
function SaveResponseLabelsButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveResponseLabelsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles = SaveResponse(handles);
% ResponseData = handles.response_data;
% save([handles.structure_name ' response data'], 'ResponseData')



function joint_angle_time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to joint_angle_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of joint_angle_time_edit as text
%        str2double(get(hObject,'String')) returns contents of joint_angle_time_edit as a double


% --- Executes during object creation, after setting all properties.
function joint_angle_time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to joint_angle_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% LABEL FRAME BUTTON PRESS
% --- Executes on button press in Display_JA_image_button.
function Display_JA_image_button_Callback(hObject, eventdata, handles)
% hObject    handle to Display_JA_image_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% try close(handles.fig); catch; end

singleStatus = get(handles.JA_radiobutton, 'Value');
seriesStatus = get(handles.TrackingButton, 'Value');

if singleStatus == 1
    % tracking single data point
    FrameTime = str2double(get(handles.joint_angle_time_edit, 'String'));
elseif seriesStatus == 1
    % joint tracking longterm
    currentFrame = str2double(get(handles.TrackingCurrentFrame, 'String'));
    FrameTime = currentFrame/300;
else
    fprintf('\n Select Single or Series tracking \n')
    return
end  
handles = TrackFrame(handles, FrameTime);
% handles = TrackFrame(handles, FrameTime);
guidata(hObject, handles);




% 
% [rep, cond, flynum] = getflyinfo(handles);
% fly_date = handles.excelfile{flynum, 1}; 
% fly_num = handles.excelfile{flynum, 2};
% fly_ID = get(handles.FlyIDEdit, 'String');  
% camlist = get(handles.cameralistbox, 'String');
% % Cam = camlist{get(handles.cameralistbox, 'Value')};
% Cam = 'B';
% fpslist = get(handles.fpslistbox, 'String');
% fps = str2double(fpslist{get(handles.fpslistbox, 'Value')});
% starttime = str2double(get(handles.joint_angle_time_edit, 'String'));
% stoptime = starttime;
% 
% % Display Frame:
% drive_root = get(handles.filerootedit, 'String');
% vid_root = [drive_root char(fly_date) '\Fly ' char(fly_num) '\Raw Video\'];
% 
% vid_name = [vid_root, fly_ID, ' R' num2str(rep), 'C', num2str(cond), ' Cam-' Cam '*'];
% FilePath = dir(vid_name);
% vid_name = FilePath.name;
% disp([vid_root, vid_name])
% % FrameRate = fps;
% % handles.fig = vidplayback([vid_root, vid_name], FrameRate, starttime, stoptime);
% % disp(vid_name) 
% handles.fig = image_still([vid_root, vid_name], starttime, 1.8);
% 
% % % Select Points:
% % [x, y] = ginput(4);
% % hold all; scatter(x,y,'filled', 'y')
% 
% % Select Points:
% % [x, y] = ginput(4);
%  
% hold all;
% for ii = 1:4
% [x(ii), y(ii)] = crosshairs(1);
% scatter(x(ii),y(ii),'filled', 'y')
% end
% 
% % coxa femur:
% a = [x(1), y(1)];
% b = [x(2), y(2)];
% c = [x(3), y(3)];
% u(1) = (b(1)-a(1));
% u(2) = (b(2)-a(2));
% u(3) = 0;
% v(1) = (c(1)-b(1));
% v(2) = (c(2)-b(2));
% v(3) = 0;
% CoFe = 180-atan2d(norm(cross(u,v)),dot(u,v));
% set(handles.CoFe_joint_edit, 'string', CoFe)
% % disp(ThetaInDegrees)
% 
% % femur tibia:
% a = [x(2), y(2)];
% b = [x(3), y(3)];
% c = [x(4), y(4)];
% u(1) = (b(1)-a(1));
% u(2) = (b(2)-a(2));
% u(3) = 0;
% v(1) = (c(1)-b(1));
% v(2) = (c(2)-b(2));
% v(3) = 0;
% FeTi = 180-atan2d(norm(cross(u,v)),dot(u,v));
% set(handles.FeTi_joint_edit, 'string', FeTi)
% % disp(ThetaInDegrees)
% 
% % update data structure:
% handles.JA_Data(flynum).angles(cond,rep).data = [x',y'];
% handles.JA_Data(flynum).CoFe{cond,rep} = CoFe;
% handles.JA_Data(flynum).FeTi{cond,rep} = FeTi;
% 
% % disp([CoFe, FeTi])
% % disp(handles.JA_Data(flynum).angles(cond,rep).data)
% % disp(handles.JA_Data(flynum).FeTi{cond,rep})
% 
% figure(handles.figure1)
% guidata(hObject, handles);
% handles.filerootedit;


% --- Executes on button press in Save_JA_button.
function Save_JA_button_Callback(hObject, eventdata, handles)
% hObject    handle to Save_JA_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = saveCurrentLabels(handles);
% handles = SaveBehavior(handles);
% handles = SavePhase(handles);

JA_Data = handles.JA_Data;
save([handles.folderpath '\' handles.structure_name ' joint angle data'], 'JA_Data')
try
tracking = handles.tracking;
save([handles.folderpath '\' handles.structure_name ' tracking data'], 'tracking')
catch
    fprintf('\n No tracking data saved\n')
end


function filerootedit_Callback(hObject, eventdata, handles)
% hObject    handle to filerootedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filerootedit as text
%        str2double(get(hObject,'String')) returns contents of filerootedit as a double


% --- Executes during object creation, after setting all properties.
function filerootedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filerootedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CoFe_joint_edit_Callback(hObject, eventdata, handles)
% hObject    handle to CoFe_joint_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CoFe_joint_edit as text
%        str2double(get(hObject,'String')) returns contents of CoFe_joint_edit as a double


% --- Executes during object creation, after setting all properties.
function CoFe_joint_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CoFe_joint_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FeTi_joint_edit_Callback(hObject, eventdata, handles)
% hObject    handle to FeTi_joint_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FeTi_joint_edit as text
%        str2double(get(hObject,'String')) returns contents of FeTi_joint_edit as a double


% --- Executes during object creation, after setting all properties.
function FeTi_joint_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FeTi_joint_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in JA_preview_button.
function JA_preview_button_Callback(hObject, eventdata, handles)
% hObject    handle to JA_preview_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% try close(handles.fig); catch; end
[rep, cond, flynum] = getflyinfo(handles);
fly_date = handles.excelfile{flynum, 1}; 
fly_num = handles.excelfile{flynum, 2};
fly_ID = get(handles.FlyIDEdit, 'String');  
camlist = get(handles.cameralistbox, 'String');
% Cam = camlist{get(handles.cameralistbox, 'Value')};
Cam = 'B';
fpslist = get(handles.fpslistbox, 'String');
fps = str2double(fpslist{get(handles.fpslistbox, 'Value')});
starttime = str2double(get(handles.joint_angle_time_edit, 'String'));
stoptime = starttime;

% Display image:
drive_root = get(handles.filerootedit, 'String');
vid_root = [drive_root char(fly_date) '\Fly ' char(fly_num) '\Raw Video\'];

vid_name = [vid_root, fly_ID, ' R' num2str(rep), 'C', num2str(cond), ' Cam-' Cam '*'];
FilePath = dir(vid_name);
vid_name = FilePath.name;
% FrameRate = fps;
handles.fig = image_still([vid_root, vid_name], starttime, 1.8);
% handles.fig = vidplayback([vid_root, vid_name], FrameRate, starttime, stoptime, 1.2);

% plot the joint position points on the figure:
x = handles.JA_Data(flynum).angles(cond,rep).data(:,1);
y = handles.JA_Data(flynum).angles(cond,rep).data(:,2);
hold all; scatter(x,y,'filled', 'y')


% --- Executes on button press in next_JA_button.
function next_JA_button_Callback(hObject, eventdata, handles)
% hObject    handle to next_JA_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NextButton_Callback(hObject, eventdata, handles)
    


% --- Executes on button press in previous_JA_button.
function previous_JA_button_Callback(hObject, eventdata, handles)
% hObject    handle to previous_JA_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PreviousButton_Callback(hObject, eventdata, handles)

% --- Executes on button press in JA_radiobutton.
function JA_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to JA_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of JA_radiobutton


% --- Executes on button press in TrackingButton.
function TrackingButton_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TrackingBuildDataStructure_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of TrackingButton


function TrackingFrameStepInput_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingFrameStepInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setTotalFrames(handles)

% --- Executes during object creation, after setting all properties.
function TrackingFrameStepInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrackingFrameStepInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TrackingLegSelectionDroplist.
function TrackingLegSelectionDroplist_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingLegSelectionDroplist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TrackingLegSelectionDroplist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TrackingLegSelectionDroplist


% --- Executes during object creation, after setting all properties.
function TrackingLegSelectionDroplist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrackingLegSelectionDroplist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton39.
function pushbutton39_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton43.
function pushbutton43_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


function TrackingCurrentFrame_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingCurrentFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~, CoFe, FeTi] = getTrackingInfo(handles);

set(handles.CoFe_joint_edit, 'String', num2str(CoFe))
set(handles.FeTi_joint_edit, 'String', num2str(FeTi))


function [LegSelection, CoFe, FeTi, TiTa] = getTrackingInfo(handles)
% get the appropriate frame labels/information
[rep, cond, flynum] = getflyinfo(handles);
LegList = {'Left_Front', 'Left_Middle', 'Left_Rear',...
           'Right_Front', 'Right_Middle', 'Right_Rear'};
selectNum = (get(handles.TrackingLegSelectionDroplist, 'Value'));       
LegSelection = LegList{selectNum};  
frameNum = str2double(get(handles.TrackingCurrentFrame, 'String'));
% offset = str2double(get(handles.StartTimeEdit, 'String'))*300;
try
    data = handles.tracking(flynum);
    CoFe = data.(LegSelection)(cond,rep).frame(frameNum).CoFe;
    FeTi = data.(LegSelection)(cond,rep).frame(frameNum).FeTi;
    TiTa = data.(LegSelection)(cond,rep).frame(frameNum).TiTa;
catch
    CoFe = [];
    FeTi = [];
    TiTa = [];
end


% --- Executes during object creation, after setting all properties.
function TrackingCurrentFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrackingCurrentFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% minVal = get(handles.StartTimeEdit, 'String');
% set()


function TrackingTotalFrames_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingTotalFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrackingTotalFrames as text
%        str2double(get(hObject,'String')) returns contents of TrackingTotalFrames as a double


% --- Executes during object creation, after setting all properties.
function TrackingTotalFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrackingTotalFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TrackingFrameUpButton.
function TrackingFrameUpButton_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingFrameUpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentFrame = str2double(get(handles.TrackingCurrentFrame, 'String'));
frameStep = str2double(get(handles.TrackingFrameStepInput, 'String'));
minFrame = str2double(get(handles.StartTimeEdit, 'String'))*300;
peakFrame = str2double(get(handles.StopTimeEdit, 'String'))*300;
% maxFrame = peakFrame-minFrame;

newFrame = currentFrame+frameStep;
% add in barriers here... i.e. nothing before 1 or greater than 600
if newFrame > peakFrame
    newFrame = peakFrame;
    warndlg('Last frame')
elseif newFrame < minFrame
    newFrame = peakFrame;
end
set(handles.TrackingCurrentFrame, 'String', num2str(newFrame))
set(handles.TrackingSlider, 'Value', newFrame)
% TrackingSlider_Callback(hObject, eventdata, handles)
try close(handles.fig); catch; end
Display_JA_image_button_Callback(hObject, eventdata, handles)


% --- Executes on button press in TrackingFrameDownButton.
function TrackingFrameDownButton_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingFrameDownButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentFrame = str2double(get(handles.TrackingCurrentFrame, 'String'));
frameStep = str2double(get(handles.TrackingFrameStepInput, 'String'));
minFrame = str2double(get(handles.StartTimeEdit, 'String'))*300;
peakFrame = str2double(get(handles.StopTimeEdit, 'String'))*300;
maxFrame = peakFrame-minFrame;

newFrame = currentFrame-frameStep;
% add in barriers here... i.e. nothing before 1 or greater than 600
if newFrame > peakFrame
    newFrame = minFrame;
elseif newFrame < minFrame
    newFrame = peakFrame;
end
set(handles.TrackingCurrentFrame, 'String', num2str(newFrame))
set(handles.TrackingSlider, 'Value', newFrame)
% TrackingSlider_Callback(hObject, eventdata, handles)
try close(handles.fig); catch; end
Display_JA_image_button_Callback(hObject, eventdata, handles)


% --- Executes on slider movement.
function TrackingSlider_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% set(hObject, 'SliderStep', [1/599, 1])
currentloc = round(get(hObject, 'Value'));
set(hObject, 'Value', currentloc)
newFrame = num2str(currentloc);
set(handles.TrackingCurrentFrame, 'String', newFrame)
% UPDATE THE JOINT INFORMATION
[~, CoFe, FeTi] = getTrackingInfo(handles);
set(handles.FeTi_joint_edit, 'String', num2str(FeTi))
set(handles.CoFe_joint_edit, 'String', num2str(CoFe))


% --- Executes during object creation, after setting all properties.
function TrackingSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrackingSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
minVal = 0.3*300;
maxVal = 1.42*300;
totalRange = maxVal-minVal;
set(hObject, 'Min', minVal, 'Max', maxVal, 'Value', minVal)%, 'SliderStep', [1,1]
set(hObject, 'SliderStep', [1/totalRange, 10/totalRange])



% --- Executes on button press in TrackingBuildDataStructure.
function TrackingBuildDataStructure_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingBuildDataStructure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if there aren't structures for the longitudinal joint angles, add it  

% list of possible fields
LegSelection = getTrackingInfo(handles);
% LegList = {'Left_Front', 'Left_Middle', 'Left_Rear',...
%            'Right_Front', 'Right_Middle', 'Right_Rear'};
% selectNum = (get(handles.TrackingLegSelectionDroplist, 'Value'));       
% LegSelection = LegList{selectNum};    
fprintf('\n Field Choice:   '); disp(LegSelection)

% Build the tracking data structure:
starttime = str2double(get(handles.StartTimeEdit, 'String'));
endtime = str2double(get(handles.StopTimeEdit, 'String'));
totalFrames = 2*300;

try tracking = handles.tracking;
catch
    fprintf('\n Building tracking structure \n')
    handles.tracking = struct;
end
 
A = isfield(handles.tracking, LegSelection); %TF leg data exists within tracking

if A == 0 % leg doesn't exist, but does tracking?
    B = exist('tracking'); 
    if B == 1 %tracking exists just not that leg 
        tracking = handles.tracking; % need to load tracking to add to it
    end
    % add the leg grouping to the tracking structure
    % Build the Tracking Structure
    for kk = 1:handles.num.flies
      for icond = 1:handles.num.conds
        for irep = 1:handles.num.reps
           for iFrames = 1:totalFrames
              tracking(kk).(LegSelection)(icond, irep).frame(iFrames) = struct('pos',[], 'CoFe',[], 'FeTi',[], 'TiTa', []);
              tracking(kk).(LegSelection)(icond, irep).FrameRange = [starttime, endtime];
           end
        end
      end
    end
    fprintf(['\n Built structure for ' LegSelection '\n'])
else %leg tracking structure already exists
    fprintf(['\n Structure exists for ' LegSelection '\n'])
end
setTotalFrames(handles)
handles.tracking = tracking;
guidata(hObject, handles);


function handles = TrackFrame(handles, FrameNumber)
[rep, cond, flynum] = getflyinfo(handles);
fly_date = handles.excelfile{flynum, 1}; 
fly_num = handles.excelfile{flynum, 2};
fly_ID = get(handles.FlyIDEdit, 'String');  
% camlist = get(handles.cameralistbox, 'String');
% Cam = camlist{get(handles.cameralistbox, 'Value')};
Cam = 'B';
% fpslist = get(handles.fpslistbox, 'String');
% fps = str2double(fpslist{get(handles.fpslistbox, 'Value')});
starttime = FrameNumber;
% stoptime = starttime;

% Display Frame:
drive_root = get(handles.filerootedit, 'String');
vid_root = [drive_root char(fly_date) '\Fly ' char(fly_num) '\Raw Video\'];

vid_name = [vid_root, fly_ID, ' R' num2str(rep), 'C', num2str(cond), ' Cam-' Cam '*'];
FilePath = dir(vid_name);
if isempty(FilePath)
    warndlg('Check video file path')
    return
end
vid_name = FilePath.name;
% disp([vid_root, vid_name]);
handles.fig = image_still([vid_root, vid_name], starttime, 1.8);

% zoom status
fig = handles.fig;
switch get(handles.TrackingZoomButton,'Value')
    case 0
        X= [0.5000 836.5000];
        Y= [0.5000 639.5000];
    case 1
        X = [159.9539 586.0724];
        Y = [25.1079 350.8132];
end
set(fig.CurrentAxes, 'Xlim', X, 'YLim', Y)


hold all;
for ii = 1:5
    [x(ii), y(ii)] = crosshairs(1);
    scatter(x(ii),y(ii),'filled', 'y')
end

% coxa femur:
a = [x(1), y(1)];
b = [x(2), y(2)];
c = [x(3), y(3)];
u(1) = (b(1)-a(1));
u(2) = (b(2)-a(2));
u(3) = 0;
v(1) = (c(1)-b(1));
v(2) = (c(2)-b(2));
v(3) = 0;
CoFe = 180-atan2d(norm(cross(u,v)),dot(u,v));
set(handles.CoFe_joint_edit, 'string', CoFe)
% disp(ThetaInDegrees)

% femur tibia:
a = [x(2), y(2)];
b = [x(3), y(3)];
c = [x(4), y(4)];
u(1) = (b(1)-a(1));
u(2) = (b(2)-a(2));
u(3) = 0;
v(1) = (c(1)-b(1));
v(2) = (c(2)-b(2));
v(3) = 0;
FeTi = 180-atan2d(norm(cross(u,v)),dot(u,v));
set(handles.FeTi_joint_edit, 'string', FeTi)
% disp(ThetaInDegrees)

% tibia tarsus:
a = [x(3), y(3)];
b = [x(4), y(4)];
c = [x(5), y(5)];
u(1) = (b(1)-a(1));
u(2) = (b(2)-a(2));
u(3) = 0;
v(1) = (c(1)-b(1));
v(2) = (c(2)-b(2));
v(3) = 0;
TiTa = 180-atan2d(norm(cross(u,v)),dot(u,v));

% UPDATE DATA structures:
singleStatus = get(handles.JA_radiobutton, 'Value');
seriesStatus = get(handles.TrackingButton, 'Value');

if singleStatus == 1
    handles.JA_Data(flynum).angles(cond,rep).data = [x',y'];
    handles.JA_Data(flynum).CoFe{cond,rep} = CoFe;
    handles.JA_Data(flynum).FeTi{cond,rep} = FeTi;
elseif seriesStatus == 1
    LegList = {'Left_Front', 'Left_Middle', 'Left_Rear',...
           'Right_Front', 'Right_Middle', 'Right_Rear'};
    selectNum = (get(handles.TrackingLegSelectionDroplist, 'Value'));       
    LegSelection = LegList{selectNum};    
    StartFrame = str2double(get(handles.TrackingCurrentFrame, 'String'));
%     StartBuffer = (str2double(get(handles.StartTimeEdit, 'String')))*300;
    Frame = StartFrame;  %+StartBuffer;
    handles.tracking(flynum).(LegSelection)(cond, rep).frame(Frame).pos = [x',y'];
    handles.tracking(flynum).(LegSelection)(cond, rep).frame(Frame).CoFe = CoFe;
    handles.tracking(flynum).(LegSelection)(cond, rep).frame(Frame).FeTi = FeTi;
    handles.tracking(flynum).(LegSelection)(cond, rep).frame(Frame).TiTa = TiTa;
    
    try
        plotNewDataPoints(handles, CoFe, FeTi, TiTa)
    catch
    end
end

figure(handles.figure1);
% handles.filerootedit;


function setTotalFrames(handles)
% update the total frame number edit in the joint angles section
frameStep = str2double(get(handles.TrackingFrameStepInput, 'String'));
minFrame = str2double(get(handles.StartTimeEdit, 'String'))*300;
peakFrame = str2double(get(handles.StopTimeEdit, 'String'))*300;
numFrames = round((peakFrame-minFrame)/frameStep);
set(handles.TrackingTotalFrames, 'String', num2str(numFrames))

set(handles.TrackingSlider, 'Min', minFrame, 'Max', peakFrame)
set(handles.TrackingSlider, 'SliderStep', [1/numFrames, 10/numFrames])
currentFrame = str2double(get(handles.TrackingCurrentFrame, 'String'));
if currentFrame > peakFrame
    currentFrame = peakFrame;
    set(handles.TrackingCurrentFrame, 'String', num2str(currentFrame))
elseif currentFrame < minFrame
    currentFrame = minFrame;
    set(handles.TrackingCurrentFrame, 'String', num2str(currentFrame))
end
set(handles.TrackingSlider, 'Value', currentFrame)


% --- Executes on button press in TrackingPlotData.
function TrackingPlotData_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingPlotData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rep, cond, flynum] = getflyinfo(handles);
LegSelection = getTrackingInfo(handles);
data = handles.tracking(flynum).(LegSelection)(cond,rep);
% extract the joint angle information
for ii = 1:600
    a = data.frame(ii).CoFe;
    b = data.frame(ii).FeTi;
    c = data.frame(ii).TiTa;
    if ~isempty(a)
        pltdata(ii,1) = a;
        pltdata(ii,2) = b;
        pltdata(ii,3) = c;
    else 
        pltdata(ii,1) = nan;
        pltdata(ii,2) = nan;
        pltdata(ii,3) = nan;
    end
end
sz = 40;
x = 1:600;
fig = getfig('', 1);  
% subplot(3,1,1)
hold all
xlim([0 600])
ylim([0 200])
scatter(x,pltdata(:,1), sz, 'blue', 'filled')
scatter(x,pltdata(:,2), sz, 'red', 'filled')
scatter(x,pltdata(:,3), sz, 'black', 'filled')
vline(150, 'g--')
vline(366, 'g--')
title({handles.structure_name;...
    ['Fly: ' num2str(flynum) ' Cond: ' num2str(cond) ' Rep: ' num2str(rep)]})
xlabel('Time (video frame)')
ylabel('Joint angle (degree)')
% legend('CoFe','FeTi', 'TiTa')
startLabeling = str2double(get(handles.StartTimeEdit, 'String'))*300;
stopLabeling = str2double(get(handles.StopTimeEdit, 'String'))*300;
vline([startLabeling, stopLabeling], 'k:')
handles.angleFig = fig;
guidata(hObject, handles);

function plotNewDataPoints(handles, CoFe, FeTi, TiTa)
% plot new data points as they are calculated from tracking manually

figure(handles.angleFig) %make fig active
% [~, CoFe, FeTi, TiTa] = getTrackingInfo(handles);%get state info
hold all
sz = 40;
x = str2double(get(handles.TrackingCurrentFrame, 'String'));
scatter(x,CoFe, sz, 'blue', 'filled')
scatter(x,FeTi, sz, 'red', 'filled')
scatter(x,TiTa, sz, 'black', 'filled')
figure(handles.figure1)


% --- Executes on button press in TrackingZoomButton.
function TrackingZoomButton_Callback(hObject, eventdata, handles)
% hObject    handle to TrackingZoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TrackingZoomButton
% zoom status
if isvalid(handles.fig)
    fig = handles.fig;
    switch get(hObject,'Value')
    case 0
        disp('Zoom OFF')
        X= [0.5000 836.5000];
        Y= [0.5000 639.5000];
    case 1
        disp('Zoom ON')
        X = [159.9539 586.0724];
        Y = [25.1079 350.8132];
    end
    set(fig.CurrentAxes, 'Xlim', X, 'YLim', Y)
else
    switch get(hObject,'Value')
    case 0
        disp('Zoom OFF')
    case 1
        disp('Zoom ON')
    end
end
