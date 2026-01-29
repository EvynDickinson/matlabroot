

function [path, location_name] = getDataPath(dataType, dataLocation, promptString)
% path = getDataPath(dataType, dataLocation, promptString)
% 
% dataType options:
% 0) user select (default)
% 1) single trial data
% 2) raw data folder
% 3) pooled trial data structures 
% 4) grouped data structures
% 5) base folder only! -- no path ending
% 6) courtship folder
% 7) NOAA Data folder
% 
%
% dataLocation options: 
% 0) user select (default)
% 1) local path
% 2) server path
% 3) removable drive
% 4) storage drive
% 5) server path 2
% 6) either server (user selected each time)
% 
% ES Dickinson, Yale University, June 2024

% %%%%%%%%%%%%%%%%%%%           EDITABLE PATH LOCATIONS        %%%%%%%%%%%%%%%%%%%
% removableDrive = 'OnTheGoData';
% single_trial = 'Trial Data/';
% raw_data = 'DATA/';
% fixedDriveLocations  = {'ACADIA','EVYNPC'}; % computer names that have permanent local drives
% acadiaServerPath = 'S:\Evyn\';
% togiakServerPath = 'S:\Evyn\';
% EvynPCServerPath = '\\svalbard.med.yale.internal\shared\Evyn\';
% EvynMacServerPath = '/Volumes/shared/Evyn/';
% acadiaLocalPath = 'D:\';
% EvynPCLocalPath = 'K:\DATA\';
% path = []; % empty so that in the end, there is something assigned
% %%%%%%%%%%%%%%%%%%        AVOID EDITING CODE BELOW HERE        %%%%%%%%%%%%%%%%%

% +++ TO UPDATE PATH NAMES EDIT THESE TWO +++
% getDataFolderOptions **edit list near the bottom**
% getPathNames **edit around line 18 where the names of the paths are**

paths = getPathNames; 
path = [];


% default values -- single trial data, user selected drive
switch nargin
    case 0
        dataType = 0;
        dataLocation = 0;
        promptString = 'Select desired type of data folder : ';
    case 1 
        dataLocation = 0;
        promptString = 'Select desired type of data folder : ';
    case 2
        promptString = 'Select desired type of data folder : ';
end


% Pull the type of data path name (e.g. final folder)
switch dataType
    case 0 % user select
        optionList = {'Raw Data', 'Single Trial', 'Grouped Trials', 'Grouped Structures','Courtship', 'NOAA'};
        idx = listdlg('PromptString', 'Select type of data desired:', 'ListString', optionList, 'ListSize',[160,200]);
        if isempty(idx)
            disp('No data selected')
            return
        end
         switch optionList{idx}
             case 'Single Trial'
                 path_end = paths.single_trial;
             case 'Raw Data'
                 path_end = paths.raw_data;
             case 'Grouped Trials'
                 path_end = paths.grouped_trials;
             case 'Grouped Structures'
                 path_end = paths.group_comparision;
             case 'Courtship'
                 path_end = paths.courtship;
             case 'NOAA'
                 path_end = paths.NOAA;
         end
    case 1 % single trial data
        path_end = paths.single_trial;
    case 2 % raw data
        path_end = paths.raw_data;
    case 3 % pool trial data
        path_end = paths.grouped_trials;
    case 4 % group structure comparison
        path_end = paths.group_comparision;
    case 5 % no path ending -- user just wants the base folder
        path_end = '';
    case 6 % courtship
        path_end = paths.courtship;
    case 7 % NOAA
        path_end = paths.NOAA;
end

% BASE DRIVE LOCATION SELECTION
% find the list of available drives
[availableDrives, availablePaths] = getDataFolderOptions;

% are they available? 
portableDrive = availableDrives.portable;
serverDrive = availableDrives.server;
serverTwoDrive = availableDrives.serverTwo;
permanentDrive = availableDrives.permanent;
storageDrive = availableDrives.storage;

% Return error if there aren't available drives
if ~permanentDrive && ~serverDrive && ~serverTwoDrive && ~portableDrive && ~storageDrive
    warndlg({'No local, portable, or server directories available';...
             'Try reestablishing a VPN connection and/or checking the OnTheGoDrive'})
    path = [];
    return
end

% if the hard-coded drive isn't available, switch to manual selection: 
if dataLocation==1 && ~permanentDrive
    dataLocation = 0;
elseif dataLocation==2 && ~serverDrive
     dataLocation = 0;
elseif dataLocation==3 && ~portableDrive
    dataLocation = 0;
elseif dataLocation==4 && ~storageDrive
    dataLocation = 0;
elseif dataLocation==5 && ~serverTwoDrive
    dataLocation = 0;
elseif dataLocation==6 && ~serverTwoDrive && ~serverDrive
    dataLocation = 0;
end

% option to select the specific server if both are present and 'just
% server' was selected (aka option 6)
if dataLocation==6 && serverTwoDrive && serverDrive
    switch questdlg('Which server do you want?','Server Selection', 'Svalbard', 'Server Two','Cancel', 'Svalbard')
        case 'Svalbard'
            dataLocation = 2;
        case 'Server Two'
            dataLocation = 5;
        case {'Cancel', ''}
            disp('No drive selected')
            dataLocation = 0;
    end
elseif dataLocation==6 && ~serverTwoDrive && serverDrive
    dataLocation = 2;
elseif dataLocation==6 && serverTwoDrive && ~serverDrive
    dataLocation = 5;
end

% Appropriate drive location selection
if dataLocation==0 % USER INPUT SELECTION
    all_locations = {'Svalbard', 'Server Two', 'Local', 'Storage', 'On The Go', 'Manually Select'};
    loc = [availableDrives.server, availableDrives.serverTwo, availableDrives.permanent,...
               availableDrives.storage, availableDrives.portable,true];
    loc_options = all_locations(loc);
    idx = listdlg("PromptString",promptString,'ListString',loc_options,'SelectionMode','single',...
        'OKString','Get It', 'CancelString','Heck No','InitialValue',1,'ListSize',[300,150]);
    if isempty(idx)
        disp('Bummer, no data drive selected')
        return
    end
    
    switch loc_options{idx}
        case 'Svalbard'
            drive_loc_sel = 2;
        case 'Server Two'
            drive_loc_sel = 5;
        case 'Local'
            drive_loc_sel = 1;
        case 'Storage'
            drive_loc_sel = 4;
        case 'On The Go'
            drive_loc_sel = 3;
        case 'Manually Select'
            drive_loc_sel = 0;
    end

    dataLocation = drive_loc_sel;
end

% CONSTRUCT THE DATA PATH
switch dataLocation
    case 0 % user select path manually:
        basePath = uigetdir;
        path = [basePath, '/' path_end];   
        location_name  = 'unknown';
    case 1 % local permanent path
        path = [availablePaths.permanentPath, path_end];
        location_name  = 'internal drive';
    case 2 % server path
        path = [availablePaths.serverPath, path_end];
        location_name  = 'svalbard';
    case 3 % portable drive
        path = [availablePaths.onthegoPath, path_end];
        location_name  = 'OnTheGo';
    case 4 % storage drive
        path = [availablePaths.storagePath, path_end];
        location_name  = 'data storage';
    case 5 % portable drive
        path = [availablePaths.serverTwoPath, path_end];
        location_name  = 'server two';
end

% disp(path)
if ~(exist('path','var')==1)
    disp('Warning: selected path not found')
    path = [];
end






