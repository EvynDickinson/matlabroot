

function path = getDataPath(rawORsingle, localORserver)
% path = getDataPath(rawORsingle, localORserver)
% 
% rawORsingle:
% 0) user select (default)
% 1) single trial data
% 2) raw data folder
%
% localORserver: 
% 0) user select (default)
% 1) local path
% 2) server path
% 3) removable drive
% 
% ES Dickinson, Yale University, June 2024

%%%%%%%%%%%%%%%%%%%           EDITABLE PATH LOCATIONS         %%%%%%%%%%%%%%%%%%%
removableDrive = 'OnTheGoData';
single_trial = 'Trial Data/';
raw_data = 'DATA/';
fixedDriveLocations  = {'ACADIA','EVYNPC'}; % computer names that have permanent local drives
acadiaServerPath = 'S:\Evyn\Trial Data\';
togiakServerPath = 'S:\Evyn\Trial Data\';
EvynPCServerPath = '\\svalbard.med.yale.internal\shared\Evyn\Trial Data\';
EvynMacServerPath = 'TBD';
acadiaLocalPath = 'D:\';
EvynPCLocalPath = 'K:\';
%%%%%%%%%%%%%%%%%%        AVOID EDITING CODE BELOW HERE        %%%%%%%%%%%%%%%%%


% Get computer name
if ismac
    [~, result] = system('hostname');
    computerName = strtrim(result);
else
    computerName = getenv('COMPUTERNAME');
end

% default values -- single trial data, user selected drive
switch nargin
    case 0
        rawORsingle = 0;
        localORserver = 0;
    case 1 
        localORserver = 0;
end


% RAW OR SINGLE DATA PATH ENDING SELECTION
switch rawORsingle
    case 0 % user select
         switch questdlg('Select desired type of data folder : ', '', 'Single Trial', 'Raw Data','Cancel','Single Trial')
             case {'Cancel',''}
                 return
             case 'Single Trial'
                 path_end = single_trial;
             case 'Raw Data'
                 path_end = raw_data;
         end
    case 1 % single trial data
        path_end = single_trial;
    case 2 % raw data
        path_end = raw_data;
end

% BASE DRIVE LOCATION SELECTION

% portable drive test
[serverDrive, portableDrive] = deal(false);
matchingDrives = findDriveByName(removableDrive); % test if there is a removable drive attached
if ~isempty(matchingDrives)
    portableDrive = true;
end
% permanent drive test
permanentDrive = any(strcmp(computerName, fixedDriveLocations)); % fixed drive yes or no
% server drive test
switch computerName
    case 'ACADIA'
        serverPath = acadiaServerPath;
        permanentPath = acadiaLocalPath;
    case 'TOGIAK'
        serverPath = togiakServerPath;
        permanentPath =[];
    case 'EVYNPC'
        serverPath = EvynPCServerPath;
        permanentPath = EvynPCLocalPath;
    case '' % TODO -- update mac name
        serverPath = EvynMacServerPath;
        permanentPath = [];
end 
if exist(serverPath, 'dir') == 7
    serverDrive = true;
end

% if the coded drive isn't available, switch to manual selection: 
if localORserver==1 && ~permanentDrive
    localORserver = 0;
elseif localORserver==2 && ~serverDrive
     localORserver = 0;
elseif localORserver==3 && ~portableDrive
    localORserver = 0;
end

% Appropriate drive location selection
if localORserver==0 % USER INPUT SELECTION
        if permanentDrive && portableDrive && serverDrive % choice of all three locations
            switch questdlg('Select data drive location:', 'Drive Location', 'Local', 'Server', 'On The Go','Local')
                case 'Local' 
                        drive_loc_sel = 1; 
                case 'Server'
                        drive_loc_sel = 2;
                case 'On The Go'
                        drive_loc_sel = 3;
                case ''
                        return
            end
        elseif permanentDrive && portableDrive && ~serverDrive % local options only
            switch questdlg('Select data drive location:', 'Drive Location', 'Local', 'On The Go', 'Cancel','Local')
                case 'Local' 
                        drive_loc_sel = 1;
                case 'On The Go'
                        drive_loc_sel = 3;
                 case {'Cancel',''}
                     return
            end
        elseif permanentDrive && ~portableDrive && serverDrive % no portable drive
            switch questdlg('Select data drive location:', 'Drive Location', 'Local', 'Server', 'Cancel','Local')
                case 'Local' 
                        drive_loc_sel = 1;
                case 'Server'
                        drive_loc_sel = 2;
                 case {'Cancel',''}
                     return
            end
        elseif ~permanentDrive && portableDrive && serverDrive % server and portable
            switch questdlg('Select data drive location:', 'Drive Location', 'On The Go', 'Server', 'Cancel', 'On The Go')
                case 'On The Go' 
                        drive_loc_sel = 3;
                case 'Server'
                        drive_loc_sel = 2;
                 case {'Cancel',''}
                     return
            end
         elseif ~permanentDrive && ~portableDrive && serverDrive % server only
                drive_loc_sel = 2;
         elseif permanentDrive && ~portableDrive && ~serverDrive % permanent local drive only
                drive_loc_sel = 1;
         elseif ~permanentDrive && portableDrive && ~serverDrive % portable drive only
                drive_loc_sel = 3;
        end
        localORserver = drive_loc_sel;
end


% CONSTRUCT THE DATA PATH
switch localORserver
    case 1 % local permanent path
        path = [permanentPath, path_end];
    case 2 % server path
        path = [serverPath, path_end];
    case 3 % portable drive
        path = [matchingDrives{:}, '/', path_end];
end

% disp(path)
if ~exist(path, 'var')
    disp('Warning: selected path not found')
    path = [];
end









