
function [availableDrives, availablePaths] = getDataFolderOptions
% path = getDataPath(rawORsingle, localORserver, promptString)
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
%
% %%%%%%%%%%%%%%%%%%%           EDITABLE PATH LOCATIONS         %%%%%%%%%%%%%%%%%%%
% paths.removableDrive = 'OnTheGoData';
% paths.single_trial = 'Trial Data/';
% paths.raw_data = 'DATA/';
% paths.fixedDriveLocations  = {'ACADIA','EVYNPC'}; % computer names that have permanent local drives
% paths.acadiaServerPath = 'S:\Evyn\';
% paths.togiakServerPath = 'S:\Evyn\';
% paths.EvynPCServerPath = '\\svalbard.med.yale.internal\shared\Evyn\';
% paths.EvynMacServerPath = '/Volumes/shared/Evyn/';
% paths.acadiaLocalPath = 'D:\';
% paths.EvynPCLocalPath = 'K:\';
% %%%%%%%%%%%%%%%%%%        AVOID EDITING CODE BELOW HERE        %%%%%%%%%%%%%%%%%
 
paths = getPathNames;

% Get computer name
if ismac
    [~, result] = system('hostname');
    computerName = strtrim(result);
    if contains(computerName,'yale')
        computerName = 'Yale VPN';
    end
else
    computerName = getenv('COMPUTERNAME');
end
paths.computerName = computerName;
% note: the computer name is a host VPN on a mac when there is a VPN connection

% BASE DRIVE LOCATION SELECTION

% portable drive test
[serverDrive, portableDrive, storageDrive, serverTwoDrive] = deal(false);
matchingDrives = findDriveByName(paths.removableDrive); % test if there is a removable drive attached
if ~isempty(matchingDrives)
    portableDrive = true;
    onthegoPath = [matchingDrives{:} '/'];
else
    onthegoPath = [];
end

% storage drive:
matchingDrives = findDriveByName(paths.storageDrive); % test if there is a removable drive attached
if ~isempty(matchingDrives)
    storageDrive = true;
    storagePath = [matchingDrives{:} '/'];
else
    storagePath = [];
end

% permanent drive test
permanentDrive = any(strcmp(computerName, paths.fixedDriveLocations)); % fixed drive yes or no
permanentPath = [];


% server drive test
[serverPath, serverTwoPath, permanentPath]  = deal([]);

switch computerName
    case 'ACADIA'
        serverPath = paths.acadiaServerPath;
        serverTwoPath = paths.acadiaServerTwoPath;
        permanentPath = paths.acadiaLocalPath;
    case 'TOGIAK'
        serverPath = paths.togiakServerPath;
        serverTwoPath = paths.togiakServerTwoPath;
    case 'EVYNPC'
        serverPath = paths.EvynPCServerPath;
        serverTwoPath = paths.EvynPCServerTwoPath;
        permanentPath = paths.EvynPCLocalPath;
    case 'Yale VPN' % VPN into Yale on Mac
        serverPath = paths.EvynMacServerPath;
        %TODO (3/1)
    case 'Evyns-M3-MacBook-Pro.local' % Mac, no VPN thus no server
        % 
    case 'DENALI'
        serverPath = paths.denaliServerPath;
    case 'MWMJ0LY4WH' %becca's computer
        serverPath = paths.chilkatServerPath;
        serverTwoPath = paths.chilkatServerTwoPath;
    case 'SLEEPINGGIANT'
        serverPath = paths.SGServerPath;
        serverTwoPath = paths.SGServerTwoPath;
end 

if isfolder(serverPath)
    serverDrive = true;
end
% TODO: sort out the time delay on the server two function call
% Reinstate this function when server2 is being actively used
% if isfolder(serverTwoPath)
    % serverTwoDrive = true;
% end
serverTwoDrive = false;

% Output folder paths for available options

availablePaths.serverPath = serverPath;
availablePaths.serverTwoPath = serverTwoPath;
availablePaths.permanentPath = permanentPath;
availablePaths.onthegoPath = onthegoPath;
availablePaths.storagePath = storagePath;

availableDrives.server = serverDrive;
availableDrives.serverTwo = serverTwoDrive;
availableDrives.permanent = permanentDrive;
availableDrives.portable = portableDrive;
availableDrives.storage = storageDrive;

% Return error if there aren't available drives
if ~permanentDrive && ~serverDrive && ~portableDrive && ~serverTwoDrive
    warndlg({'No local, portable, or server directories available';...
             'Try reestablishing a VPN connection and/or checking the OnTheGoDrive'})
end











