



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
else
    computerName = getenv('COMPUTERNAME');
end
paths.computerName = computerName;
% note: the computer name is a host VPN on a mac when there is a VPN connection

% BASE DRIVE LOCATION SELECTION

% portable drive test
[serverDrive, portableDrive] = deal(false);
matchingDrives = findDriveByName(paths.removableDrive); % test if there is a removable drive attached
if ~isempty(matchingDrives)
    portableDrive = true;
    onthegoPath = [matchingDrives{:} '/'];
else
    onthegoPath = [];
end

% permanent drive test
permanentDrive = any(strcmp(computerName, paths.fixedDriveLocations)); % fixed drive yes or no
permanentPath = [];

% server drive test
serverPath = [];
switch computerName
    case 'ACADIA'
        serverPath = paths.acadiaServerPath;
        permanentPath = paths.acadiaLocalPath;
    case 'TOGIAK'
        serverPath = paths.togiakServerPath;
        permanentPath =[];
    case 'EVYNPC'
        serverPath = paths.EvynPCServerPath;
        permanentPath = paths.EvynPCLocalPath;
    case 'vpn1722513182.its.yale.internal' % VPN into Yale on Mac
        serverPath = paths.EvynMacServerPath;
        permanentPath = [];
    case 'Evyns-M3-MacBook-Pro.local' % Mac, no VPN thus no server
        serverPath = [];
        permanentPath = [];
end 
if exist(serverPath, 'dir') == 7
    serverDrive = true;
end

% Output folder paths for available options

availablePaths.serverPath = serverPath;
availablePaths.permanentPath = permanentPath;
availablePaths.onthegoPath = onthegoPath;

availableDrives.server = serverDrive;
availableDrives.permanent = permanentDrive;
availableDrives.portable = portableDrive;

% Return error if there aren't available drives
if ~permanentDrive && ~serverDrive && ~portableDrive
    warndlg({'No local, portable, or server directories available';...
             'Try reestablishing a VPN connection and/or checking the OnTheGoDrive'})
end











