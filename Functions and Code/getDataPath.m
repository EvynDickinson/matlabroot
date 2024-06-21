

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
% 1) server path
% 2) local path
% 3) removable drive
% 
% ES Dickinson, Yale University, 6/2024

% path names (edit these if they change)
removableDrive = 'OnTheGoData';
single_trial = 'Trial Data/';
raw_data = 'DATA/';

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
             case 'Cancel'
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
matchingDrives = findDriveByName(removableDrive); % test if there is a removable drive attached
if ~isempty(matchingDrives)
    moveDrive = true;
end

switch localORserver



switch getenv('COMPUTERNAME')
    case 'ACADIA'
        path = 'D:\';
    case 'TOGIAK'
        warndlg('manually select data path')
        path = uigetdir;
    case 'EVYNPC'
        path = 'K:\';
    case '' %shows up as empty on the mac %TODO: find one to test here for how it shows up
        warndlg('manually select data path')
        path = uigetdir;
end 

if nargin >= 1
    switch folderOption
        case 1 % single trial folders
            basePath = [path, 'Trial Data/'];
        case 2 % raw data folder
            basePath = [path, '/DATA/'];
    end
end











