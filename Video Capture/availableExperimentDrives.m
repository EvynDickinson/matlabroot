

function [drive, expProtocol] = availableExperimentDrives
% drive = availableExperimentDrives
% Gives the letter drive for where the data will be saved
% And the selected temperature protocol 
% 
% ES Dickinson 2025

% check the amount of free space on the two drives to select where the data will be
% stored: 
F_free = round(java.io.File("F:").getFreeSpace()*1e-12,1); % answer in Terabytes
E_free = round(java.io.File("E:").getFreeSpace()*1e-12,1); % answer in Terabytes
H_free = round(java.io.File("H:").getFreeSpace()*1e-12,1); % answer in Terabytes

matchingDrives = findDriveByName('Data Storage');
if ~isempty(matchingDrives)
    Storage = matchingDrives{:};
    Storage_free = round(java.io.File(Storage).getFreeSpace()*1e-12,1); % answer in Terabytes
else
    Storage_free = 0;
end

% query what temp protocol is being used: 
a = questdlg('Select your temp protocol:','','F LRR 25-17', 'LTS 35-15','Cancel','F LRR 25-17');
% select drive & space:
switch a
    case 'F LRR 25-17'
        reqSpace = 1.5; % how many terrabytes of free space required
        expProtocol = 'courtship_F_LRR_25-17';
    case 'LTS 35-15'
        expProtocol = 'high_res_LTS_35-15';
        reqSpace = 3.4; % how many terrabytes of free space required
    case 'Cancel'
        warndlg('No drive selected')
        drive = nan;
        return
end

F_available = F_free>=reqSpace;
E_available = E_free>=reqSpace;
H_available = H_free>=reqSpace;
storage_available = Storage_free>=reqSpace;


options = {['Storage: ' num2str(Storage_free)], ['F: ' num2str(F_free)], ['E: ' num2str(E_free)], ['H: ' num2str(H_free)]};
available = [storage_available, F_available, E_available, H_available];
list_opt = options(available);
full_List = ['Available:', list_opt,'---------- ', 'Not enough space:', options(~available)];

a = listdlg('liststring', full_List, 'promptstring', 'Select from available drives:','InitialValue',2);
drive_base = full_List{a};
b = strsplit(drive_base,':');
switch b{1}
    case 'Storage'
        if isempty(Storage)
            drive = uigetdir;
            drive = drive(1:end-1);
        else
            drive = Storage;
        end
        
    case 'E'
        drive = 'E:';
    case 'F'
        drive = 'F:';
    case 'H'
        drive = 'H:';
    case {'---------- ', 'Available', 'Not enough space'}
        disp('no valid drive selected')
        drive = nan;
        return
end

if strcmp(drive, 'Cancel')
    warndlg('No drive selected')
    drive = nan;
end






