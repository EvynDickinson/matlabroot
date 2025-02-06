

function [drive, expProtocol] = availableExperimentDrives
% drive = availableExperimentDrives
% Gives the letter drive for where the data will be saved
% And the selected temperature protocol 
% 
% ES Dickinson 2025

% check the amount of free space on the two drives to select where the data will be
% stored: 
F_free = (java.io.File("F:").getFreeSpace()*1e-12); % answer in Terabytes
E_free = (java.io.File("E:").getFreeSpace()*1e-12); % answer in Terabytes

% query what temp protocol is being used: 
a = questdlg('Select your temp protocol:','','F LRR 25-17', 'LTS 35-15','Cancel','F LRR 25-17');
% select drive & space:
switch a
    case 'F LRR 25-17'
        reqSpace = 1.5; % how many terrabytes of free space required
        expProtocol = 'courtship_F_LRR_25-17';
    case 'LTS 35-15'
        expProtocol = 'high_res_LTS_35-15';
        reqSpace = 3.8; % how many terrabytes of free space required
    case 'Cancel'
        warndlg('No drive selected')
        drive = nan;
        return
end
sizestr = ['E: ' num2str(E_free) ' (TB) | F: ' num2str(F_free) ' (TB)'];

F_available = F_free>=reqSpace;
E_available = E_free>=reqSpace;
if F_available && ~E_available
    str = 'F';
    autostr = 'F';
elseif E_available && ~F_available
    str = 'E';
    autostr = 'E';
elseif E_available && F_available
    str = 'E and F';
    autostr = 'E';
else 
    str = 'no';
    autostr = 'Cancel';
end

fullstr = ['Sufficent space available on ' str ' drive(s)'];

drive = questdlg({fullstr; ' '; sizestr;' '; 'Select drive for saving the experiment:'},'','E', 'F', 'Cancel', autostr);
if strcmp(drive, 'Cancel')
    warndlg('No drive selected')
    drive = nan;
end
