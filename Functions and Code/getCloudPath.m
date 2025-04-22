
function [path, folder] = getCloudPath(folderOption)
% path = getCloudPath(DateSelOpt)
% folderOption --> list and select a folder from the basepath
%                  1) basepath & folder selection merged
%                  2) folder is just the folder name
%                  3) folder is the path + 'structures' folder
%                  4) outvar_1 = path to data folder, outvar2 = date folder, outvar3 = exp name
% if no folderOption is not included, just the environment is given
% folder --> folder/path selected from path
% 
% 
% ES Dickinson, Yale University, Aug 2021

switch getenv('COMPUTERNAME')
    case 'DENALI'
        path = 'G:\My Drive\Jeanne Lab\DATA\';
    case 'ACADIA'
        path = 'G:\My Drive\Jeanne Lab\DATA\';
    case 'TOGIAK'
        path = 'G:\My Drive\Jeanne Lab\DATA\';
    case 'EVYNPC'
        path = 'G:\My Drive\Jeanne Lab\DATA\';
    case 'MWMJ0LY4WH'
        path = 'G:\My Drive\Jeanne Lab\DATA\';
    case 'SLEEPINGGIANT'
        path = 'G:\My Drive\Jeanne Lab\DATA\';
%     case '' %shows up as empty on the mac
% %         disp('Evyn''s Macbook');
%         path = '/Volumes/GoogleDrive/My Drive/Jeanne Lab/DATA/';
    case ''
        path = '/Users/evyndickinson/Library/CloudStorage/GoogleDrive-evyn.dickinson@yale.edu/My Drive/Jeanne Lab/DATA/';
end 

if nargin == 1
    if folderOption==3 %no folder, just data path
        folder = [path 'Data structures/'];
    else
        dirc = dir(path);
        dirc = flip(dirc);
%         dirc = flip(dirc(cellfun(@isfolder,{dirc(:).name})));
%         dirc = flip(dirc(find(~cellfun(@isdir,{dirc(:).name}))));
        folderNames = ['Today', {dirc(:).name}];
        indx = listdlg('ListString', folderNames, 'SelectionMode', 'Single');
        if isempty(indx)
            fprintf('\n No folder selected\n')
            folder = '';
            return
        end
        if strcmpi(folderNames{indx}, 'Today')==true
            dir_sel = strrep(datestr(datetime,'mm-dd-yyyy'),'-','.');
        else
            dir_sel = folderNames{indx};
        end
        switch folderOption
            case 1
                folder = fullfile(path, dir_sel);
            case 2
                folder = dir_sel;

        end  
    end
elseif nargin == 0
    folder = [];
end


end