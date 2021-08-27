
function [path, folder] = getCloudPath(folderOption)
% path = getCloudPath(DateSelOpt)
% folderOption --> list and select a folder from the basepath
%                  1) basepath & folder selection merged
%                  2) folder is just the folder name
% folder --> folder/path selected from path
% 
% 
% ES Dickinson, Yale University, Aug 2021

    switch getenv('COMPUTERNAME')
        case 'DENALI'
            path = 'G:\My Drive\Jeanne Lab\DATA\';
        case 'TOGIAK'
            path = 'G:\My Drive\Jeanne Lab\DATA\';
        case 'EVYNPC'
            path = 'G:\My Drive\Jeanne Lab\DATA\';
    end
    
    if nargin == 1
        dirc = dir(path);
        dirc = flip(dirc(find(~cellfun(@isdir,{dirc(:).name}))));
        folderNames = ['Today', {dirc(:).name}];
        indx = listdlg('ListString', folderNames, 'SelectionMode', 'Single');
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


end