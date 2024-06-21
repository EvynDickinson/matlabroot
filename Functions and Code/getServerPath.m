

function basePath = getServerPath(folderOption)
% path = getServerPath(DateSelOpt)
% folderOption --> list and select a folder from the basepath
%                  1) basepath & folder selection merged
%                  2) folder is just the folder name
%                  3) folder is the path + 'structures' folder
% if no folderOption is not included, just the environment is given
% folder --> folder/path selected from path
% 
% 
% ES Dickinson, Yale University, Aug 2021

% eventually have it auto-detect the location, but for now:

switch getenv('COMPUTERNAME')
    % case 'DENALI'
    %     path = 'G:\My Drive\Jeanne Lab\DATA\';
    case 'ACADIA'
        path = 'S:\Evyn\Trial Data\';
    case 'TOGIAK'
        path = 'S:\Evyn\Trial Data\';
    case 'EVYNPC'
        path = '\\svalbard.med.yale.internal\shared\Evyn\Trial Data\';
    case '' %shows up as empty on the mac
%         disp('Evyn''s Macbook');
        path = '/Volumes/GoogleDrive/My Drive/Jeanne Lab/DATA/';
end 
 % working here TODO
if nargin >= 1
    switch folderOption
        case 1 % single trial folders
            basePath = [path, 'Trial Data/'];
        case 2 % raw data folder
            basePath = [path, '/DATA/'];
    end
end

if exist(path, 'dir') == 7
    serverDrive = true;
else
    warndlg('Connect to VPN and try again')
    return
end






% if nargin == 1
%     if folderOption==3 %no folder, just data path
%         folder = [path 'Data structures/'];
%     else
%         dirc = dir(path);
%         dirc = flip(dirc);
% %         dirc = flip(dirc(cellfun(@isfolder,{dirc(:).name})));
% %         dirc = flip(dirc(find(~cellfun(@isdir,{dirc(:).name}))));
%         folderNames = ['Today', {dirc(:).name}];
%         indx = listdlg('ListString', folderNames, 'SelectionMode', 'Single');
%         if isempty(indx)
%             fprintf('\n No folder selected\n')
%             folder = '';
%             return
%         end
%         if strcmpi(folderNames{indx}, 'Today')==true
%             dir_sel = strrep(datestr(datetime,'mm-dd-yyyy'),'-','.');
%         else
%             dir_sel = folderNames{indx};
%         end
%         switch folderOption
%             case 1
%                 folder = fullfile(path, dir_sel);
%             case 2
%                 folder = dir_sel;
% 
%         end  
%     end
% elseif nargin == 0
%     folder = [];
% end


end