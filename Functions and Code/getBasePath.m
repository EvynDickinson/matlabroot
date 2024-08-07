

function basePath = getBasePath(folderOption)
% path = getBasePath(DateSelOpt)
% folderOption 1 (default) : trial data folder
% folderOption 2 : raw data folder 
% 
% ES Dickinson, Yale University, Aug 2021

% TODO: set something to scan for a particular drive

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