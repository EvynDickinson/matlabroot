

function [folderSelected, fullPath] = selectFolder(path,mulitselect,promptString,default_search_names)
% [folderSelected, fullPath] = selectFolder(path,mulitselect,promptString,default_search_names)
% 
% PURPOSE
% select folder name(s) from the given directory or folder
% 
% INPUTS 
%   'path' :  filepath to the current folder that has the desired folder paths 
%   'multiselect' : (optional) allow for selection of mulitple folders 
%           'single' -- only can select one folder
%           'multiple' -- can select multiple folder options
%   'promptString' : (optional) text displayed in the folder selection box
%   'default_search_names'  : (optional) value to start selected if present 
%           e.g. if the folder names are {'A', 'B', 'C'} and you want the
%           inital selection option highlighted to be 'B' set default_search_names to 'B'
% 
% OUTPUTS
%   'folderSelected' : cell with string name of the selected folder or cell with
%           multiple folders names
%   'fullPath' : full path string to the input folder
%           fullPath can only be returned if multiselect if set to
%           'single' otherwise it returns the input root 'path'
%
% ES DICKINSON, 2024

%%
switch nargin 
    case 1
        promptString = 'Select folders';
        mulitselect = 'Single'; % default (other option: multiple)
    case 2
        promptString = 'Select folders';
end

%TODO: screen for only folders? 

dirc = dir(path);

% Remove files from the list of directory contents
folders = dirc([dirc.isdir]);

% Get the names of the folders
folderNames = {folders.name}';
folderNames(1:2) = [];
if isempty(folderNames)
    fprintf('\n No available folders \n')
    folderSelected = '';
    return
end

% if start search values, look for those: 
if nargin==4 && ~isempty(default_search_names)
    start_loc = [];
    for i = 1:length(default_search_names)
        start_loc(i) = find(strcmpi(default_search_names{i},folderNames)); %#ok<AGROW>
    end
else 
    start_loc = 1;
end

indx = listdlg('ListString', folderNames, 'SelectionMode', mulitselect,'PromptString', promptString, 'ListSize',[450 600],'InitialValue',start_loc);
nFolders = length(indx);

if isempty(indx)
    fprintf('\n No folder selected\n')
    folderSelected = '';
    return
end

dir_sel = cell(1,nFolders);
for i = 1:nFolders
    dir_sel{i} = folderNames{indx(i)};
end

folderSelected = dir_sel;

if strcmpi(mulitselect, 'single')
    if strcmp(path(end),'/') ||  strcmp(path(end),'\')
        fullPath = [path folderSelected{1} '/'];
    else
        fullPath = [path '/' folderSelected{1} '/'];
    end
else 
    fullPath = path;
end









