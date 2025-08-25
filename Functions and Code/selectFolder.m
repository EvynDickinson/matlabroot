

function folderSelected = selectFolder(path,mulitselect,promptString,default_search_names)
% folderSelected = selectFolder(path,mulitselect,promptString,default_search_names)
% default_search_names = {'name 1', 'name 2'} search for these names and
% highlight them as the default start values...
% select a folder from an input directory with the auto added 

% path_A = 'C:\Users\evynd\OneDrive\Desktop\Testing End Folder\'

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

