

function folderSelected = selectFolder(path,mulitselect,promptString)
% folderSelected = selectFolder(path)
% select a folder from an input directory with the autoadded 

% path_A = 'C:\Users\evynd\OneDrive\Desktop\Testing End Folder\'

switch nargin 
    case 1
        promptString = 'Select folders';
        mulitselect = 'Single'; % default (other option: multiple)
    case 2
        promptString = 'Select folders';
end

git 
dirc = dir(path);
folderNames = {dirc(:).name};
folderNames(1:2) = [];
indx = listdlg('ListString', folderNames, 'SelectionMode', mulitselect,'PromptString', promptString, 'ListSize',[300 500]);
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

