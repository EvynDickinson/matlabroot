
function folderPath = createFolder(folderPath)
% folderPath = createFolder(folderPath)
% check if a directory path exists as a folder and if not, create one.
%
%
% ES Dickinson, 2025

if ~exist(folderPath, 'dir')
    mkdir(folderPath)
end