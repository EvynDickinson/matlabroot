
function success = SyncFolderandSubfolder(startPath, targetPath, version_control)
% success = SyncFolderandSubfolder(startPath, targetPath)
% 
%
% ES DICKINSON

% % testing area: 
% clear; clc
% targetPath = 'C:\Users\evynd\Desktop\Testing Folder End\';
% startPath = 'C:\Users\evynd\Desktop\Testing Folder Start\';
% version_control = false;


% OPTION 1: no data folder exists, full copy the folder
if ~isfolder(targetPath) % no current folder exists
    success = copyfile(startPath, targetPath);
    if ~success
        warndlg(['Error in copying ' startPath])
        return
    end
% OPTION 2: folder exists locally, check for mismatch between folder contents    
else 
    % Copy any files from the main folder to the target destination
    if version_control % update any files to the new version
        updateToNewerFile(startPath,targetPath);
    end
    success = copyMissingFiles(startPath, targetPath);
    if ~success
        disp(['Error copying: ' targetPath])
    end          
    % determine if any contents within 1 layer of subFOLDERS are missing
    subfolder_contents = dir(startPath);
    subfolderList = {subfolder_contents(:).name};
    % get list of folders
    subfolderNames = subfolderList([subfolder_contents(:).isdir]);

    % fold_idx = [subfolder_contents(:).isdir]; 
    % subfolders_idx = find(fold_idx==true);
    for jj = 1:length(subfolderNames)
        subfolder_name = subfolderNames{jj};
        % skip the empty folder holders
        if any(strcmpi(subfolder_name, {'.', '', '..'})) %TODO:  see how this works on mac 
            continue
        end
        start_subfolder = [startPath subfolder_name '/'];
        end_subfolder = [targetPath subfolder_name '/'];
        if version_control % update any files to the newest version
            updateToNewerFile(start_subfolder,end_subfolder);
        end
        success = copyMissingFiles(start_subfolder, end_subfolder); %copy any missing files in the folder
        if ~success
            disp(['Error copying: ' end_subfolder])
        end
    end
end

