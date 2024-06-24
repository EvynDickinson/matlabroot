
clear
clc

%% download raw data from server for tracking and analysis


serverPath = getDataPath(2,2);
disp('Select local destination for raw data')
localPath = getDataPath(2,0);

% select the data files to download:
folderSelected = selectFolder(serverPath,'mulitple'); 
nDir = length(folderSelected);


% check that the folder and data doesn't already exist in the target drive:
for i = 1:nDir
    targetPath = [localPath, folderSelected{i} '/'];
    if ~isfolder(targetPath) % no current folder exists
        disp(['Copying '  folderSelected{i} ' ...'])
        % copy the full folder contents to the local destination
        copyStatus = copyfile([serverPath, folderSelected{i} '/'], targetPath);
        if ~copyStatus
            warndlg(['Error in copying ' folderSelected{i}])
            return
        end

    else % folder exists locally, check for mismatch between folder contents
        local_contents = dir(targetPath); % TODO HERE...
        % copy content differences
        % check that copy was successful
    end
end

