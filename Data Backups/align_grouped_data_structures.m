

%% Comparison between grouped data structures


clear
clc

%%  Select the data to be synchronized between the drives (can be between any two drives)


 
serverPath = getDataPath(2,2);
localPath = getDataPath(2,0, 'Select the LOCAL data folder:');


% Synchronization should occur based on the most recent file timestamp,
% if two copies exist -- ask for everything so that nothing gets copied
% that may have been targeted for deletion



switch questdlg('Select the direction you want to align the data:','Raw Data Alignment',...
                             'Download from Server', 'Upload to Server', 'Cancel', 'Download from Server')
    case 'Download from Server'
                FROM_base = serverPath;
                TO_base = localPath;
    case 'Upload to Server'
                FROM_base = localPath;
                TO_base = serverPath;
    case {'Cancel',''}
        disp('Data synchronization canceled')
        return
end

% select the data files to download:
folderSelected = selectFolder(FROM_base,'mulitple'); 
nDir = length(folderSelected);

%% upload any new data from the local raw data folder to the server folder

% check that the folder and data doesn't already exist in the target drive:
for i = 1:nDir
    disp(['Aligning data from '  folderSelected{i} ' ...'])
    targetPath = [TO_base, folderSelected{i} '/'];
    startPath = [FROM_base, folderSelected{i} '/'];
    if ~isfolder(targetPath) % no current folder exists
        % copy the full folder contents to the local destination
        copyStatus = copyfile(startPath, targetPath);
        if ~copyStatus
            warndlg(['Error in copying ' folderSelected{i}])
            return
        end
    else % folder exists locally, check for mismatch between folder contents
        % Copy any files from the main folder to the target destination
        success = copyMissingFiles(startPath, targetPath);
        if ~success
            disp(['Error copying: ' targetPath])
        end

        % determine if any contents within subFOLDERS are missing
        subfolder_contents = dir(startPath);
        subfolderList = {subfolder_contents(:).name};
        fold_idx = [subfolder_contents(:).isdir]; 
        subfolders_idx = find(fold_idx==true);
        for jj = 1:length(subfolders_idx)
            subfolder_name = subfolderList{subfolders_idx(jj)};
            if any(strcmpi(subfolder_name, {'.', '', '..'})) %skip the empty folder holders
                continue
            end
            start_subfolder = [startPath subfolder_name '/'];
            end_subfolder = [targetPath subfolder_name '/'];
            success = copyMissingFiles(start_subfolder, end_subfolder); %copy any missing files in the folder
            if ~success
                disp(['Error copying: ' end_subfolder])
            end
        end
    end
end

disp('Completed all data transfers')

        

    


