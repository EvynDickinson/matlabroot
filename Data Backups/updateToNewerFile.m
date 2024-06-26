


%% Compare time stamps on items found in BOTH locations to copy the newest versions to the destination folder

function results = updateToNewerFile(startPath, destPath) %TODO -- add a switch for bidirectional?

destPath = uigetdir();
startPath = uigetdir();

% get the names of the folders in the single trial portion of the local drive:
destination_contents = dir(destPath);
destFolders = {destination_contents(:).name};

% get the names of the folders in the single trial portion of the svalbard server:
startList = dir(startPath);
startFolders = {startList(:).name};

% 1) copy anything that doesn't exist in the first location
% 2) check dates for anything that exists in both locations

% determine if any FILES are missing on the local folder in the main folder:
missingdestinationFiles = setdiff(startFolders,destFolders); % what is in the start folder but not in the dest

[overlappingdestinationFiles, idxstart, idxDest] = intersect(startFolders,destFolders); % found in both locations
 
% get timestamps for the files that are overlapping


    

                % Define the source folders and the target folder
                % sourceFolder1 = 'path_to_folder1';
                % sourceFolder2 = 'path_to_folder2';
                % targetFolder = 'path_to_target_folder';
                
                % Get the list of files in both source folders
                files1 = dir(startPath);
                files2 = dir(destPath);
                
                % Remove directories from the list of files
                files1 = files1(~[files1.isdir]);
                files2 = files2(~[files2.isdir]);
                
                % Get the names of the files
                fileNames1 = {files1.name};
                fileNames2 = {files2.name};
                
                % Find the overlapping files
                [overlappingFiles, idx1, idx2] = intersect(fileNames1, fileNames2);
                
                % Iterate through the overlapping files
                for i = 1:length(overlappingFiles)
                    file1 = files1(idx1(i));
                    file2 = files2(idx2(i));
                    
                    % Compare the modification dates
                    if file1.datenum > file2.datenum
                        % Copy file1 to the target folder if it's newer
                        copyfile(fullfile(sourceFolder1, file1.name), fullfile(targetFolder, file1.name));
                    else
                        % Copy file2 to the target folder if it's newer
                        copyfile(fullfile(sourceFolder2, file2.name), fullfile(targetFolder, file2.name));
                    end
                end
                
                % Inform the user that the operation is complete
                disp('File comparison and copying complete.');




% determine if any FILES are missing on the local folder in the main folder:
missingdestinationFiles = setdiff(startFolders,destFolders); % what is in the start folder but not in the dest

if ~isempty(missingdestinationFiles)
    statusList = false(length(missingdestinationFiles),1); % fill in blanks with copy status
    for idx = 1:length(missingdestinationFiles)
        start_location = [startPath missingdestinationFiles{idx}];
        if any(strcmpi(start_location, {'.', '', '..'})) %skip the empty folder holders
            statusList(idx) = true;
            continue
        end
        statusList(idx) = copyfile(start_location, destPath);
    end
    if ~all(statusList)
        warndlg({'Copy error moving: '; ['From: ' start_location]; ['To: ' destPath]})
        results = false;
    else
        results = true;
    end
else % no missing files
    results = true;
end

