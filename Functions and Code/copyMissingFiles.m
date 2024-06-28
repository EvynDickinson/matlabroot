

function results = copyMissingFiles(startPath, destPath)
% results = copyMissingFiles(startPath, destPath)
% copy any missing contents between the start directory and destination
% directory
%
% ES Dickinson, Yale University 2024


% get the names of the folders in the single trial portion of the local drive:
destination_contents = dir(destPath);
destFolders = {destination_contents(:).name};

% get the names of the folders in the single trial portion of the svalbard server:
startList = dir(startPath);
startFolders = {startList(:).name};

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
        % if the file is a folder, make sure that it is copying the
        % contents of the folder into a new folder of the same name! 
        if isfolder(start_location)
            statusList(idx) = copyfile(start_location, [destPath missingdestinationFiles{idx}]);
        else
            statusList(idx) = copyfile(start_location, destPath);
        end
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

