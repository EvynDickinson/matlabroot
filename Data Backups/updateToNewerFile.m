


%% Compare time stamps on items found in BOTH locations to copy the newest versions to the destination folder

function results = updateToNewerFile(startPath, destPath, bidirectional) 
% results = updateToNewerFile(startPath, destPath, bidirectional)
% Update files of the same name and type to the newest version
% input the paths to the folders with the files to compare
%
% bidirectional default is false; 
% if bidirectional is true, then the newest version in either directory is
% copied over the older version, if it is false, then only a newer file in the
% start directory is written into the destination folder
% 
% ES Dickinson

if nargin<3 
    bidirectional = false;
end

% Get the list of files in both source folders
files1 = dir(startPath);
files2 = dir(destPath);

% Remove directories from the list of files
files1 = files1(~[files1.isdir]);
files2 = files2(~[files2.isdir]);

% Get the names of the files
fileNames1 = {files1.name}';
fileNames2 = {files2.name}';

% Find the overlapping files
[overlappingFiles, idx1, idx2] = intersect(fileNames1, fileNames2);
fail_count = false;
% Iterate through the overlapping files
for i = 1:length(overlappingFiles)
    file1 = files1(idx1(i)); % start files
    file2 = files2(idx2(i)); % destrination files
    
    % Compare the modification dates
    if file1.datenum > file2.datenum
        % Copy file1 to the target folder if it's newer
        success = copyfile(fullfile(startPath, file1.name), fullfile(destPath, file1.name));
        if ~success
            disp(['Error updating ' file1.name ' from ' startPath])
            fail_count = [fail_count; true];
        end
    % Copy file2 to the target folder if it's newer    
    elseif file1.datenum < file2.datenum && bidirectional 
        success = copyfile(fullfile(destPath, file2.name), fullfile(startPath, file2.name));
        if ~success
            disp(['Error updating ' file2.name ' from ' destPath])
            fail_count = [fail_count; true];
        end
    end
end

% check for number of failed updating:
n_fail = sum(fail_count);
if any(fail_count)
    results = false;
    disp(["Number of failed updates: " num2str(n_fail)])
else
    results = true;
end

% Inform the user that the operation is complete
disp('File comparison and updating complete.');





