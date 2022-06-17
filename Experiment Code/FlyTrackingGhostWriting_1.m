

% Find missing tracked files from a folder:

clear
% select the folder to analyze 
[basePath, folder] = getCloudPath(2);
vidPath = [basePath folder '\'];

% Pull the full list of movies that should have been tracked:
list_dirs = dir([vidPath, '*.avi']); %only matlab files
videoNames = {list_dirs(:).name};
videos = cellfun(@(x) strsplit(x, '.'), videoNames, 'UniformOutput', false);
videos = vertcat(videos{:}); % To remove nesting of cell array newA
videos = videos(:,1);
% Pull list of SLP tracking files:
list_dirs = dir([vidPath, '*.slp']); %only matlab files
trackNames = {list_dirs(:).name};
tracking = cellfun(@(x) strsplit(x, '.'), trackNames, 'UniformOutput', false);
tracking = vertcat(tracking{:}); % To remove nesting of cell array newA
tracking = tracking(:,1);
% Pull list of prediction files:
list_dirs = dir([vidPath, '*.h5']); %only matlab files
predictionNames = {list_dirs(:).name};
predictions = cellfun(@(x) strsplit(x, '.'), predictionNames, 'UniformOutput', false);
if isempty(predictions)
    % first round of tracking didn't even finish...
    % move all videos that haven't been predicted move over to a new folder
end % TODO

for ii = 1:size(predictions,2)
    new_predictions{ii,1} = predictions{1,ii}{1}; % To remove nesting of cell array
end

% Find and reprint the missing files:
trackFail = find(ismember(videos, tracking)==false);
precictFail = find(ismember(videos, new_predictions)==false);
if ~isempty(trackFail) || ~isempty(precictFail)
    fid = fopen(codeBlockPath,'w');
    disp('Missing files found')
    fprintf(fid, 'import os\nimport subprocess\n');
    
    if ~isempty(trackFail)
        %print out new files to be run
        for ii = 1:length(trackFail)
        vidName = videoNames{trackFail(ii)};
        trackStr = [tracking_start vidName tracking_end '\n'];
        fprintf(fid, trackStr);
        end
    end
    if ~isempty(precictFail)
    %   print missing prediction files
        for ii = 1:length(precictFail)
            vidName = videoNames{precictFail(ii)};
            predName = [vidName '.predictions.slp'];
            analName = [vidName(1:end-4) '.h5'];
            convertStr = [convert_start analName convert_mid predName convert_end '\n']; %prediciton file 
            fprintf(fid, convertStr);
        end
    end
    fclose(fid);
else
    disp(['No missing tracking files in ' arena])
end





%% Tracking models and code copy:
clear
% select the folder to analyze 
[basePath, folder] = getCloudPath(1);
trackingDir = [basePath 'Tracking'];
copyfile(trackingDir, folder)

% TODO: update to partition across multiple videos & find missing videos,
% move those into a new folder and track them...

% % how many partitions? 
% num_it = str2double(cell2mat(inputdlg('How many tracking partitions?')));
% % how many videos total in the folder?
% dirc = dir([folder '\*avi']);
% nvideos = length(dirc);
% groupSize = floor(nvideos/num_it);
% % copy videos into new folders and move models into each
% for ii = 1:num_it
%     loc = ii:ii+groupSize;
%     % make new folder

% import os
% import subprocess
% # Run tracking
% subprocess.call('(for %i in (*.avi) do sleap-track ''%i'' --tracking.tracker simple --tracking.similarity centroid -m "centered_instance_model" -m "centroid_model")', shell=True) ##--batch_size 4
% # Convert tracking to analysis file
% subprocess.call('(for %i in (*.slp) do sleap-convert -o %i.h5 --format analysis ''%i'')', shell=True)

%% Find missing tracking and prediction files and call them
clear

tracking_start = 'subprocess.call(''(sleap-track "';
tracking_end = '" --tracking.tracker simple --tracking.similarity centroid -m "centered_instance_model" -m "centroid_model")'', shell=True)';

convert_start = 'subprocess.call(''(sleap-convert -o "';  
convert_mid = '"  --format analysis "';
convert_end = '")'', shell=True)';


% select the folder to analyzes
[~, folder] = getCloudPath(1);
codeBlockPath = [folder '\TrackingCodeBlock.txt'];
vidPath = [folder '\'];
% Pull the full list of movies that should have been tracked:
list_dirs = dir([vidPath, '*.avi']); %only matlab files
videoNames = {list_dirs(:).name};
videos = cellfun(@(x) strsplit(x, '.'), videoNames, 'UniformOutput', false);
videos = vertcat(videos{:}); % To remove nesting of cell array newA
videos = videos(:,1);
% Pull list of SLP tracking files:
list_dirs = dir([vidPath, '*.slp']); %only matlab files
trackNames = {list_dirs(:).name};
tracking = cellfun(@(x) strsplit(x, '.'), trackNames, 'UniformOutput', false);
tracking = vertcat(tracking{:}); % To remove nesting of cell array newA
tracking = tracking(:,1);
% Pull list of prediction files:
list_dirs = dir([vidPath, '*.h5']); %only matlab files
predictionNames = {list_dirs(:).name};
predictions = cellfun(@(x) strsplit(x, '.'), predictionNames, 'UniformOutput', false);
for ii = 1:size(predictions,2)
    new_predictions{ii,1} = predictions{1,ii}{1}; % To remove nesting of cell array
end

% Find and reprint the missing files:
trackFail = find(ismember(videos, tracking)==false);
precictFail = find(ismember(videos, new_predictions)==false);
if ~isempty(trackFail) || ~isempty(precictFail)
    fid = fopen(codeBlockPath,'w');
    disp('Missing files found')
    fprintf(fid, 'import os\nimport subprocess\n');
    
    if ~isempty(trackFail)
        %print out new files to be run
        for ii = 1:length(trackFail)
        vidName = videoNames{trackFail(ii)};
        trackStr = [tracking_start vidName tracking_end '\n'];
        fprintf(fid, trackStr);
        end
    end
    if ~isempty(precictFail)
    %   print missing prediction files
        for ii = 1:length(precictFail)
            vidName = videoNames{precictFail(ii)};
            predName = [vidName '.predictions.slp'];
            analName = [vidName(1:end-4) '.h5'];
            convertStr = [convert_start analName convert_mid predName convert_end '\n']; %prediciton file 
            fprintf(fid, convertStr);
        end
    end
    fclose(fid);
else
    disp(['No missing tracking files in ' arena])
end


%% Rename the tracking file names: TODO

% select the folder to analyzes
[~, folder] = getCloudPath(1);
codeBlockPath = [folder '\TrackingCodeBlock.txt'];
vidPath = [folder '\'];
% Pull the full list of movies that should have been tracked:
list_dirs = dir([vidPath, '*.avi']); %only matlab files
videoNames = {list_dirs(:).name};
videos = cellfun(@(x) strsplit(x, '.'), videoNames, 'UniformOutput', false);
videos = vertcat(videos{:}); % To remove nesting of cell array newA
videos = videos(:,1);
% Pull list of SLP tracking files:
list_dirs = dir([vidPath, '*.slp']); %only matlab files
trackNames = {list_dirs(:).name};
tracking = cellfun(@(x) strsplit(x, '.'), trackNames, 'UniformOutput', false);
tracking = vertcat(tracking{:}); % To remove nesting of cell array newA
tracking = tracking(:,1);
% Pull list of prediction files:
list_dirs = dir([vidPath, '*.h5']); %only matlab files
predictionNames = {list_dirs(:).name};
predictions = cellfun(@(x) strsplit(x, '.'), predictionNames, 'UniformOutput', false);
for ii = 1:size(predictions,2)
    new_predictions{ii,1} = predictions{1,ii}{1}; % To remove nesting of cell array
end





























