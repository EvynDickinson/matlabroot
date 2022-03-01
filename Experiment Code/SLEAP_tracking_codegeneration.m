

%% Code to both TRACK and CONVERT to analysis file
% parameter inputs: 
clear
defaultModel = 'G:\My Drive\Jeanne Lab\DATA\Tracking\models';
model_1 = 'centered_instance_model';
model_2 = 'centroid_model';
% select the folder to analyze 
[~, folder] = getCloudPath(1);

% select which arenas you want to track:
rootDir = [folder '\'];
codeBlockPath = [rootDir 'TrackingCodeBlock.txt'];

% copy the current models into that folder if they don't already exist
if ~isfolder([rootDir '\' model_1]) && ~isfolder([rootDir '\' model_2])
    copyfile([defaultModel '\' model_1],[rootDir model_1])
    copyfile([defaultModel '\' model_2],[rootDir model_2])
end

% % text to set the cd to the desired folder: 
% setCD = ['cd ' rootDir];
% setCD = strrep(setCD, '\', '/');
% fid = fopen(codeBlockPath,'w');  
% fprintf(fid, setCD);
% fclose(fid);

%Select the complete experiments to process
list_dirs = dir([rootDir '*.avi']); %only matlab files
videoNames = {list_dirs(:).name};

% generate new text block for each selected video: 
fid = fopen(codeBlockPath,'a+'); 
for ii = 1:length(videoNames)
    vidName = videoNames{ii};
    trackStr = ['\n' 'sleap-track --tracking.tracker simple '...
                ' -m "centered_instance_model" -m "centroid_model" '...
                '"' vidName '"'];
    fprintf(fid, trackStr);
    predName = [vidName '.predictions.slp'];
    analName = [vidName(1:end-4) '.h5'];
    convertStr = ['\n' 'sleap-convert --format analysis -o '...
                '"' analName '" '... % analysis file
                '"' predName '"']; %prediciton file
    fprintf(fid, convertStr);
end
fclose(fid);
fprintf('\nDone\n')

%% Find and print out list of missing tracking or prediction files:
clear

% select the folder to analyze s
[~, folder] = getCloudPath(1);
arenaList = {'Arena A', 'Arena B', 'Arena C', 'Arena D'};
arenaSel = listdlg('ListString', arenaList);

for tt = 1:length(arenaSel)
    arena = arenaList{arenaSel(tt)};
    rootDir = [folder, '\' arena, '\'];
    codeBlockPath = [rootDir 'TrackingCodeBlock.txt'];
    vidPath = rootDir;
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
    predictions = vertcat(predictions{:}); % To remove nesting of cell array newA
    predictions = predictions(:,1);
    
    % Find and reprint the missing files:
    trackFail = find(ismember(videos, tracking)==false);
    precictFail = find(ismember(videos, predictions)==false);
    if ~isempty(trackFail) || ~isempty(precictFail)
        fid = fopen(codeBlockPath,'w');
        disp(['Missing files found in ' arena])
        fprintf(fid, '\n\nMISSED TRACKING/PREDICTION FILES:');
        if ~isempty(trackFail)
            %print out new files to be run
            for ii = 1:length(trackFail)
            vidName = videoNames{trackFail(ii)};
            trackStr = ['\n' 'sleap-track --tracking.tracker simple '...
                        ' -m "centered_instance_model" -m "centroid_model" '...
                        '"' vidName '"'];
            fprintf(fid, trackStr);
            end
        end
        if ~isempty(precictFail)
        %   print missing prediction files
            for ii = 1:length(precictFail)
                vidName = videoNames{precictFail(ii)};
                predName = [vidName '.predictions.slp'];
                analName = [vidName(1:end-4) '.h5'];
                convertStr = ['\n' 'sleap-convert --format analysis -o '...
                            '"' analName '" '... % analysis file
                            '"' predName '"']; %prediciton file 
                fprintf(fid, convertStr);
            end
        end
        fclose(fid);
    else
        disp(['No missing tracking files in ' arena])
    end
end


