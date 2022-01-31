
%      fig = cropCheck;


%% Code to both TRACK and CONVERT to analysis file
% parameter inputs: 
clear
defaultModel = 'G:\My Drive\Jeanne Lab\DATA\Tracking\models';
model_1 = 'centered_instance_model';
model_2 = 'centroid_model';
% select the folder to analyze 
[~, folder] = getCloudPath(1);

% select which arenas you want to track:
arenaList = {'Arena A', 'Arena B', 'Arena C', 'Arena D'};
arenaSel = listdlg('ListString', arenaList);

for tt = 1:length(arenaSel)
    arena = arenaList{arenaSel(tt)};
    rootDir = [folder '\' arena '\'];
    codeBlockPath = [rootDir 'TrackingCodeBlock.txt'];
    
    % Check the ROI fit before commiting to tracking
    cropCheck(folder, arena);
    temp = questdlg('Arena match?');
    if strcmp(temp,'No')
        disp(['skipped ' arena ' for bad ROI fit'])
        continue
    elseif strcmp(temp,'Cancel')
        return 
    end
    fig = gcf;
    close(fig)
    % copy the current models into that folder if they don't already exist
    if ~isfolder([rootDir '\' model_1]) && ~isfolder([rootDir '\' model_2])
        copyfile([defaultModel '\' model_1],[rootDir model_1])
        copyfile([defaultModel '\' model_2],[rootDir model_2])
    end
    
    % text to set the cd to the desired folder: 
    setCD = ['cd ' rootDir];
    setCD = strrep(setCD, '\', '/');
    fid = fopen(codeBlockPath,'w');  
    fprintf(fid, setCD);
    fclose(fid);
    
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
    fprintf(['\nFinished ' arena])
end
fprintf('\nDone\n')

%% Find and print out list of missing tracking or prediction files:
clear

% select the folder to analyze 
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


%% Select data file to make new folders for
%get base folder pathway
[baseFolder, folder] = getCloudPath(2); 

% Make a new video folder:
videosDirA = [baseFolder folder '\Arena A\'];
if ~isfolder(videosDirA); mkdir(videosDirA); end

videosDirB = [baseFolder folder '\Arena B\'];
if ~isfolder(videosDirB); mkdir(videosDirB); end

videosDirC = [baseFolder folder '\Arena C\'];
if ~isfolder(videosDirC); mkdir(videosDirC); end

videosDirD = [baseFolder folder '\Arena D\'];
if ~isfolder(videosDirD); mkdir(videosDirD); end



%% Old Code to both TRACK and CONVERT to analysis file
% % parameter inputs:
% clear
% defaultModel = 'G:\My Drive\Jeanne Lab\DATA\Tracking\models';
% model_1 = 'centered_instance_model';
% model_2 = 'centroid_model';
% 
% % select the folder to analyze 
% [~, folder] = getCloudPath(1);
% codeBlockPath = [folder '\Split Videos\TrackingCodeBlock.txt'];
% 
% % copy the current models into that folder
% switch questdlg('Copy models into selected folder?')
%     case 'Yes'
%         copyfile([defaultModel '\' model_1],[folder '\Split Videos\' model_1])
%         copyfile([defaultModel '\' model_2],[folder '\Split Videos\' model_2])
%     case 'No'
%     case 'Cancel'
%         return
% end
% 
% % text to set the cd to the desired folder: 
% setCD = ['cd ' folder '/Split Videos'];
% setCD = strrep(setCD, '\', '/');
% fid = fopen(codeBlockPath,'w');
% fprintf(fid, setCD);
% fclose(fid);
% 
% % pull video names from within the folder:
% 
% %Select the complete experiments to process
% list_dirs = dir([folder, '\Split Videos\*.avi']); %only matlab files
% videoNames = {list_dirs(:).name};
% indx = listdlg('ListString', videoNames, 'SelectionMode', 'Multiple');
% 
% % generate new text block for each selected video:
% fid = fopen(codeBlockPath,'a+');
% for ii = 1:length(indx)
%     vidName = videoNames{indx(ii)};
%     trackStr = ['\n' 'sleap-track --tracking.tracker simple '...
%                 ' -m "centered_instance_model" -m "centroid_model" '...
%                 '"' vidName '"'];
%     fprintf(fid, trackStr);
%     predName = [vidName '.predictions.slp'];
%     analName = [vidName(1:end-4) '.h5'];
%     convertStr = ['\n' 'sleap-convert --format analysis -o '...
%                 '"' analName '" '... % analysis file
%                 '"' predName '"']; %prediciton file
%     fprintf(fid, convertStr);
% end
% fclose(fid);


%% OLD Find and print out list of missing tracking or prediction files:
% clear
% 
% % select the folder to analyze 
% [~, folder] = getCloudPath(1);
% codeBlockPath = [folder '\Split Videos\TrackingCodeBlock.txt'];
% vidPath = [folder '\Split Videos\'];
% % Pull the full list of movies that should have been tracked:
% list_dirs = dir([vidPath, '*.avi']); %only matlab files
% videoNames = {list_dirs(:).name};
% videos = cellfun(@(x) strsplit(x, '.'), videoNames, 'UniformOutput', false);
% videos = vertcat(videos{:}); % To remove nesting of cell array newA
% videos = videos(:,1);
% % Pull list of SLP tracking files:
% list_dirs = dir([vidPath, '*.slp']); %only matlab files
% trackNames = {list_dirs(:).name};
% tracking = cellfun(@(x) strsplit(x, '.'), trackNames, 'UniformOutput', false);
% tracking = vertcat(tracking{:}); % To remove nesting of cell array newA
% tracking = tracking(:,1);
% % Pull list of prediction files:
% list_dirs = dir([vidPath, '*.h5']); %only matlab files
% predictionNames = {list_dirs(:).name};
% predictions = cellfun(@(x) strsplit(x, '.'), predictionNames, 'UniformOutput', false);
% predictions = vertcat(predictions{:}); % To remove nesting of cell array newA
% predictions = predictions(:,1);
% 
% % Find and reprint the missing files:
% trackFail = find(ismember(videos, tracking)==false);
% precictFail = find(ismember(videos, predictions)==false);
% if ~isempty(trackFail) || ~isempty(precictFail)
%     fid = fopen(codeBlockPath,'a+');
%     disp('Missing files found')
%     fprintf(fid, '\n\nMISSED TRACKING/PREDICTION FILES:');
%     if ~isempty(trackFail)
%         %print out new files to be run
%         for ii = 1:length(trackFail)
%         vidName = videoNames{trackFail(ii)};
%         trackStr = ['\n' 'sleap-track --tracking.tracker simple '...
%                     ' -m "centered_instance_model" -m "centroid_model" '...
%                     '"' vidName '"'];
%         fprintf(fid, trackStr);
%         end
%     end
%     if ~isempty(precictFail)
%     %   print missing prediction files
%         for ii = 1:length(precictFail)
%             vidName = videoNames{precictFail(ii)};
%             predName = [vidName '.predictions.slp'];
%             analName = [vidName(1:end-4) '.h5'];
%             convertStr = ['\n' 'sleap-convert --format analysis -o '...
%                         '"' analName '" '... % analysis file
%                         '"' predName '"']; %prediciton file
%             fprintf(fid, convertStr);
%         end
%     end
%     fclose(fid);
% else
%     disp('No missing tracking files')
% end




%% Write code to track flies 
% parameter inputs:
clear
defaultModel = 'G:\My Drive\Jeanne Lab\DATA\Tracking\models';
model_1 = 'centered_instance_model';
model_2 = 'centroid_model';

% select the folder to analyze 
[~, folder] = getCloudPath(1);
codeBlockPath = [folder '\TrackingCodeBlock.txt'];

% copy the current models into that folder
switch questdlg('Copy models into selected folder?')
    case 'Yes'
        copyfile([defaultModel '\' model_1],[folder '\' model_1])
        copyfile([defaultModel '\' model_2],[folder '\' model_2])
    case 'No'
    case 'Cancel'
        return
end

% % text to set the cd to the desired folder: 
setCD = ['cd ' folder];
setCD = strrep(setCD, '\', '/');
fid = fopen(codeBlockPath,'w');
fprintf(fid, setCD);
fclose(fid);

% pull video names from within the folder:

% Select the complete experiments to process
list_dirs = dir([folder, '\*.avi']); %only matlab files
videoNames = {list_dirs(:).name};
indx = listdlg('ListString', videoNames, 'SelectionMode', 'Multiple');

% generate new text block for each selected video:
fid = fopen(codeBlockPath,'a+');
for ii = 1:length(indx)
    vidName = videoNames{indx(ii)};
    trackStr = ['\n' 'sleap-track --tracking.tracker simple '...
                ' -m "centered_instance_model" -m "centroid_model" '...
                '"' vidName '"'];
    fprintf(fid, trackStr);
end
fclose(fid);



%% Code to convert the prediction files to data files:
% terminal name adjustment .avi.predictions.slp
% sleap-convert --format analysis -o "session1.h5" "session1.predictions.slp"
clear
% Select a folder and check for prediction files, then convert those
[~, folder] = getCloudPath(1);

% Select the complete experiments to process
list_dirs = dir([folder, '\*.slp']); %only matlab files
predictions = {list_dirs(:).name};
indx = listdlg('ListString', predictions, 'SelectionMode', 'Multiple', 'ListSize', [300 300]);

% Generate code text for each predicted file: 
codeBlockPath = [folder '\ConversionCodeBlock.txt'];
fid = fopen(codeBlockPath,'w');
fprintf(fid, 'Conversion code \n');
fclose(fid);

% Add conversion code line by line:
fid = fopen(codeBlockPath,'a+');
for ii = 1:length(indx)
    predName = predictions{indx(ii)};
    newstr = split(predName, '.');
    analName = [newstr{1} '.h5'];
    trackStr = ['\n' 'sleap-convert --format analysis -o '...
                '"' analName '" '... % analysis file
                '"' predName '"']; %prediciton file
    fprintf(fid, trackStr);
end
fclose(fid);




















