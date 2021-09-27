


%% Code to both TRACK and CONVERT to analysis file
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

% text to set the cd to the desired folder: 
setCD = ['cd ' folder];
setCD = strrep(setCD, '\', '/');
fid = fopen(codeBlockPath,'w');
fprintf(fid, setCD);
fclose(fid);

% pull video names from within the folder:

%Select the complete experiments to process
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
    predName = [vidName '.predictions.slp'];
    analName = [vidName(1:end-4) '.h5'];
    convertStr = ['\n' 'sleap-convert --format analysis -o '...
                '"' analName '" '... % analysis file
                '"' predName '"']; %prediciton file
    fprintf(fid, convertStr);
end
fclose(fid);



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




















