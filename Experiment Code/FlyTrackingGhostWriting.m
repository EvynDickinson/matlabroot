
%% Write code to track flies 
% parameter inputs:
defaultModel = 'G:\My Drive\Jeanne Lab\DATA\Tracking\models';
model_1 = 'centered_instance_model';
model_2 = 'centroid_model';
codeBlockPath = 'G:\My Drive\Jeanne Lab\DATA\Tracking\TrackingCodeBlock.txt';

% select the folder to analyze 
[baseFolder, folder] = getCloudPath(1);

% copy the current models into that folder
switch questdlg('Copy models into selected folder?')
    case 'Yes'
        copyfile(defaultModel,[folder '\models'])
    case 'No'
    case 'Cancel'
        return
end


% text to set the cd to the desired folder: 
setCD = ['cd ' folder];
setCD = strrep(setCD, '\', '/');

writematrix(setCD, codeBlockPath, 'WriteMode', 'overwrite')


% pull video names from within the folder:

% Select the complete experiments to process
list_dirs = dir([folder, '\*.avi']); %only matlab files
videoNames = {list_dirs(:).name};
indx = listdlg('ListString', videoNames, 'SelectionMode', 'Multiple');

% generate new text block for each selected video:
codeBlock = [];
for ii = 1:length(indx)
    vidName = videoNames{indx(ii)};
    trackStr = ['sleap-track "' vidName '" --tracking.tracker simple --model "' ...
                'models/' model_1 '" --models "models/' model_2 '"'];
    writematrix(trackStr, codeBlockPath, 'WriteMode', 'append')
end








% 
% sleap-track "PlantYeastChoice_N1_2.avi" 
% --frames 1-10 
% --tracking.tracker simple 
% --model "models\centered_instance_model" 
% --model "models\centroid_model"