
%% Write code to track flies 

defaultModel = 'G:\My Drive\Jeanne Lab\DATA\Tracking\models';

% select the folder to analyze 
[baseFolder, folder] = getCloudPath(1);

% copy the current models into that folder
copyfile(defaultModel,[folder '\models'])




% text to set the cd to the desired folder: 
setCD = ['cd ' folder];
setCD = strrep(setCD, '\', '/');


% pull video names from within the folder:

% Select the complete experiments to process
list_dirs = dir([folder, '\*.avi']); %only matlab files
videoNames = {list_dirs(:).name};
indx = listdlg('ListString', videoNames, 'SelectionMode', 'Multiple');

for ii = 1:length(indx)
    
    
end
% set tracking text:
trackStr = ['sleap-track "' 






% 
sleap-track "PlantYeastChoice_N1_2.avi" 
--frames 1-10 
--tracking.tracker simple 
--model "models\centered_instance_model" 
--model "models\centroid_model"