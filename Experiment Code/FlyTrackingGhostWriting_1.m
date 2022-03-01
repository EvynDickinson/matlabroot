
%      fig = cropCheck;
%% Tracking models and code copy:
clear
% select the folder to analyze 
[basePath, folder] = getCloudPath(1);

% parameter inputs: 
defaultModel = 'G:\My Drive\Jeanne Lab\DATA\Tracking\models';
model_1 = 'centered_instance_model';
model_2 = 'centroid_model';
batchCode = 'batch.py';
trackingDir = [basePath 'Tracking\'];

% Copy model and program into selected folder
copyfile([trackingDir model_1], folder)
copyfile([trackingDir model_2], folder)
copyfile([trackingDir batchCode], folder)




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




