

%% SELECT EXPERIMENT TO TRANSFER
clear  

%raw data folder: 
date_today = strrep(datestr(datetime,'mm-dd-yyyy'),'-','.');
start_dir = 'F:\Evyn\DATA\'; %'C:\Users\jeannelab\Documents\Evyn\DATA\';

%get base folder pathway  
baseFolder = 'S:\Evyn\Raw Data\DATA'; % G:\My Drive\Jeanne Lab\DATA\';  

%select folder date    
list_dirs = dir(start_dir);
for i = 3:length(list_dirs)
    folderNames{i-2} = list_dirs(i).name;
end
folderNames = ['Today', flip(folderNames)];

indx = listdlg('ListString', folderNames, 'SelectionMode', 'Single');
if strcmpi(folderNames{indx}, 'Today')==true
    dir_sel = date_today;
else
    dir_sel = folderNames{indx};
end 
folder = fullfile(start_dir, dir_sel);
targetDir = [baseFolder dir_sel];

% COPY TEMPERATURE LOG FOR SECOND EXPERIMENT FILE
% 1) get the list of experiments in the date folder
videoList = dir([folder '\*_1.avi']);
[exp_names, start_times] = deal([]);
for i = 1:length(videoList)
    exp_names{i} =  videoList(i).name(1:end-6);
    exp_start_time{i} = videoList(i).date(end-7:end-3);
end

% 2) get the list of temperature log files
tempLogList = dir([folder '\*_RampLog.csv']);
% get experiment names that have a temp log file
tempLog_names = [];
for i = 1:length(tempLogList)
    tempLog_names{i} =  tempLogList(i).name(1:end-12);
end

% 3) Find the experiments that do not have a temperature log
idx = 0;
for i = 1:length(exp_names)
    if ~any(strcmp(exp_names{i}, tempLog_names)) %no matches
       %note the unmatched experiment name
       idx = idx + 1;
       no_log_list{idx} = exp_names{i};
    end
end

% 4) Find any experiments that match the start time of the missing log experiments:
if idx > 0
    for i = 1:idx %for each unmatched experiment
        sel_start_time =  exp_start_time{strcmp(no_log_list{i},exp_names)}; % start time of the selected experiment
        search_loc = find(~strcmp(no_log_list{i},exp_names)); % index of other experiments
        for ii = search_loc
              if str2double(strrep(sel_start_time,':','')) == str2double(strrep(exp_start_time{ii},':','')) % start time match
                    source_file = [folder '\' exp_names{ii} '_RampLog.csv'];
                    destination_file = [folder '\' no_log_list{i} '_RampLog.csv'];
                    copyStatus = copyfile(source_file, destination_file);
                    if ~copyStatus
                        h = warndlg([exp_names{ii} ' not copied']);
                        uiwait(h)
                    end
                    continue
              end
        end
    end
end
                 
%% Move videos from computer drive to the server

% Move folders to google drive:
copyfile(folder, targetDir)

% Make new video folders for each Arena:
videosDirA = [targetDir '\Arena A\'];
if ~isfolder(videosDirA); mkdir(videosDirA); end

videosDirB = [targetDir '\Arena B\'];
if ~isfolder(videosDirB); mkdir(videosDirB); end

videosDirC = [targetDir '\Arena C\'];
if ~isfolder(videosDirC); mkdir(videosDirC); end

videosDirD = [targetDir '\Arena D\'];
if ~isfolder(videosDirD); mkdir(videosDirD); end

% Copy tracking models over to new folder: 
trackingDir = [baseFolder 'Tracking'];
copyfile(trackingDir, targetDir)

% Add second & third batch processing folder
batch2Dir = [targetDir '\Batch 2\'];
if ~isfolder(batch2Dir); mkdir(batch2Dir); end
copyfile(trackingDir, batch2Dir)

% Add second & third batch processing folder
batch3Dir = [targetDir '\Batch 3\'];
if ~isfolder(batch3Dir); mkdir(batch3Dir); end
copyfile(trackingDir, batch3Dir)


fprintf('done')


























% %% Copy folders for models over to selected location
% % Move the 
% clear
% % select the folder to analyze 
% [basePath, folder] = getCloudPath(1);
% 
% % parameter inputs: 
% model_1 = 'centered_instance_model';
% model_2 = 'centroid_model';
% batchCode = 'batch.py';
% trackingDir = [basePath 'Tracking\'];
% 
% % Copy model and program into selected folder
% copyfile([trackingDir model_1], folder)
% copyfile([trackingDir model_2], folder)
% copyfile([trackingDir batchCode], folder)
% 
% fprintf('done')





