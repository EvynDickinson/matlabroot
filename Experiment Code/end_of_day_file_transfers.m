

%% Move videos from computer drive to google drive
clear
%raw data folder: 
h = warndlg('Turn off IR lights on rig');
uiwait(h)
date_today = strrep(datestr(datetime,'mm-dd-yyyy'),'-','.');
start_dir = 'C:\Users\jeannelab\Documents\Evyn\DATA\';

%get base folder pathway  
switch getenv('COMPUTERNAME')
    case 'DENALI'
        baseFolder = 'E:\My Drive\Jeanne Lab\DATA\';
    case 'TOGIAK'
        baseFolder = 'G:\My Drive\Jeanne Lab\DATA\';
    case 'ACADIA'
        baseFolder = 'G:\My Drive\Jeanne Lab\DATA\';
end

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

fprintf('done')



%% Copy folders for models over to selected location
% Move the 
clear
% select the folder to analyze 
[basePath, folder] = getCloudPath(1);

% parameter inputs: 
model_1 = 'centered_instance_model';
model_2 = 'centroid_model';
batchCode = 'batch.py';
trackingDir = [basePath 'Tracking\'];

% Copy model and program into selected folder
copyfile([trackingDir model_1], folder)
copyfile([trackingDir model_2], folder)
copyfile([trackingDir batchCode], folder)

fprintf('done')





