
% COPY PROCESSED DATA OVER TO THE SERVER FROM EXTERNAL STORAGE DRIVES

% 1) copy full trial data folder over to the trial data folder on the server
% 2) Create a date folder in the Courtship Videos folder
% 3) For each folder in the date folder: 
    % a. copy all folders into the new date folder 
            % **double check that there is not a 'raw' folder if there is, flag for deletion **
    % b. move all video files into a new folder called 'Videos'
    % c. copy all loose files into the new date folder on the server

%%  Get the original data folder

% Pull list of files that haven't been (according to the excel sheet)
% copied to the server
clear; clc

[excelfile, Excel, xlFile] = load_HighResExperiments;
paths = getPathNames;

% Find files that can be run
[StartRoot, StartName] = getDataPath(5, 0, 'Where is the raw data?'); 

% Where is the data going?
[EndRoot, EndName] = getDataPath(5, 6, 'Where do you want the processed data to go?');

% what files still haven't been transfered to the server from the selected location?
if strcmp(StartName,'data storage')
    basestr = 'storage';
end
rows1 = ~cellfun('isempty', regexp(string(excelfile(:, Excel.storagedrive)), basestr, 'ignorecase'));
rows0 = strcmpi(string(excelfile(:, Excel.groupready)), 'Y');
% list of files that haven't been uploaded to the server yet
rows2 = cellfun(@(x) isempty(x) || (isnumeric(x) && isnan(x)), excelfile(:, Excel.processed_data_on_server));
% list of files that have been processed far enough to be ready for processed data transfer
rows3 = cellfun(@(x) isempty(x) || (isnumeric(x) && isnan(x)), excelfile(:, Excel.trialID));
% find overlapping files that need to be transferred to the server: 
loc = find(rows1 & rows2 & ~rows3 & rows0);
expList = excelfile(loc, Excel.trialID);

% Auto limit options to those on the drive that was selected: 
currOptions = dir([StartRoot paths.courtship 'Trial Data/']); % pull list of folder names that are in the selected start drive
currOptions = {currOptions(3:end).name};
[SharedTrials,~,ia] = intersect(currOptions, expList); % find the folders that are present in both locations
loc = loc(ia); % restrict excel locations to those are on the selected drive

% cancel out if there are no available options: 
if isempty(loc)
    warndlg(['There are no files that are ready for data transfer at this point on the drive you selected.' ...
         'Try checking that the High Res excel sheet has updated values in the storage drive and group ready columns'...
         'aka go add Y!']);
    return
end

% offer up the list for approval
idx = listdlg('PromptString','Select experiments to copy:','ListString',SharedTrials,...
                    'SelectionMode','multiple','InitialValue',1:length(loc),'ListSize',[350,300]);
sel_loc = loc(idx); % final selected list of experiments to copy over

%% Copy data: 

% copy Trial Data to server
for i = 1:length(sel_loc)
    idx = sel_loc(i);
    folderName = excelfile{idx, Excel.trialID};
    originalFolder = [StartRoot paths.courtship 'Trial Data/' folderName];
    destinationFolder = [EndRoot paths.courtship 'Trial Data/' folderName];
    copyfile(originalFolder, destinationFolder)
end

% Copy regular Date folder files
for i = 1:length(sel_loc)
    idx = sel_loc(i);
    folderDate = excelfile{idx, Excel.date};
    folderName = excelfile{idx, Excel.expID};
    originalFolder = [StartRoot paths.courtship folderDate '\' folderName '\'];
    destinationFolder = createFolder([EndRoot paths.courtship folderDate  '\' folderName '\']);
    
    % check if there is raw data still in existence:
    if isfolder([originalFolder 'raw']) || ~isempty(dir([originalFolder '*file*.mat']))
        h = warndlg({'Please delete the raw data folder in '; originalFolder});
        uiwait(h)
    end

    % get list of all folders in the base folder: 
    folderList = dir(originalFolder);
    folderList(1:2) = [];
    folder_idx = find([folderList(:).isdir]); 
    folderNames = {folderList(folder_idx).name};
    folderNames(strcmp(folderNames,'AVI_Files')) = []; % skip the processed AVI file folder
    folderNames(strcmp(folderNames,'raw')) = []; % skip the raw data folder if it wasn't already deleted

    for f = 1:length(folderNames)
        subfolderName = folderNames{f};
        copyfile([originalFolder subfolderName], [destinationFolder subfolderName])
    end
    
    % move all video files into their own folder
    avifiles = dir([originalFolder '*.avi']);
    if ~isempty(avifiles)
        aviFolder = createFolder([originalFolder 'AVI_Files']);
        movefile(fullfile(originalFolder, '*.avi'), aviFolder); 
    end

    % move all data files into their own folder
    dataTypes = {'*.mat', '*.slp','*.h5', '*.py', '*.png','*.csv'};
    dataFolder = createFolder([originalFolder 'data']);
    for f = 1:length(dataTypes)
        if ~isempty(dir(fullfile(originalFolder, dataTypes{f})))
            movefile(fullfile(originalFolder, dataTypes{f}), dataFolder);
        end
    end

    % copy data files to the new folder
    copyfile(dataFolder, destinationFolder);

    % write update to processed data in excel:
     if ~isExcelFileOpen(xlFile,true)
        try writecell({EndName},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.processed_data_on_server) num2str(idx)]);
        catch 
            disp('couldn''t write to excel:  manually update')
        end
    else
        disp('Couldn''t write to excel sheet for:')
        disp(originalFolder)
    end
disp(['successfully transferred: ', dataFolder])
end































