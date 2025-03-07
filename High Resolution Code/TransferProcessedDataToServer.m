
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
% list of files that haven't been uploaded to the server yet
rows2 = cellfun(@(x) isempty(x) || (isnumeric(x) && isnan(x)), excelfile(:, Excel.processed_data_on_server));
% list of files that have been processed far enough to be ready for processed data transfer
rows3 = cellfun(@(x) isempty(x) || (isnumeric(x) && isnan(x)), excelfile(:, Excel.trialID));
% find overlapping files that need to be transferred to the server: 
loc = find(rows1 & rows2 & ~rows3);
expList = excelfile(loc, Excel.trialID);
% offer up the list for approval
idx = listdlg('PromptString','Select experiments to copy:','ListString',expList,...
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
    for f = 1:length(folder_idx)
        subfolderName = folderList(folder_idx(f)).name;
        copyfile([originalFolder subfolderName], [destinationFolder subfolderName])
    end
    
    % move all video files into their own folder
    aviFolder = createFolder([originalFolder 'AVI_Files']);
    movefile(fullfile(originalFolder, '*.avi'), aviFolder); %TODO: set this up as a quick function for 'file1.mat' folders

    % move all data files into their own folder
    dataFolder = createFolder([originalFolder 'data']);
    movefile(fullfile(originalFolder, '*.mat'), dataFolder);
    movefile(fullfile(originalFolder, '*.slp'), dataFolder);
    movefile(fullfile(originalFolder, '*.h5'), dataFolder);
    movefile(fullfile(originalFolder, '*.csv'), dataFolder);
    movefile(fullfile(originalFolder, '*.py'), dataFolder);
    movefile(fullfile(originalFolder, '*.png'), dataFolder);

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

end






























