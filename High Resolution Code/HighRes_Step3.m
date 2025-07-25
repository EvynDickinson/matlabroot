

%% Load tracking points for courtship
clear; clc;
baseFolder = getDataPath(6,0);

% Find files that can be run
[excelfile, Excel, xlFile] = load_HighResExperiments;

% Find trials that have nothing in their Proofed Column & have been tracked: 
loc = cellfun(@isnan,excelfile(2:end,Excel.tracked));
loc = ~loc; 
rownums = find(loc)+1; 
eligible_files = excelfile([false;loc],[Excel.date, Excel.expID, Excel.ramp, Excel.proofed]);
FileNames = format_eligible_files(eligible_files);%,filler_string);

fileIdx = listdlg('ListString', FileNames,'ListSize',[350,450],'promptstring', 'Select data to process');
if isempty(fileIdx)
    disp('No trials selected')
    return
end
% pull the list of dates and arenas to be 
List.date = eligible_files(fileIdx,1);
List.expID = eligible_files(fileIdx,2); 
ntrials = length(fileIdx);

initial_var = who;
initial_var{end+1} = 'trial';
initial_var{end+1} = 'initial_var';

%% Load and check for proofing needs
for trial = 1:ntrials
    baseDir = [baseFolder, List.date{trial} '\', List.expID{trial} '\'];
    fileList = dir([baseDir '*alignment table.mat']);
    load([baseDir, fileList(1).name]) % load the parameters and temp table
    nvids = parameters.nVids; % number of videos
    nBP = 5; % num of body parts

    % load the video tracks
    fileList = dir([baseDir '*.h5']);
    data = [];
    % TODO: quick check to make sure this matches the expected number of videos...
    for vid = 1:nvids
        filePath = [baseDir, 'compiled_video_' num2str(vid) '.avi.predictions.slp.h5'];
        data(vid).occupancy_matrix = h5read(filePath,'/track_occupancy');
        % data(vid).tracks = h5read(filePath,'/tracks');
        % data(vid).node_names = h5read(filePath, '/node_names');
        % data(vid).track_names = h5read(filePath, '/track_names');
        dataIn = true; % log successful data load
    end

    % Quality control
    ntracks = [];
    for vid = 1:parameters.nVids
        ntracks(vid) = size(data(vid).occupancy_matrix,1);
    end

%%%%%%%%%%%%%% MAYBE HERE CLEAR EXTRA TRACKS FOR FROM 15C HR EXPS %%%%%%%%%%%%%%%%%%
    % Open the file for writing
    outputFile = fullfile(baseDir, 'proof_videos.txt');
    fid = fopen(outputFile, 'w'); % 'w' for writing (overwrites if the file exists)
    if fid == -1
        error('Could not open file for writing.');
    end
    % Identify videos to proof
    vidloc = find(~(ntracks == 2));
    fprintf('Proof the following videos: \n');
    disp(baseDir);
    fprintf(fid, 'Proof the following videos: \n');
    fprintf(fid, '%s\n', baseDir);
    % Loop through the videos and output to file and console
    for i = vidloc
        outputStr = sprintf('video %d | current tracks %d\n', i, ntracks(i));
        fprintf('%s', outputStr); % Display in command window
        fprintf(fid, '%s', outputStr); % Write to file
    end
    % Close the file
    fclose(fid);
    disp(['Output written to: ' outputFile]);

    % Copy the original files to the backup folder: 
    response = questdlg('Copy files to backup folder?');
    if strcmpi(response, 'Yes')
        a = [baseDir 'Tracking backup'];
        filelist = dir([a '*']);
        foldercount = length(filelist)+1;
        enddir = [a ' ' num2str(foldercount)];
        if ~exist(enddir,"dir")
            mkdir(enddir)
        end
        for i = vidloc
            targetfile = [baseDir 'compiled_video_' num2str(i) '.avi.predictions.slp'];
            copyfile(targetfile,enddir) % slp file
            copyfile([targetfile '.h5'],enddir) % h5 file
        end
    end 

    % write 'R' for ready into the excel sheet: 
    if ~isExcelFileOpen(xlFile,true)
        writecell({'R'},xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.proofed) num2str(fileIdx(trial)+1)]);
    else 
        disp('Couldn''t write to excel sheet for:')
        disp([List.date{trial} ' ' List.expID{trial}])
    end
disp('  ')
end

disp('All finished')

































