

%% Find pathways to extract only the specific information needed from existing trial data for model analysis

% From grouped data structure: 
% (one value that needs to be expanded to fit all the data points in time)
% -----------------------------------------------------------
% fly id --> data(temp protocol type).T.ExperimentID(trials)
% temp protocol --> data(temp protocol type).T.TempProtocol(trials)
% arena --> data(temp protocol type).T.Arena(trials)
% date --> data(temp protocol type).T.Date(trials)


% Data that already exists across all of time
% -----------------------------------------------------------
% time in experiment
% temperature
% flies in outer ring
% flies in inner food quad
% flies sleeping
% fly speed

%%  Load in selected data? Use datastore to find all the options and select from there?

warning off 
format shortG
clear; close all; clc
paths = getPathNames; % get the appropriate file path names
baseFolder = getDataPath(3,0,'Select where you want to find the grouped data structures');
structFolder = [baseFolder 'PID Modeling\'];
fileOptions = dir([structFolder '*.mat']);
fileNames = {fileOptions(:).name};
% select the trials that you want to load: 
fileIdx = listdlg('ListString', fileNames, 'SelectionMode', 'multiple','ListSize',[500,700]);
fileNames = fileNames(fileIdx);
nExp = length(fileIdx);

trialTables = cell(nExp, 1); % initialize cell struct for all experiment tables

for ii = 1:(nExp)
    tic 
    filename = [structFolder fileNames{ii}];
    temp = load([structFolder fileNames{ii}]);
   
    % find default matrix sizes
    ntrials = size(temp.all.speed.all,2);
    npoints = size(temp.all.speed.all,1);

    % create a unique trial ID
    [~, fname] = fileparts(filename);
    TrialID = compose('%s_t%02d', fname, (1:ntrials)); % make unique trial IDs
    % TrialID = repmat(string(TrialID), [npoints, 1]);
    TrialID = repmat(TrialID, [npoints, 1]);
    TrialID = TrialID(:);
    
    % extract timeseries data from the different experiment types: 
    InnerFoodQuad = temp.all.innerquad.food.all;
    EscapeRing = temp.all.ring.all;
    Sleep = temp.all.sleep.all;
    Speed = temp.all.speed.all;

    % data that is by trial and not by time: 
    a = temp.all.T.Date';
    Date = repmat(a, [npoints, 1]);
    a = temp.all.T.TempProtocol';
    TempProtocol = repmat(a, [npoints, 1]);
    % a = temp.all.T.foodName';
    % FoodName = repmat(a, [npoints, 1]);
    % a = temp.all.T.Genotype';
    % Genotype = repmat(a, [npoints, 1]);

    % single data trial information
    Time = temp.all.time;
    Temperature = temp.all.real_temp;

    % reorganize the data in a time-series (concatenate)
    InnerFoodQuad = InnerFoodQuad(:);
    EscapeRing = EscapeRing(:);
    Sleep = Sleep(:);
    Speed = Speed(:);
    Time = repmat(Time,[ntrials,1]);
    Temperature = repmat(Temperature, [ntrials, 1]);
    Date = Date(:);
    TempProtocol = TempProtocol(:);
    
    % Combine these into a giant table: 
    T = table(TrialID, Date, TempProtocol, Time, Temperature, InnerFoodQuad, EscapeRing, Sleep, Speed);
    trialTables{ii} = T;
    fprintf('\n Finished %s\n', fileNames{ii})

    % % save table to main data folder?
    % save([structFolder, fname ' table.mat'], 'T', '-v7.3');
    % % save(['D:\PID Modeling Data\', fname ' table.mat'], 'T', '-v7.3');

    % save as a flattened parquet file: 
    parquetwrite([structFolder, fname ' table'], T);
    
    toc
end

T = vertcat(trialTables{:});

figure; hold on
subplot(2,1,1); plot(T.Temperature)
subplot(2,1,2); plot(smooth(T.InnerFoodQuad, 360, 'moving')) % 2 min moving filter


%% Save the concatenated file: 


% save as a giant file: 
dataFile = [structFolder 'Berlin Caviar all data.parquet'];
tic
parquetwrite(dataFile, T);
toc

% 
% % saving created data structure -- 
% save(['D:\PID Modeling Data\', 'Berlin Caviar dynamic temps flat_extracted.mat'], 'T', '-v7.3');
% % structFolder

%% Load static data or dynamic or both into the workspace

% % point to the data ....
% dynamic_path = 'D:\Evyn Lab Data\Data structures\PID Modeling Data\Berlin Caviar dynamic temps flat_extracted.mat';
% static_path = "D:\Evyn Lab Data\Data structures\PID Modeling Data\Berlin Caviar static temps flat_extracted.mat";
% % load the data ...
% temp = load(dynamic_path);
% Td = temp.T;
% temp = load(static_path);
% Ts = temp.T;
% % combine the data into a giant structure
% T = [Td; Ts];


%% Load from parquet style flattened data: 
clear
close all
clc

% look for available parquet data files to load and compile: 
paths = getPathNames; % get the appropriate file path names
baseFolder = getDataPath(3,0,'Select where you want to find the grouped data structures');
structFolder = [baseFolder 'PID Modeling\'];
% select the trials that you want to load: 
fileOptions = dir([structFolder '*.parquet']);
nExp = length(fileOptions);

trialTables = cell(nExp, 1); % initialize cell struct for all experiment tables

for ii = 1:nExp
    dataFile = [structFolder, fileOptions(ii).name];
    tic
    trialTables{ii} = parquetread(dataFile);
    toc
    fprintf('Loaded %s', fileOptions(ii).name);
end

% compiile the data files into a large table: 
T = vertcat(trialTables{:});

% save compiled table ... 
fileSaveName = 'Berlin caviar all';
parquetwrite([structFolder, fileSaveName], T);

% parquetwrite("K:\DATA\Data structures\PID Modeling\Berlin caviar all", T);



% % read only specific columns -- very fast, skips the rest
% T = parquetread('data.parquet', 'SelectedVariableNames', {'TrialID', 'Temperature', 'EscapeRing'})


%%  convert existing .mat table files to parquet style 

warning off 
format shortG
clear; close all; clc
paths = getPathNames; % get the appropriate file path names
baseFolder = getDataPath(3,0,'Select where you want to find the grouped data structures');
structFolder = [baseFolder 'PID Modeling\'];
fileOptions = dir([structFolder '*table.mat']);
fileNames = {fileOptions(:).name};
% select the trials that you want to load: 
fileIdx = listdlg('ListString', fileNames, 'SelectionMode', 'multiple','ListSize',[500,700]);
fileNames = fileNames(fileIdx);
nExp = length(fileIdx);

for ii = 1:(nExp)
    tic 
    filename = [structFolder fileNames{ii}];
    load(filename, 'T'); % load the 'T' table
   
    % save as a flattened parquet file: 
    parquetwrite(filename(1:end-4), T);
    disp(fileNames{ii})
    toc
end






    
    
