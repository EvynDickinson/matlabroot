

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
% select the trials that you want to load: 
fileOptions = dir([structFolder '*.mat']);

for ii = 1:length(fileOptions)
    filename = [structFolder fileOptions(ii).name];

    % temp = load([structFolder fileOptions(ii).name]);
    % 
    % % extract timeseries data from the different experiment types: 
    % InnerFoodQuad = temp.all.innerquad.food.all;
    % EscapeRing = temp.all.ring.all;
    % Sleep = temp.all.sleep.all;
    % Speed = temp.all.speed.all;
    % 
    % % single data trial information
    % Time = temp.all.time;
    % Temperature = temp.all.real_temp;
    % 
    % % data that is by trial and not by time: 
    % npoints = length(Time);
    % a = temp.all.T.Date';
    % Date = repmat(a, [npoints, 1]);
    % a = temp.all.T.TempProtocol';
    % TempProtocol = repmat(a, [npoints, 1]);
    % a = temp.all.T.foodName';
    % FoodName = repmat(a, [npoints, 1]);
    % a = temp.all.T.Genotype';
    % Genotype = repmat(a, [npoints, 1]);
    % 
    % % reorganize the data in a time-series? or not yet, since if we want to
    % % do a temperature signal transformation we want to know the time
    % % series for each so they can operate independently? (e.g. not pull
    % % past time history from a different trial into the current one because
    % % they are concatenated) 
    % 
    % % Combine these into a giant table: 
    % T = table(Date, TempProtocol, Genotype, FoodName, Time, Temperature, InnerFoodQuad, EscapeRing, Sleep, Speed);

end

%% STEP 2: build the raw datastore
ds_raw = fileDatastore([structFolder '*.mat'], ...
    'ReadFcn', @extractTrialData, ...
    'FileExtensions', '.mat');

% extract once and save -- don't re-run extraction every time
T_raw = readall(ds_raw);
T_raw = vertcat(T_raw{:});   % fileDatastore returns a cell, so concatenate
save([structFolder 'flat_extracted.mat'], 'T_raw', '-v7.3');


%% 

function T = extractTrialData(filename)
    temp = load(filename);

    % --- shared time series (n points) ---
    Time = temp.all.time(:);
    Temperature = temp.all.real_temp(:);
    npoints = length(Time);

    % --- trial-varying data (n x m matrices) ---
    InnerFoodQuad = temp.all.innerquad.food.all;   % n x m
    EscapeRing = temp.all.ring.all;              % n x m
    Sleep = temp.all.sleep.all;             % n x m
    Speed = temp.all.speed.all;             % n x m
    ntrials = size(EscapeRing, 2);            % m varies per file

    % --- trial-level metadata (one value per trial, m entries) ---
    Date = temp.all.T.Date(:);         % m x 1
    TempProtocol = temp.all.T.TempProtocol(:); % m x 1
    FoodName  = temp.all.T.foodName(:);     % m x 1
    Genotype  = temp.all.T.Genotype(:);     % m x 1

    % --- file identity ---
    [~, fname] = fileparts(filename);

    % --- build one table per trial, then concatenate ---
    trialTables = cell(ntrials, 1);
    for m = 1:ntrials
        TrialID  = repmat(sprintf('%s_t%02d', fname, m), [npoints, 1]);
        TrialTimeIndex = (1:npoints)';

        % expand scalar trial metadata to match time series length
        trialTables{m} = table(...
            repmat(string(TrialID), [npoints, 1]), ...
            TrialTimeIndex, ...
            repmat(Date(m), [npoints, 1]), ...
            repmat(TempProtocol(m), [npoints, 1]), ...
            repmat(Genotype(m), [npoints, 1]), ...
            repmat(FoodName(m), [npoints, 1]), ...
            Time, ...
            Temperature, ...
            InnerFoodQuad(:, m), ...
            EscapeRing(:, m), ...
            Sleep(:, m), ...
            Speed(:, m), ...
            'VariableNames', {'TrialID', 'TrialTimeIndex', 'Date', ...
                              'TempProtocol', 'Genotype', 'FoodName', ...
                              'Time', 'Temperature', 'InnerFoodQuad', ...
                              'EscapeRing', 'Sleep', 'Speed'});
    end

    T = vertcat(trialTables{:});
end




    
    
    
