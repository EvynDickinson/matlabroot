

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
fileOptions = dir([structFolder 'Berlin hold*.mat']);
nExp = length(fileOptions);

trialTables = cell(nExp, 1); % initialize cell struct for all experiment tables

for ii = 1:(nExp)
    tic 
    filename = [structFolder fileOptions(ii).name];

    temp = load([structFolder fileOptions(ii).name]);
   
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
    fprintf('\n Finished %s\n', fileOptions(ii).name)

    % save table to main data folder?
    save([structFolder, fname ' table.mat'], 'T', '-v7.3');
    toc
end

T = vertcat(trialTables{:});

figure; hold on
subplot(2,1,1); plot(T.Temperature)
subplot(2,1,2); plot(smooth(T.InnerFoodQuad, 360, 'moving')) % 2 min moving filter



%% Save the concatenated file: 
save([structFolder 'Berlin Caviar dynamic temps flat_extracted.mat'], 'T', '-v7.3');


%% 





    
    
    
