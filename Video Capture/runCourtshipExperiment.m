
clear; clc; imaqreset; parameters = [];

% PARAMETERS THAT MAY CHANGE EACH TIME
expName = 'Berlin_courtship_F_LRR_caviar_ramp3';% Experiment name
parameters.protocol = 'courtship_F_LRR_25-17'; % Temperature Protocol
parameters.experimenter = 'Becca';  % Experimenter
parameters.ArenaA.genotype = 'Berlin'; % Genotype
parameters.ArenaA.well_1 = 'Caviar'; % Well 1
parameters.ArenaA.well_2 = 'Empty'; % Well 2
parameters.ArenaA.well_3 = 'Empty'; % Well 3
parameters.ArenaA.well_4 = 'Empty';% Well 4
parameters.ArenaA.sex = 'mixed'; % Sex of flies

% Set up the experiment parameters
dateStr = char(datetime('now','Format','MM.dd.yyyy'));
baseFolder = ['E:\Evyn\Courtship Tracking\' dateStr '\'];

if ~exist(baseFolder)
    mkdir(baseFolder)
end
cd(baseFolder) % set current path to video folder

start_pause = 5; % min delay before recording should start
totalLength = 16*4; % 16 min pre, 16 down, 16 up, 16 post
hz = 30; % camera FPS 
trial = 5; % fragment recording durations (sec)
numSections = 1;
ITI = 48;

nsamples = ceil((totalLength*60)/trial);

% Save all the selected parameters 
parameters = [];
parameters.nSections = numSections;
parameters.ITI = 5;
parameters.FPS = hz;
parameters.recordingduration = totalLength;
parameters.fragmentduration = trial;
parameters.numFrag = nsamples; % number of video fragments
parameters.day_night = 'D'; %'N' 
parameters.date = dateStr;
parameters.expID = expName;
parameters.videoName = expName;
parameters.ArenaA.starved_hours = 0;


% save parameter file
if isfile([baseFolder expName 'dataMat.mat']) %prevent overwriting the file
    warndlg('Check the experiment name to not overwrite files')
    return
end
save([expName ' parameters.mat'],"parameters")

% Make a temperature log file:
tempLogPath = [parameters.videoName '_RampLog.csv'];
fclose(fopen(tempLogPath, 'w+'));

disp('Parameters saved')

%% Initialize camera
vid = videoinput('pointgrey', 1, 'F7_Raw8_2048x2048_Mode0');
% partial_ROI = [0 0 2048 2048];
src = getselectedsource(vid);
[src, vid] = initialize_CourtshipCamera(src,vid,hz);
% % camera parameters
% src.Brightness = 29;
% src.Exposure  = 1.5648;
% src.FrameRate = hz;
% src.Gain = 1.752491; 
% src.Gamma = 1.5338; 
% src.Shutter = 11.6188; 
% vid.FramesPerTrigger = inf;
% vid.ROIPosition = partial_ROI;
% disp('Purple BACK camera initialized')

h = preview(vid);
    
disp('Ready to start the experiment')
    
%% Run the experiment   
disp('Start temp protocol')
tempLogStart = [];

% start delay on the experiment
pause(start_pause*60)
% pause(10)
    
for n = 1:numSections
    % save start temp log information    
    tempLogStart(n,:) = [n, logTempNow(tempLogPath)];
    disp(tempLogStart(n,:))
    save([expName ' RampStartTempLog.mat'],'tempLogStart')

    % acquire videos
    nFrames = trial*hz*nsamples;
    get_samples_v3(nFrames,trial*hz,expName,tempLogPath,hz);
    disp('recording finished')
    
    % pause for ITI (if more than 1 section)
    if numSections>1 && n<numSections
        pause(ITI*60) % pause for the duration of each video
    elseif n==numSections
        pause(start_pause*60)
    end
end

disp(tempLogStart(n,:))

save([expName ' RampStartTempLog.mat'],'tempLogStart')

fprintf(['\n' 'Experiment Done' '\n']) 




































