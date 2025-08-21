
clear; clc
imaqreset

%% Set up the experiment parameters
parameters = [];
dateStr = char(datetime('now','Format','MM.dd.yyyy'));
[drive, parameters.protocol] = availableExperimentDrives; % get drive location for saving
if isnan(drive); return; end
baseFolder = createFolder([drive '\Courtship Videos\' dateStr '\']);
cd(baseFolder) % set current path to video folder
genotype = getGenotype; % pull genotype information
[food, wellLoc] = getFoodType(1);% pull food well option

switch parameters.protocol
    case 'courtship_F_LRR_25-17'
        totalLength = 64; % 16 min pre, 16 down, 16 up, 16 post
        rampNum = questdlg('Ramp number?','','1','2','3','1'); %todo -- make this more dynamic and search for # in folder
        expName = [genotype '_' parameters.protocol '_' food '_ramp' rampNum];

    case 'high_res_LTS_35-15'
        totalLength = 485; %
        expName = [genotype '_' parameters.protocol '_' food];

    case {'Hold15C', 'Hold20C', 'Hold25C', 'Hold30C', 'Hold35C'}
        totalLength = 485; %
        expName = [genotype '_' parameters.protocol '_' food];

end

% parameters.protocol = 'Hold35C';
% expName = [genotype '_' parameters.protocol '_' food];

% auto-select and fill experiment information: 
% expName = 'Berlin_courtship_F_LRR_caviar_ramp1';
start_pause = 5; % min delay before recording should start
hz = 30; % camera FPS 
trial = 5; % fragment recording durations (sec)
numSections = 1;
ITI = 48;
nsamples = ceil((totalLength*60)/trial);

% Save all the selected parameters 
parameters.nSections = numSections;
parameters.ITI = 5;
parameters.FPS = hz;
parameters.recordingduration = totalLength;
parameters.fragmentduration = trial;
% parameters.protocol = 'courtship_F_LRR_25-17'; %TODO dynamic
parameters.numFrag = nsamples; % number of video fragments
parameters.experimenter = questdlg('Select experimenter:','','Becca', 'Evyn', 'Cancel', 'Becca');
parameters.day_night = 'D'; 
parameters.date = dateStr;

% Experiment:
parameters.expID = expName;
parameters.videoName = expName;
parameters.ArenaA.genotype = genotype;
parameters.ArenaA.sex = 'mixed';
parameters.ArenaA.starved_hours = 0;
for i = 1:4
    parameters.ArenaA.(['well_' num2str(i)]) = 'Empty';
end
parameters.ArenaA.(['well_' wellLoc]) = food;

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

% for n = 1:numSections
n = 1;
% save start temp log information    
tempLogStart(n,:) = [n, logTempNow(tempLogPath)];
disp(tempLogStart(n,:))
save([expName ' RampStartTempLog.mat'],'tempLogStart')

% acquire videos
nFrames = trial*hz*nsamples;
get_samples_v3(nFrames,trial*hz,expName,tempLogPath,hz);
disp('recording finished')

% % pause for ITI (if more than 1 section)
% if numSections>1 && n<numSections
%     pause(ITI*60) % pause for the duration of each video
% elseif n==numSections
pause(start_pause*60)
%     end
% end

disp(tempLogStart(n,:))

save([expName ' RampStartTempLog.mat'],'tempLogStart')

fprintf(['\n' 'Experiment Done' '\n']) 






%%
% % USE THIS TO PULL OUT THE TIME AND TEMP DATA FOR AN EXPERIMENT
% basepath = 'E:\Evyn\Courtship Tracking\02.10.2025\Dummy_high_res_LTS_35-15_Caviar_1\';
% 
% timestamps = NaT(5820,1);
% tempLogs = NaN(5820,5);
% 
% for i = 1:5820
%     dummy = load([basepath 'file' num2str(i)],'tempLog','timestamp');
%     timestamps(i,1) = dummy.timestamp;
%     tempLogs(i,:) = dummy.tempLog;
%     disp([num2str(i) ' / 5820'])
% end
% 
% save('LTS test timing data.mat')
 



























