
% courtship step 1
clear; clc;

%% Load the data

% TODO -- update this to read from the Courtship Experiments.xlsx file
rootDir = getDataPath(5, 2, 'Select location for data');
% rootDir = uigetdir;
paths = getPathNames;
dateDir = selectFolder([rootDir, paths.courtship]);
rampName = selectFolder([rootDir paths.courtship  dateDir{1}]);
% A = strsplit(rampName{:}, '_');
% expName = strjoin(A(1:end-1), '_');
expName = rampName{:};

baseDir = [rootDir paths.courtship  dateDir{1} '\' rampName{:} '\'];

% Load pertinent data
load([baseDir expName ' parameters.mat'],'parameters'); % load exp params to workspace
tempLog = readmatrix([baseDir expName '_RampLog']); % TECA LOG full timecourse temp log

% Create a 'compile' list
compileLength = 3; %  minutes
try fragDur = parameters.fragmentlength;
catch
    fragDur = parameters.fragmentduration;
end
compile_size = floor((compileLength*60)/fragDur);
vid_rois = 0:compile_size:parameters.numFrag;
vid_rois(1) = 1;
if vid_rois(end) < parameters.numFrag
    vid_rois = [vid_rois, parameters.numFrag];
end
nVids = length(vid_rois)-1; % how many compiled videos there will be
vROI = [vid_rois(1), vid_rois(2)];
for i = 2:nVids
    vROI(i,:) = [vid_rois(i)+1, vid_rois(i+1)];
end
parameters.nVids = nVids; % add the number of videos to the parameters
parameters.vROI = vROI; % add the file index for the compiled videos

% Pull the temp & time information together: 
fileList = dir([baseDir 'file*.mat']);
dummy = [];
for i = 1:length(fileList)
    dummy(i).data = load([baseDir, 'file' num2str(i) '.mat'],'timestamp','tempLog');
end
% need to load the start temp too: 
load([baseDir, expName ' RampStartTempLog.mat'],'tempLogStart');
% Create the temp start-stop matrix for the data: 
for i = 1:parameters.numFrag
    tempLogStart(i+1,2:end) = dummy(i).data.tempLog;
end

% frames are LOCKED in the recording so we can use this to create a duration matrix
total_frames = parameters.numFrag * fragDur * parameters.FPS;
[vidNums, vidFrame, temperature, tempWork,fragNum] = deal(nan(total_frames,1));
frames = fragDur * parameters.FPS; % frames per fragment
frame = (1:total_frames)'; % total frame count

% For each of the fragments, get the information to concatenate
currFrame = 1;
for frag = 1:parameters.numFrag
    fROI = currFrame:currFrame+frames-1;
    % Frag number 
    fragNum(fROI) = frag*ones(frames,1);
   
    % Video number
    vidIdx = (find(frag>=vROI(:,1) & frag<=vROI(:,2)));
    vidNums(fROI) = vidIdx*ones(frames,1);

    % frame num in video (NOT fragment)
    if frag == vROI(vidIdx,1) %for the first frame in a chunk
        vidCount = 1;
    end
    vidFrame(fROI) = (vidCount:vidCount+frames-1)';
    vidCount = vidCount+frames;

    % temperature log
    logROI = [];
    logROI(1) = find(tempLog(:,1)==tempLogStart(frag,2));
    logROI(2) = find(tempLog(:,1)==tempLogStart(frag+1,3));
    tempCourse = tempLog(logROI(1):logROI(2),2);
    x = round(linspace(1, frames, length(tempCourse)));
    fullTempList = interp1(x,tempCourse,1:frames,'spline'); 
    temperature(fROI) = fullTempList'; 

     % temp plate work log
    workCourse = tempLog(logROI(1):logROI(2),4);
    x = round(linspace(1, frames, length(workCourse)));
    fullWorkList = interp1(x,workCourse,1:frames,'spline');   
    tempWork(fROI) = fullWorkList'; 

    % up count index
    currFrame = currFrame+frames;
end

% % Time count
time = (linspace(1, (frame(end)/parameters.FPS)/60, frame(end)))';
% figure; plot(time,temperature,'color',Color('black'))

% Data table with continuous variables:
T = table(frame, time, temperature, tempWork, vidNums, vidFrame, fragNum);

save([baseDir, expName ' alignment table'],'T','parameters','expName')

%% Find difference in time between the video saves: 

timestamps = [];
for i = 1:length(dummy)
    [h,m,s] = hms(dummy(i).data.timestamp);
    timestamps(i) = seconds(duration(h,m,s));
end

y = diff(timestamps);

fig = getfig('',false,[1123 601]); 
    scatter(1:length(y),y,50,Color('black'),'filled')
    h_line(5,'red','--',1.5)
    h_line(mean(y),'k','-',1)
    formatFig(fig)
    xlabel('Video Number')
    ylabel('Video record and save time (s)')

save_figure(fig, [baseDir, 'Video Durations'],'-png');

%% Alignment plot

kolor = Color('teal');
LW = 1.5;
r = 3;
c = 1;
fig = getfig('',1);
subplot(r,c,1)
    plot(T.time, T.temperature,'color', kolor,'linewidth', LW)
    ylabel('temp (\circC)')
subplot(r,c,2)
    plot(T.time, T.vidNums,'color', kolor,'linewidth', LW)
    ylabel('Compiled Video (#)')
subplot(r,c,3)
    plot(T.time, T.fragNum,'color', kolor,'linewidth', LW)
    ylabel('Fragement (#)')
    xlabel('Time (min)')

formatFig(fig,false,[r,c]);
subplot(r,c,1)
set(gca, 'xcolor', 'none')
subplot(r,c,2)
set(gca, 'xcolor', 'none')

% for i = 1:3
%     subplot(r,c,i)
%     xlim([30,37])
% end

save_figure(fig, [baseDir, 'Video Alignments'],'-png');
















