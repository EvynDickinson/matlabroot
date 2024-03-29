

%% Pixel distance between wells in the new and old camera / lens setups:

% open image and measure distance between points...
clear

% Select video
path = uigetdir;
vid_opt = dir([path '\*.avi']);
vid_opt = {vid_opt(:).name};
idx = listdlg('ListString', vid_opt,'ListSize', [250, 400]);
video_path = [path '\' vid_opt{idx}];

% Load frame 1
movieInfo = VideoReader(video_path); 
demoImg = rgb2gray(read(movieInfo,1));

% Click the inside points of the top right image wells, repeat 3x for
% 'accuracy'
% Zoom in?
[W12, W3, W6, W9] = deal([]);
for i = 1:4
    wellcenters = readPoints(demoImg,4); % get locations of the wells from the image
    W12(i,:) = wellcenters(:,1)';
    W3(i,:) = wellcenters(:,2)';
    W6(i,:) = wellcenters(:,3)';
    W9(i,:) = wellcenters(:,4)';
end

% Find the average pixel distance between the 4 wells
Vertical_avg = sqrt((mean(W12(:,1))-mean(W6(:,1)))^2 + (mean(W12(:,2))-mean(W6(:,2)))^2);
Horizontal_avg = sqrt((mean(W3(:,1))-mean(W9(:,1)))^2 + (mean(W3(:,2))-mean(W9(:,2)))^2);

D = [];
D(1) = sqrt((mean(W12(:,1))-mean(W3(:,1)))^2 + (mean(W12(:,2))-mean(W3(:,2)))^2);
D(2) = sqrt((mean(W3(:,1))-mean(W6(:,1)))^2 + (mean(W3(:,2))-mean(W6(:,2)))^2);
D(3) = sqrt((mean(W6(:,1))-mean(W9(:,1)))^2 + (mean(W6(:,2))-mean(W9(:,2)))^2);
D(4) = sqrt((mean(W9(:,1))-mean(W12(:,1)))^2 + (mean(W9(:,2))-mean(W12(:,2)))^2);


disp(vid_opt{idx})
disp(['Mean inter-well distance: ' num2str(mean(D)) ' +/- ' num2str(std(D))  ' pixels'])
disp('   ')

% blue: 291.5144 2.6718
% purp: 292.7555  2.2574
% older: 290.4455 0.45493

% Plot graph of the inter-well distance across setups for Jaime

fig = getfig('',1,[422 680]); hold on
kolor1 = Color('dodgerblue');
errorbar(1, 291.5144, 2.6718,'Color', kolor1, 'linewidth', 1.5,'Marker','o', 'MarkerFaceColor',kolor1);
kolor2 = Color('mediumpurple');
errorbar(2, 292.7555,2.2574,'Color', kolor2, 'linewidth', 1.5,'Marker','o', 'MarkerFaceColor',kolor2);
kolor3 = Color('white');
errorbar(3, 290.4455,0.45493,'Color', kolor3, 'linewidth', 1.5,'Marker','o', 'MarkerFaceColor',kolor3);
xlim([0,4]); ylim([0,300])
ylabel('Distance (pixels)')
formatFig(fig,true);
set(gca, 'xcolor', 'k')

%% Video save timing

% open image and measure distance between points...
clear

% Get blue camera data:
disp('Select blue video data mat')
[file, path] = uigetfile; 
BM = load([path, file]);
% Pull video information: 
dummy = strsplit(file, ' ');
Blue_root = dummy{1};
movieInfo = VideoReader([path, Blue_root '_1.avi']); 
disp(movieInfo)

% Get purple camera data:
disp('Select purple video data mat')
file = uigetfile([path '\*dataMat.mat']); 
PM = load([path, file]);
% Pull video information: 
dummy = strsplit(file, ' ');
Purple_root = dummy{1};

% Load temperature log
disp('select temperature log')
file = uigetfile([path '\*RampLog.csv']);
tempLog = readmatrix([path, file]);


% Check for differences in start/stop of purple & blue cameras

% TIMEPOINT BASED:
startDiff = PM.tempLogStart(:,3)-BM.tempLogStart(:,3);
stopDiff = PM.tempLogEnd(:,3)-BM.tempLogEnd(:,3);
T = [(1:BM.parameters.numVids)', startDiff, stopDiff];
% Display information about loss/etc
disp(['Frame rate: ' num2str(BM.parameters.FPS)])
disp(['Vid length: ' num2str(BM.parameters.video_length)])
disp(['Number of videos: ' num2str(BM.parameters.numVids)])
table(T);
disp(T)


% INDEX BASED:
startDiff = PM.tempLogStart(:,2)-BM.tempLogStart(:,2);
stopDiff = PM.tempLogEnd(:,2)-BM.tempLogEnd(:,2);
T = [(1:BM.parameters.numVids)', startDiff, stopDiff];
% Display information about loss/etc
% disp(['Frame rate: ' num2str(BM.parameters.FPS)])
table(T);
disp(T)


% What are the expected times between video samples? 
% video length / 15, yeah?

video_durations_B = ((diff(BM.tempLogStart(:,2)))*15)./60;
video_durations_P = ((diff(PM.tempLogStart(:,2)))*15)./60;
figure; hold on
x = shuffle_data(linspace(1,2,length(video_durations_B)));
scatter(x,video_durations_B,35, 'b', 'filled')
x = shuffle_data(linspace(2,3,length(video_durations_P)));
scatter(x,video_durations_P,35, Color('purple'), 'filled')


%%

% How many frames are in each of the videos? ie are they saving the same
% duration?
[FramesPerSec, FrameCount] = deal([]);
for i = 1:BM.parameters.numVids

    movieInfo = VideoReader([path, Blue_root '_' num2str(i) '.avi']); 
    FrameCount(i,1) = movieInfo.NumFrames;
    FramesPerSec(i,1) = movieInfo.FrameRate;

    movieInfo = VideoReader([path, Purple_root '_' num2str(i) '.avi']); 
    FrameCount(i,2) = movieInfo.NumFrames;
    FramesPerSec(i,2) = movieInfo.FrameRate;

end

FrameNum = median(FrameCount(:,1));
FPS = median(FramesPerSec(:,1));
% 
% (FrameNum./FPS)./60
% 
totalTime = (FrameCount/FPS)/60;

fig = figure;
hold on
scatter(1.25*ones(1,length(totalTime)),totalTime(:,1),35,'b','filled')
plot([1.15,1.35],[BM.parameters.video_length, BM.parameters.video_length],'color', 'w','linewidth', 2)
scatter(1.75*ones(1,length(totalTime)),totalTime(:,1),35,Color('purple'),'filled')
plot([1.65,1.85],[BM.parameters.video_length, PM.parameters.video_length],'color', 'w','linewidth', 2)
xlim([1,2])
ylim([9,12])
formatFig(fig,true);
set(gca, 'xcolor', 'k')



% disp(['Video recording FPS: ' num2str(median(FramesPerSec))])
% disp(['Assigned FPS: ' num2str(BM.parameters.FPS(:)) ' ' num2str(PM.parameters.FPS)])
% 
% disp(['Video recording length: ' num2str()])
% 
% std(FrameCount(:,2))
% 
% FrameCount(:,2)-FrameCount(:,1)
% 
% 
% BM.parameters.FPS
































