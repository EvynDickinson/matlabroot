
% Single fly manual tracking through time

%% General Setup

clear; clc

% Load file
load('G:\My Drive\Jeanne Lab\DATA\09.19.2022\Arena C\caviar_recovery_rampC timecourse data.mat');

baseFolder = getCloudPath;
figDir = [baseFolder 'Manual Tracking\'];
expDate = expData.parameters.date;
expArena = data.name;


% Pull ramp start times
TP = getTempTurnPoints(expData.parameters.protocol);

% 'Buffer' before/after ramps
buffer = 5; %minutes

% Ramps 1-4
vidROI = [];
vidROI(:,1) = TP.down(:,1)-(buffer*3*60);
vidROI(:,2) = TP.up(:,2)+(buffer*3*60);

% Do the regions make sense? 
fig = figure; plot(data.T.time,data.T.temperature,'color', 'w')
    vline(data.T.time(vidROI(:,1)),'b')
    vline(data.T.time(vidROI(:,2)),'r')
    xlabel('time (min)')
    ylabel('temperature (\circC)')
    formatFig(fig,true);

save_figure(fig,[figDir, 'time ROI selection'],'-png');


% Find video list that aligns for each ramp
Vid_starts = data.T.vidNums(vidROI(:,1));
Vid_ends = data.T.vidNums(vidROI(:,2));


%% Get video cropping control

movieInfo = VideoReader([baseFolder expDate '/' expName '_1.avi']); %read in video
demoImg = rgb2gray(read(movieInfo,1));

fig = figure;
    imshow(demoImg); hold on
    x = data.centre(1);
    y = data.centre(2);
    r = data.r;
    scatter(x,y,50,'r','filled')
    x_lims = [x-(r+50),x+(r+50)];
    y_lims = [y-(r+50),y+50+r];
    xlim(x_lims); ylim(y_lims)

% Make sure the image zoom fits the whole arena!


%% Set parameters and variables

nTrack = 3; % number of flies to track
nPast = 3; % number of previous tracked points to display

% structure to load position data for appropriate frames...
MT = nan([size(data.T,1),2]);







rames = 3;
for arena = 1:4
    arenaData(arena).nflies = excelfile{XLrow(arena),Excel.numflies};
    nflies(arena) = arenaData(arena).nflies;
end
movieInfo = VideoReader([baseFolder folder '/' expName '_1.avi']); %read in video
ii = randi(size(data(1).occupancy_matrix,2),[3,1]); %random selection of frames to count fly number
demoImg = rgb2gray(read(movieInfo,1));
if any(isnan(nflies)) 
    nflies = [];
    % manual count of flies
    fprintf('\nCount the number of flies in the picture by clicking them\n then hit ENTER\n')
    T = true;
   while T % get the number of flies in the arena
        for jj = 1:nframes
            demoImg = rgb2gray(read(movieInfo,ii(jj)));
            PT(jj).frame = readPoints(demoImg);
            for arena = 1:4
                centre = arenaData(arena).centre;
                X = PT(jj).frame(1,:);
                Y = PT(jj).frame(2,:);
                % find points within each arena sphere
                loc = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5)<=r;
                nflies(arena,jj) = sum(loc);
            end
        end






demoImg = rgb2gray(read(movieInfo,1));
radii = 165; %well surround regions


























