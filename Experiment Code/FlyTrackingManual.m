% 
% % Single fly manual tracking through time
% 
% %% General Setup
% 
% clear; clc
% 
% % Load file
% load('G:\My Drive\Jeanne Lab\DATA\09.19.2022\Arena C\caviar_recovery_rampC timecourse data.mat');
% 
% baseFolder = getCloudPath;
% figDir = [baseFolder 'Manual Tracking\'];
% expDate = expData.parameters.date;
% expArena = data.name;
% 
% 
% % Pull ramp start times
% TP = getTempTurnPoints(expData.parameters.protocol);
% 
% % 'Buffer' before/after ramps
% buffer = 20; %minutes
% 
% % Ramps 1-4
% vidROI = [];
% vidROI(:,1) = TP.down(:,1)-(buffer*3*60);
% vidROI(:,2) = TP.up(:,2)+(buffer*3*60);
% 
% % Do the regions make sense? 
% fig = figure; plot(data.T.time,data.T.temperature,'color', 'w')
%     vline(data.T.time(vidROI(:,1)),'b')
%     vline(data.T.time(vidROI(:,2)),'r')
%     xlabel('time (min)')
%     ylabel('temperature (\circC)')
%     formatFig(fig,true);
% 
% save_figure(fig,[figDir, 'time ROI selection'],'-png');
% 
% 
% % Find video list that aligns for each ramp
% Vid_starts = data.T.vidNums(vidROI(:,1));
% Vid_ends = data.T.vidNums(vidROI(:,2));
% 
% 
% %% Get video cropping control
% 
% movieInfo = VideoReader([baseFolder expDate '/' expName '_1.avi']); %read in video
% demoImg = rgb2gray(read(movieInfo,1));
% 
% fig = figure;
%     imshow(demoImg); hold on
%     x = data.centre(1);
%     y = data.centre(2);
%     r = data.r;
%     scatter(x,y,50,'r','filled')
%     x_lims = [x-(r+50),x+(r+50)];
%     y_lims = [y-(r+50),y+50+r];
%     xlim(x_lims); ylim(y_lims)
% 
% % Make sure the image zoom fits the whole arena!
% 
% 
% %% Set parameters and variables
% 
% nTrack = 3; % number of flies to track
% nPast = 3; % number of previous tracked points to display
% 
% % structure to load position data for appropriate frames...
% MT = nan([size(data.T,1),2]);
% 
% 
% 
% 
% 
% 
% 
% rames = 3;
% for arena = 1:4
%     arenaData(arena).nflies = excelfile{XLrow(arena),Excel.numflies};
%     nflies(arena) = arenaData(arena).nflies;
% end
% movieInfo = VideoReader([baseFolder folder '/' expName '_1.avi']); %read in video
% ii = randi(size(data(1).occupancy_matrix,2),[3,1]); %random selection of frames to count fly number
% demoImg = rgb2gray(read(movieInfo,1));
% if any(isnan(nflies)) 
%     nflies = [];
%     % manual count of flies
%     fprintf('\nCount the number of flies in the picture by clicking them\n then hit ENTER\n')
%     T = true;
%    while T % get the number of flies in the arena
%         for jj = 1:nframes
%             demoImg = rgb2gray(read(movieInfo,ii(jj)));
%             PT(jj).frame = readPoints(demoImg);
%             for arena = 1:4
%                 centre = arenaData(arena).centre;
%                 X = PT(jj).frame(1,:);
%                 Y = PT(jj).frame(2,:);
%                 % find points within each arena sphere
%                 loc = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5)<=r;
%                 nflies(arena,jj) = sum(loc);
%             end
%         end
% 
% 
% 
% 
% 
% 
% demoImg = rgb2gray(read(movieInfo,1));
% radii = 165; %well surround regions
% 
% 
% 
% 
% 
% %%
% fig = gcf;
% frame = 32761;
% save_figure(fig,['G:\My Drive\Jeanne Lab\DATA\Manual Tracking\tracking progress ' num2str(frame)],'-png');
% 
% %%
%%
load('G:\My Drive\Jeanne Lab\DATA\Manual Tracking\caviar_recovery_ramp manual tracks.mat')

% calculate avg distance from food well each time... or avg fly position...

x = mean(trackPoints(1,:,:),3,'omitnan');
y = mean(trackPoints(2,:,:),3,'omitnan');


figure
scatter(x,y)

x = (trackPoints(1,:,:));
y = mean(trackPoints(2,:,:),3,'omitnan');

% x_dist = 

temp = data.occupancy.temp;


%% Manually tracked data loading
clear
%load initial data set
load('G:\My Drive\Jeanne Lab\DATA\Manual Tracking\Tracking Data by Parts\caviar_recovery_ramp manual tracks - Part A.mat')
data = trackPoints;
% Repeat for B-F
load('G:\My Drive\Jeanne Lab\DATA\Manual Tracking\Tracking Data by Parts\caviar_recovery_ramp manual tracks - Part F.mat')
loc = (~isnan(trackPoints(1,:,1)));
data(:,loc,:) = trackPoints(:,loc,:);

loc = (~isnan(data(1,:,1)));
figure; plot(loc)
% save data
trackPoints = data;
save('G:\My Drive\Jeanne Lab\DATA\Manual Tracking\9.19.22 Arena C individual fly tracking\manual tracks.mat','trackPoints');

%% Load data:
figFolder = 'G:\My Drive\Jeanne Lab\DATA\Manual Tracking\9.19.22 Arena C individual fly tracking\';
load([figFolder 'manual tracks.mat']);
load([figFolder 'caviar_recovery_rampC timecourse data.mat']);


%% Distance to each well over the course of the experiment...
pix2mm = 12.8;
loc = ~isnan(trackPoints(1,:,1));
ncountedFrames = sum(~isnan(trackPoints(1,:,1)));
X = squeeze(trackPoints(1,:,:)); % all X positions per frame
Y = squeeze(trackPoints(2,:,:)); % all Y positions per frame    
manual_x = mean(trackPoints(1,:,:),3,'omitnan');
manual_y = mean(trackPoints(2,:,:),3,'omitnan');
temperature = data.occupancy.temp;
time = data.occupancy.time;

% distance from well 
well_loc = 1; %manually adjust for this trial...
well_C = data.wellcenters(:,well_loc);
c1 = well_C(1);
c2 = well_C(2);

% Find distance to well center
dist2well = sqrt((X-c1).^2 + (Y-c2).^2)./pix2mm;
avg_dist = mean(dist2well,2,'omitnan');
loc = ~isnan(avg_dist);


% FIGURE
rows = 4;
cols = 1;
sb(1).idx = 1;
sb(2).idx = 2:4;

fig = figure; 
subplot(rows,cols,sb(1).idx)
    plot(time, temperature,'Color','w','LineWidth',1.5)
    ylabel('T (\circC)')

subplot(rows,cols,sb(2).idx); hold on
    plot(time,data.occupancy.dist2wells(:,well_loc),'Color','w','linewidth',1)
    plot(time(loc), avg_dist(loc),'Color',Color('red'),'linewidth',1)
    ylabel('distance to food (mm)')
    xlabel('time (min)')

formatFig(fig,true,[rows,cols],sb);

subplot(rows,cols,sb(1).idx)
set(gca,'XColor','k')
set(gca,'TickDir','out')

subplot(rows,cols,sb(2).idx)
set(gca,'TickDir','out','ydir','reverse')
ylabel('proximity to food (mm)')


save_figure(fig,'G:\My Drive\Jeanne Lab\DATA\Manual Tracking\9.19.22 Arena C\manual vs SLEAP','-pdf',true,false)




%%  FLY EATING MANUAL COUNT
clear
load('G:\My Drive\Jeanne Lab\DATA\Manual Tracking\caviar_recovery_ramp manual tracks.mat')

temp = squeeze(trackPoints(1,:,:));
flyNum = sum(~isnan(temp),2);

fig = figure;
scatter(1:length(flyNum),flyNum)














