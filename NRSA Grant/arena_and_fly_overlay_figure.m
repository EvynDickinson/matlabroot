


clear
close all
data_path = "G:\My Drive\Jeanne Lab\DATA\03.13.2024\Arena B\C2_F_LRR_25-17B timecourse data.mat";
load(data_path);
T = data.T;

% load('G:\My Drive\Jeanne Lab\DATA\09.19.2022\Arena C\caviar_recovery_rampC timecourse data.mat'); %arena C, slow, ramp 2?
% load('G:\My Drive\Jeanne Lab\DATA\02.01.2023\Arena B\linear_recovery_F_caviarB timecourse data.mat'); %ramp 2, Fast, Arena B

% % parameters
% ramp = 2;
% buffer = 10; % how many minutes before the start of the ramp do we want to use?
% x_speed = 300; % how many times faster than the fps to write the video (e.g. 20x)
% 
% 
% % Get timing starts for ramp
% 
% tempPoints = getTempTurnPoints(expData.parameters.protocol);
% rampDurations = [(tempPoints.down(:,1))/180, (tempPoints.up(:,2))/180]; %in minutes
% 
% 
% % find the videos that are required to cover the selected ramps
% startVid.loc = tempPoints.down(ramp,1)-(buffer*180);
% startVid.num = T.vidNums(startVid.loc);
% startVid.frame = T.vidFrame(startVid.loc);
% 
% endVid.loc = tempPoints.up(ramp,2)+(buffer*180);
% endVid.num = T.vidNums(endVid.loc);
% endVid.frame = T.vidFrame(endVid.loc);
% 
% videos = startVid.num:endVid.num;
% 
% vidDuration = ((endVid.loc-startVid.loc)/180) / x_speed;
% fprintf(['\n Video duration: ' num2str(vidDuration) ' mins\n'])
% videoFPS = 3*x_speed;
% 
% total_time = T.time(startVid.loc:endVid.loc) - T.time(startVid.loc);
% total_temp = T.temperature(startVid.loc:endVid.loc);
% time_limits = [-2,ceil(total_time(end))];
% temp_limits = [floor(min(total_temp)),ceil(max(total_temp))];


% %% Write the video
% disp(datetime)
baseFolder = getCloudPath;
vidDir = [baseFolder expData.parameters.date '/' expName '_'];
dataDir = [baseFolder 'FlyRampVideos\'];
% vidName = [dataDir expData.parameters.date ' ' expName ' ' data.name  ' Ramp ' num2str(ramp) '.avi'];
% 
% % Set axis limits for the selected arena
% x = data.centre(1);
% y = data.centre(2);
% r = data.r;
% xlimit = [x-(r+50),x+(r+50)];
% ylimit = [y-(r+50),y+50+r];
% 
% skipSize = 1;
% 
% r = 4; c = 1;
% sb(1).idx = 1;
% sb(2).idx = 2:4;
% 
% 
% % Open video object to which you will write new frames
% v = VideoWriter(vidName,'Motion JPEG AVI');
% v.FrameRate = videoFPS;
% v.Quality = 85;
% open(v);
% 
% % Get the start frame info & image data
% oldVid = T.vidNums(startVid.loc);
% vid_frame = T.vidFrame(startVid.loc);
% videoPath = [vidDir num2str(oldVid) '.avi'];
% movieInfo = VideoReader(videoPath); %read in video
% 
% fig = figure; set(fig,'pos',[-1030 279 772 1009],'color','k'); %'Visible','Off'
% for frame = startVid.loc:skipSize:endVid.loc
%     vid = T.vidNums(frame);
%     vid_frame = T.vidFrame(frame);
%     if vid>oldVid
%         videoPath = [vidDir num2str(vid) '.avi'];
%         movieInfo = VideoReader(videoPath); %read in video
%         oldVid = vid;
%     end
%     currentImg = rgb2gray(read(movieInfo,vid_frame));
% 
%     time = T.time(startVid.loc:frame) - T.time(startVid.loc);
%     temp = T.temperature(startVid.loc:frame);
% 
%     % PLOT:
%     set(fig,'color','k')
%     subplot(r,c,sb(2).idx)
%     % display current image
%     imshow(currentImg)
%     xlim(xlimit); ylim(ylimit); 
%     % display temperature at current frame
%     subplot(r,c,sb(1).idx);
%     plot(time, temp, 'linewidth', 2,'Color','w') 
%     hold on
%     scatter(time(end), temp(end),50,'red','filled')
%     xlim(time_limits)
%     ylim(temp_limits)
%     set(gca,'color','k','xcolor','w','ycolor','w','LineWidth',2)
%     set(gca,'FontName','arial','FontSize', 15,'box','off')
%     xlabel('time (min)')
%     ylabel('temp (\circC)')
%     hold off
% 
%     % save the current figure as new frame in the video
%     f = getframe(fig);
%     writeVideo(v, f)
% %     pause(0.001)
% end
% disp('Done')
% disp(datetime)
% 
% close(v) %close the movie writer object
% close(fig) %close the figure


%% Make figure of the whole data set and then the selected time frames
blkbgd = false;
[foreColor,backColor] = formattingColors(blkbgd);

r = 3; c = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;

fig = getfig('',1,[866 796]);
subplot(r,c,sb(1).idx)
    plot(T.time,T.temperature,'color',foreColor, 'LineWidth',2)
    ylabel('temp (\circC)')
subplot(r,c,sb(2).idx)
    wellLoc = data.foodwell;
    y = T.occ_P(:,wellLoc).*100;
    plot(T.time,y,'color',foreColor, 'LineWidth',1)
    % plot(T.time,T.dist2food,'color','w', 'LineWidth',1)

    ylabel('occupancy of food region (%)')
    xlabel('time (min)')
    h_line(14.4,'grey',':',1) 
formatFig(fig,blkbgd,[r,c], sb);

% plot the video range
subplot(r,c,sb(1).idx); hold on
y = rangeLine(fig,0);
timeROI = [T.time(startVid.loc), T.time(endVid.loc)];
plot(timeROI, [y,y], 'color',Color('gold'))


% Save figure
save_figure(fig, [vidName(1:end-4) '.png'],'-png',true,true)

%% 

time_point = 247.015;
tp = getTempTurnPoints(expData.parameters.protocol);
[~,exp_idx]  = min(abs((time_point-T.time)));

vid_num = T.vidNums(exp_idx);
frame_num = T.vidFrame(exp_idx);

% load the video 
baseFolder = getDataPath(2,2);
vidDir = [baseFolder expData.parameters.date '/' expName '_'];
videoPath = [vidDir num2str(vid_num) '.avi'];
movieInfo = VideoReader(videoPath); %read in video

% Set axis limits for the selected arena
x = data.centre(1);
y = data.centre(2);
r = data.r;
xlimit = [x-(r+50),x+(r+50)];
ylimit = [y-(r+50),y+50+r];

% load the image
currentImg = rgb2gray(read(movieInfo,frame_num));


% PLOT image:
fig = getfig('',1);
set(fig,'color','k')
imshow(currentImg)
xlim(xlimit); ylim(ylimit); 
hold on
imcontrast

% a = imadjust(currentImg,[53/255,149/255]); % adjust the contrast
    
% plot a yellow point for each fly
x = data.x_loc(exp_idx,:);
y = data.y_loc(exp_idx,:);
scatter(x,y,50,'y','filled')
scatter(x,y,10,'r','filled')

% plot the food region of interest
centre = data.centre;
viscircles(centre', r, 'color', 'w');

well_center = data.wellcenters(:,data.foodwell);
viscircles(well_center', data.radii, 'color', 'w');

save_figure(fig, [getCloudPath, 'FlyRampVideos/' expName ' still frame ' num2str(exp_idx)], '-pdf')

%% USING DATA POST 4.2  -- Make figure of the whole data set and then the selected time frames
% blkbgd = false;
[foreColor,backColor] = formattingColors(blkbgd);
exp = 1;
trial = 1;

T = data(exp).data(trial).data.T;

r = 3; c = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;


fig = getfig('',1,[866 796]);
subplot(r,c,sb(1).idx)
    plot(T.time,T.temperature,'color',foreColor, 'LineWidth',2)
    ylabel('temp (\circC)')
subplot(r,c,sb(2).idx)
    wellLoc = data(exp).data(trial).data.foodwell;
    y = T.occ_P(:,wellLoc).*100;
    plot(T.time,y,'color',foreColor, 'LineWidth',1)
    % plot(T.time,T.dist2food,'color','w', 'LineWidth',1)

    ylabel('occupancy of food region (%)')
    xlabel('time (min)')
    h_line(14.4,'grey',':',1) 
formatFig(fig,blkbgd,[r,c], sb);

% select time point from exp aboves 

time_point = 247.015;
[~,exp_idx]  = min(abs((time_point-T.time)));

vid_num = T.vidNums(exp_idx);
frame_num = T.vidFrame(exp_idx);

% load the video 
% vidDir = [getDataPath(2,2) data(exp).T.Date{trial} '/' data(exp).T.ExperimentID{trial} '_'];
vidDir = [getCloudPath data(exp).T.Date{trial} '/' data(exp).T.ExperimentID{trial} '_'];
videoPath = [vidDir num2str(vid_num) '.avi'];
movieInfo = VideoReader(videoPath); %read in video

% Set axis limits for the selected arena
x = data(exp).data(trial).data.centre(1);
y = data(exp).data(trial).data.centre(2);
r = data(exp).data(trial).data.r;
xlimit = [x-(r+50),x+(r+50)];
ylimit = [y-(r+50),y+50+r];

% load the image
currentImg = rgb2gray(read(movieInfo,frame_num));


% PLOT image:
fig = getfig('',1);
set(fig,'color','k')
imshow(currentImg)
xlim(xlimit); ylim(ylimit); 
hold on
imcontrast
    
% plot a yellow point for each fly
x = data(exp).data(trial).data.x_loc(exp_idx,:);
y = data(exp).data(trial).data.y_loc(exp_idx,:);
scatter(x,y,50,'y','filled')
scatter(x,y,10,'r','filled')

% plot the food region of interest
centre = data(exp).data(trial).data.centre;
viscircles(centre', r, 'color', 'w');

well_center = data(exp).data(trial).data.wellcenters(:,data(exp).data(trial).data.foodwell);
viscircles(well_center', data(exp).data(trial).data.radii, 'color', 'w');

% Show 50% edge circle: 
R = data(exp).data(trial).data.r; % radius of the arena in pixels
innerR = R/sqrt(2); % radius of the inner 50% occupancy space
% innerR = R*sqrt(3/4);
viscircles(centre', innerR, 'color', 'w');

save_figure(fig, [getCloudPath, 'FlyRampVideos/' data(exp).T.ExperimentID{trial} ' still frame ' num2str(exp_idx)], '-pdf')








