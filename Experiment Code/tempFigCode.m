% Make video of each of the temperature profiles
clear
close all
baseFolder = getCloudPath;
data_path = [baseFolder '09.19.2022/Arena C/caviar_recovery_rampC timecourse data.mat'];
load(data_path);
% "03.13.2024\Arena B\C2_F_LRR_25-17B timecourse data.mat";
% '05.08.2024/Arena A/C1_F_LRR_25-17_n2A timecourse data.mat'



% load('09.19.2022\Arena C\caviar_recovery_rampC timecourse data.mat'); %arena C, slow, ramp 2?
% load('DATA\02.01.2023\Arena B\linear_recovery_F_caviarB timecourse data.mat'); %ramp 2, Fast, Arena B

% parameters
% ramp = 2;
% buffer = 10; % how many minutes before the start of the ramp do we want to use?
x_speed = 300; % how many times faster than the fps to write the video (e.g. 20x)


% Get timing starts for ramp
T = data.T;
tempPoints = getTempTurnPoints(expData.parameters.protocol);
% rampDurations = [(tempPoints.down(:,1))/180, (tempPoints.up(:,2))/180]; %in minutes

% % Run all videos but the last one
% % find the videos that are required to cover the selected ramps
% startVid.loc = 1; %tempPoints.down(ramp,1)-(buffer*180);
% startVid.num = 1; %T.vidNums(startVid.loc);
% startVid.frame = 1; %T.vidFrame(startVid.loc);
%
% endVid.loc = 133080; %tempPoints.up(ramp,2)+(buffer*180);
% endVid.num = 37; %T.vidNums(endVid.loc);
% endVid.frame = 3597; %T.vidFrame(endVid.loc);
%
% videos = startVid.num:endVid.num;
%
% vidDuration = ((endVid.loc-startVid.loc)/180) / x_speed;
% fprintf(['\n Video duration: ' num2str(vidDuration) ' mins\n'])
videoFPS = 3*x_speed;

% total_time = T.time(startVid.loc:endVid.loc) - T.time(startVid.loc);
% total_temp = T.temperature(startVid.loc:endVid.loc);
% time_limits = [-2,ceil(total_time(end))];
% temp_limits = [floor(min(total_temp)),ceil(max(total_temp))];

%% Write the video
disp(datetime)
% vidDir = [baseFolder expData.parameters.date '/' expName '_'];
% dataDir = [baseFolder 'FlyRampVideos\'];
vidName = ['/Users/evyndickinson/Data/Data Movies/' expData.parameters.protocol '.avi'];

% Open video object to which you will write new frames
v = VideoWriter(vidName,'Motion JPEG AVI');
v.FrameRate = videoFPS;
v.Quality = 85;
open(v);
roi = 1:T.frame(end);

fig = figure; set(fig,'pos',[-1030 279 772 200],'color','k','Visible','Off'); %'Visible','Off'
for frame = 1:roi(end)

    plot(T.time(roi),T.temperature(roi),'linewidth', 2,'Color','w')
    hold on
    v_line(T.time(frame),'red','-',2)
    set(gca,'color','k','xcolor','w','ycolor','w','LineWidth',2)
    set(gca,'FontName','arial','FontSize', 15,'box','off')
    xlabel('time (min)')
    ylabel('temp (\circC)')
    hold off

    % save the current figure as new frame in the video
    f = getframe(fig);
    writeVideo(v, f)
    if rem(frame,1000)==0
        disp(['Frame ' num2str(frame)])
    end
end
disp('Done')
disp(datetime)

close(v) %close the movie writer object
close(fig) %close the figure



%% Write the video
disp(datetime)
vidDir = [baseFolder expData.parameters.date '/' expName '_'];
dataDir = [baseFolder 'FlyRampVideos\'];
vidName = [dataDir 'Fast Temp Ramp.avi'];

% % Set axis limits for the selected arena
% x = data.centre(1);
% y = data.centre(2);
% r = data.r;
% xlimit = [x-(r+50),x+(r+50)];
% ylimit = [y-(r+50),y+50+r];

skipSize = 1;

% r = 4; c = 1;
% sb(1).idx = 1;
% sb(2).idx = 2:4;


% Open video object to which you will write new frames
v = VideoWriter(vidName,'Motion JPEG AVI');
v.FrameRate = videoFPS;
v.Quality = 85;
open(v);

% Get the start frame info & image data
oldVid = T.vidNums(startVid.loc);
vid_frame = T.vidFrame(startVid.loc);
videoPath = [vidDir num2str(oldVid) '.avi'];
movieInfo = VideoReader(videoPath); %read in video

fig = figure; set(fig,'pos',[-1030 279 772 200],'color','k'); %'Visible','Off'
for frame = startVid.loc:skipSize:endVid.loc
    vid = T.vidNums(frame);
    vid_frame = T.vidFrame(frame);
    if vid>oldVid
        videoPath = [vidDir num2str(vid) '.avi'];
        movieInfo = VideoReader(videoPath); %read in video
        oldVid = vid;
    end
    currentImg = rgb2gray(read(movieInfo,vid_frame));

    time = T.time(startVid.loc:frame) - T.time(startVid.loc);
    temp = T.temperature(startVid.loc:frame);

    % PLOT:
    set(fig,'color','k')
    subplot(r,c,sb(2).idx)
    % display current image
    imshow(currentImg)
    xlim(xlimit); ylim(ylimit);
    % display temperature at current frame
    subplot(r,c,sb(1).idx);
    plot(time, temp, 'linewidth', 2,'Color','w')
    hold on
    scatter(time(end), temp(end),50,'red','filled')
    xlim(time_limits)
    ylim(temp_limits)
    set(gca,'color','k','xcolor','w','ycolor','w','LineWidth',2)
    set(gca,'FontName','arial','FontSize', 15,'box','off')
    xlabel('time (min)')
    ylabel('temp (\circC)')
    hold off

    % save the current figure as new frame in the video
    f = getframe(fig);
    writeVideo(v, f)
%     pause(0.001)
end
disp('Done')
disp(datetime)

close(v) %close the movie writer object
close(fig) %close the figure


%% Make figure of the whole data set and then the selected time frames

r = 3; c = 1;
sb(1).idx = 1;
sb(2).idx = 2:3;

fig = getfig('',1,[866 796]);
subplot(r,c,sb(1).idx)
    plot(T.time,T.temperature,'color','w', 'LineWidth',2)
    ylabel('temp (\circC)')
subplot(r,c,sb(2).idx)

    plot(T.time,T.dist2wells(:,1),'color','w', 'LineWidth',1)
    % plot(T.time,T.dist2food,'color','w', 'LineWidth',1)

    ylabel('proximity to food (mm)')
    xlabel('time (min)')
    set(gca, 'ydir', 'reverse')
    % h_line(18.1,'grey',':',1) %36.2
formatFig(fig,true,[r,c], sb);

% plot the video range
subplot(r,c,sb(1).idx); hold on
y = rangeLine(fig,0);
timeROI = [T.time(startVid.loc), T.time(endVid.loc)];
plot(timeROI, [y,y], 'color','y')


% Save figure
save_figure(fig, [vidName(1:end-4) '.png'],'-png',true,true)
