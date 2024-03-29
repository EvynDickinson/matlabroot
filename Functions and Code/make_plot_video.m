


%% Video Parameters (edit these)

video_name = 'VideoName';   % name of the video
video_folder = uigetdir;    % select the folder where you want to save the video
video_fps = 10;             % how many frames per second do you want the video playback?

%% Generate dummy data to plot: (for an example)

x = 1:100;
y = randi(100,1,100);
% figure;
% plot(x,y)

%% Video extracted properties  | plotting formats

% get the mins/maxs for the axes before making video so that the axes don't
% change size throughout the video
xlimit = [min(x),max(x)];
ylimit = [min(y),max(y)]; 

video_path = [video_folder '\' video_name '.avi'];
num_frames = size(x,2); % how many frames are in the video?

% plot formatting options
kolor = [0  0.5020  0.5020];    % teal plotting color
font_selection = 'Arial';       % plotting font
font_size = 15;                 % axes font sizes
LW = 1.5;                       %plot line width

%% Write video

% Open video object to which you will write new frames
v = VideoWriter(video_path, 'Uncompressed AVI');
v.FrameRate = video_fps;
open(v);

% open the figure from which you will grab frames to save
fig = figure; set(fig,'color','w');

for i = 1:num_frames
    % plot the important stuff
    % -----------------------------------------------------
    plot(x(1:i), y(1:i), 'color', kolor,'LineWidth', LW)

    %format figure
    xlim(xlimit);
    ylim(ylimit);
    xlabel('time (s)')
    ylabel('dummy')
    set(gca,'FontName',font_selection,'FontSize',font_size)
    box off
    % -----------------------------------------------------

    % save the current figure as new frame in the video
    f = getframe(fig);
    writeVideo(v, f)
end

close(v) %close the movie writer object
close(fig) %close the figure

