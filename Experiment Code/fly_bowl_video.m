clear
imaqreset

%% Path selection: 
basepath = 'C:\Users\jeannelab\Documents\Evyn\DATA\';
dirName = strrep(datestr(datetime, 'mm-dd-yyyy'), '-', '.');
video_path = [basepath, dirName, '\'];
if ~isfolder(video_path); mkdir(video_path); end

% vid parameters
full_ROI = [0 0 2048 2048];
partial_ROI = [493 693 1241 1222];
rectangular_ROI = [552 586 628 891]; %srectangular arena space

%% Load in the camera / open preview

vid = videoinput('pointgrey', 1, 'F7_Raw8_2048x2048_Mode0');
src = getselectedsource(vid);

% camera parameters
src.Brightness = 11.127068;
src.Exposure  = 0;
src.FrameRate = 20;
src.Gain = 1.75;
src.Gamma = 1.4795;
src.Shutter = 7.649936;

vid.FramesPerTrigger = inf;
% vid.ROIPosition = rectangular_ROI;
%   vid.ROIPosition = full_ROI;
  vid.ROIPosition = partial_ROI;

preview(vid);

%% Acquire and save video!

closepreview

vid.LoggingMode = 'disk';   %'disk&memory'; 

diskLogger = VideoWriter([video_path 'blackbackgroundarena.avi'], 'Motion JPEG AVI');  %'Grayscale AVI'
diskLogger.Quality = 85;
% diskLogger = VideoWriter([video_path 'testingNOcompression.avi'], 'Grayscale AVI');
diskLogger.FrameRate = 20; 
vid.DiskLogger = diskLogger;

vid_dur = 10;   % duration of video recording

start(vid)
pause(vid_dur)  
stop(vid)

fprintf('All wrapped')
 
 
 %% Quick save a series of videos
 N_repeats = 10;
 
for i = 1:N_repeats
    
    diskLogger = VideoWriter([video_path 'DummyVid_15C_' num2str(i) '.avi'], 'Motion JPEG AVI');
    diskLogger.Quality = 85;
    diskLogger.FrameRate = 20; 
    vid.DiskLogger = diskLogger;  
    vid_dur = 60;   % duration of video recording

    start(vid)
    pause(vid_dur)  
    stop(vid)

    fprintf(['\n' num2str(i) ' / ' num2str(N_repeats)  '\n'])

end
 

%% Background picture (approx 5 images to average together)

vid.LoggingMode = 'disk';   %'disk&memory'; 

diskLogger = VideoWriter([video_path 'background _vid.avi'], 'Grayscale AVI');
diskLogger.FrameRate = 5; 
vid.DiskLogger = diskLogger;

vid_dur = 1;   % duration of video recording

start(vid)
pause(vid_dur)  
stop(vid)

fprintf('All wrapped')
 
 


