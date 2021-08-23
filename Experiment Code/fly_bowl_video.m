clear
imaqreset
warning off
%% Video and camera parameters: 
basepath = 'C:\Users\jeannelab\Documents\Evyn\DATA\';
dirName = strrep(datestr(datetime, 'mm-dd-yyyy'), '-', '.');
video_path = [basepath, dirName, '\'];
if ~isfolder(video_path); mkdir(video_path); end

% vid parameters
full_ROI = [0 0 2048 2048];
partial_ROI = [493 693 1241 1222];
rectangular_ROI = [552 586 628 891]; %srectangular arena space

% Load in the camera / open preview
vid = videoinput('pointgrey', 1, 'F7_Raw8_2048x2048_Mode0');
src = getselectedsource(vid);

% camera parameters
src.Brightness = 11.127068;
src.Exposure  = 0;
src.FrameRate = 3;
src.Gain = 4.0;
src.Gamma = 1.5;
src.Shutter = 5;

vid.FramesPerTrigger = inf;
% vid.ROIPosition = rectangular_ROI;
%   vid.ROIPosition = full_ROI;
  vid.ROIPosition = partial_ROI;

preview(vid);

%% Acquire SINGLE save video!

% ------- video params ----------------
VID_LENGTH = 2; 
name_modifier = 'testvid';
% -------------------------------------

closepreview
vid.LoggingMode = 'disk';   %'disk&memory'; 
diskLogger = VideoWriter([video_path name_modifier '.avi'], 'Motion JPEG AVI');  
diskLogger.Quality = 85;
% diskLogger = VideoWriter([video_path 'testingNOcompression.avi'], 'Grayscale AVI');
diskLogger.FrameRate = 3; 
vid.DiskLogger = diskLogger;

start(vid)
pause(VID_LENGTH)  
stop(vid)

fprintf('\n All wrapped\n')

preview(vid)
 
%% DEMO VIDEO QUALITY FOR LATER ANALYSIS
vidPath = [basepath dirName '\' name_modifier '.avi'];
frame = 1;

% video parameters
 
movieInfo = VideoReader(vidPath); %read in video

% extract parameters from video
nframes = movieInfo.Duration * movieInfo.FrameRate;
height = movieInfo.Height;
width = movieInfo.Width;
occCumSum = zeros(height, width); %setup blank occupation
n_tot = nframes;

% default analysis params:
pixel_thresh = 60;
cluster_thresh = 5;
Img = read(movieInfo,frame);
demoImg = processOccImg(Img, false, pixel_thresh, cluster_thresh);

f = figure; imshowpair(Img, demoImg,'montage'); title('Thresholding test')
          

imtool(demoImg)




%% Acquire SERIES of videos 
 
% ------- video params ----------------
N_repeats = 10;
VID_LENGTH = 60;
name_modifier = 'MovementTest_N19_20C';
% -------------------------------------

% find the current video count for this name modifier and use as offset number:
list_dirs = dir([folder, '\*' name_modifier '*.avi']); %only videos
offset = length(list_dirs);

closepreview    
for i = offset+1:offset+N_repeats
    
    diskLogger = VideoWriter([video_path name_modifier '_' num2str(i) '.avi'], 'Motion JPEG AVI');
    vid.LoggingMode = 'disk';   %'disk&memory'; 
    diskLogger.Quality = 65;
    diskLogger.FrameRate = 20; 
    vid.DiskLogger = diskLogger;  
    vid_dur = 60;   % duration of video recording

    start(vid)
    pause(vid_dur)  
    stop(vid)

    fprintf(['\n' num2str(i) ' / ' num2str(offset+N_repeats)  '\n']) 

end
preview(vid)


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
 


 


