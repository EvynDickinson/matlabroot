
%% Initialize camera
vid = videoinput('pointgrey', 1, 'F7_Raw8_2048x2048_Mode0');
    partial_ROI = [0 0 2048 2048];
    src = getselectedsource(vid);
    % camera parameters
    src.Brightness = 29; %11.127068;
    src.Exposure  = 1.5648;
    src.FrameRate = 30;
    src.Gain = 1.752491; %6;
    src.Gamma = 1.5338; % 1.502002;
    src.Shutter = 11.6188; %7.420205;
    vid.FramesPerTrigger = inf;
    vid.ROIPosition = partial_ROI;
    disp('Purple BACK camera initialized')

%% Open folder and run trials

baseFolder = 'F:\Evyn\DATA\09.23.2024\Video Testing\';
% save videos
cd(baseFolder) % set path to video folder

hz = 30; %FPS 
trial = 5; %recording length
totalLength = 16;
nsamples = ceil((totalLength*60)/trial);

get_samples_v3(trial*hz*nsamples,(trial*hz)) % vid duration, frames per vid