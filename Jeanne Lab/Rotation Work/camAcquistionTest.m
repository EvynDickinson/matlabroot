
%%  Test video acquisition for the hamamatsu camera
% propinfo(vid)     --> properties of the camera
% triggerinfo(vid)  --> trigger options for the camera
% propinfo(src)     --> camera source properties
vid_length = 2;
fps = 60;

vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');
vid.LoggingMode = 'disk&memory';         % save the data to memory then disk
vid.FramesPerTrigger = vid_length*fps;          % frames per input trigger
vid.ROIPosition = [0 0 512 512];  % x offset, y offset, width, height [0 0 2048 2048]


triggerconfig(vid, 'manual') % set camera to manual triggermode
triggerconfig(vid, 'hardware', 'RisingEdge', 'EdgeTrigger')
triggerconfig(vid, 'immediate', 'none', 'none') 

imaqmontage(vid)


% intensity for preview options
imaqmex('feature', '-previewFullBitDepth', true);
imaqmex('feature', '-previewFullBitDepth', false);

preview(vid)

closepreview(vid)
imaqreset

src = getselectedsource(vid);
src.ExposureTime = 0.01;    % cam exposure
% src.TriggerConnector = 'bnc';


preview(vid)

imaqreset



%% test record video...

%think about the best way to do this...

video_name = 'TestVideo3';
diskLogger = VideoWriter(video_name, 'Grayscale AVI'); % 'Motion JPEG 2000' <--8 bit not 16
diskLogger.FrameRate = fps;      %fps for the camera 
% diskLogger.FrameCount 

vid.DiskLogger = diskLogger;  

tic
% aquire video:
start(vid)
pause(3)
%ensure that all frames 
while (vid.FramesAcquired ~= vid.DiskLoggerFrameCount) 
    pause(.1)
end
stop(vid)

toc


start(vid)
%triggers
queueOutputData(s, data) 
startBackground(s) 

pause(5)
while (vid.FramesAcquired ~= vid.DiskLoggerFrameCount) 
    pause(.1)
end
stop(vid)






 %% Demo the mj2 video file format
 
 
videoFReader = vision.VideoFileReader('E:\Evyn\TestVideo1.mj2');

videoPlayer = vision.VideoPlayer;

while ~isDone(videoFReader)
  videoFrame = videoFReader();
  videoPlayer(videoFrame);
  pause(0.1)
end
 
 
 
 
 
 
 
 
 
%% Basler version:
% Setup source and vid objects for Basler Camera
vid = videoinput('gentl', cam_spec.num, 'Mono8');
triggerconfig(vid, 'hardware');                 % set trigger to come from hardware i.e. line in
vid.LoggingMode = 'disk'; %'disk&memory';
vid.FramesPerTrigger = cam_spec.FramesPerTrigger;
vid.ROIPosition = cam_spec.ROI;            % fly-only video acq.
Basler_src = vid.Source;
% Basler_src = getselectedsource(vid);
Basler_src.TriggerMode = 'Off';

% Set up for input wires from Basler Cam
Basler_src.LineSelector = 'Line4';             % brings up settings for line4
Basler_src.LineMode = 'output';                % should be 'output'; Basler cam info
Basler_src.LineInverter = 'False';             % should be 'False'
Basler_src.LineSource = 'ExposureActive';      % send out signal of when the exposure was active on basler
Basler_src.LineSelector = 'Line3';             % brings up settings for line3
Basler_src.LineMode = 'input';                 % should be 'output'; sends trig info to cam
Basler_src.TriggerSelector = 'FrameStart';     % start frame with trigger
Basler_src.TriggerMode = 'Off';
Basler_src.LineSelector = 'Line3';             % brings up settings for line3
Basler_src.TriggerActivation = 'RisingEdge';
Basler_src.TriggerMode = 'On';
Basler_src.GainAuto = 'Off';
Basler_src.ExposureTime = cam_spec.ExposureTime;                % Exposure setting for Basler
Basler_src.Gain = cam_spec.Gain;          % Lighting gain
Basler_src.Gamma = cam_spec.Gamma;         % White enhancement on video
