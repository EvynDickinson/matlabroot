
function [src, vid] = InitializeCamera(FPS)
    % image cropping parameters
    full_ROI = [0 0 2048 2048];
    % always want to keep the ROI size 1920 x 1920 if possible. The offsets
    % can change to match a moving plate if needed
    partial_ROI = [44 80 1920 1920];
%                   [72 94 1901 1918]; 
%                   [31 71 1949 1949]
    % Load in the camera / open preview
    vid = videoinput('pointgrey', 1, 'F7_Raw8_2048x2048_Mode0');
    src = getselectedsource(vid);

    % camera parameters
    src.Brightness = 11.127068;
    src.Exposure  = 2;
    src.FrameRate = FPS;
    src.Gain = 6.0;
    src.Gamma = 1.5;
    src.Shutter = 5;

    vid.FramesPerTrigger = inf;
    % vid.ROIPosition = rectangular_ROI;
    %   vid.ROIPosition = full_ROI;
    vid.ROIPosition = partial_ROI;
    disp('Camera initialized')
end