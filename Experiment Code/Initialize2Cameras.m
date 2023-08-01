function [src, vid] = Initialize2Cameras(FPS,cam_sel)
    % Initialize the cameras that will be used for the experiment
    %  "purple camera 1"
%     if cam_sel==1
        imaqreset
        
        vid1 = videoinput('pointgrey', 1, 'F7_Raw8_2048x2048_Mode0');
        partial_ROI = [66 44 1920 1940];
        src1 = getselectedsource(vid1);
        % camera parameters
        src1.Brightness = 11.127068;
        src1.Exposure  = 2;
        src1.FrameRate = 3; %FPS;
        src1.Gain = 6;
        src1.Gamma = 1.502002;
        src1.Shutter = 7.420205;
        vid1.FramesPerTrigger = inf;
        vid1.ROIPosition = partial_ROI;
        disp('Camera initialized')

        preview(vid1)
        
    
        % "BLUE CAMERA"
        vid2 = videoinput('pointgrey', 2, 'F7_Raw8_2048x2048_Mode0');
        partial_ROI = [66 72 1920 1940];
        src2 = getselectedsource(vid2);
        % camera parameters
        src2.Brightness = 15;
        src2.Exposure  = 2;
        src2.FrameRate = 3; %FPS;
        src2.Gain = 6;
        src2.Gamma = 1.502002;
        src2.Shutter = 7.420205;
        vid2.FramesPerTrigger = inf;
        vid2.ROIPosition = partial_ROI;
        disp('Camera initialized')

        preview(vid2)

    
    
    
    
    
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