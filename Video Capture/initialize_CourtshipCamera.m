
function [src, vid] = initialize_CourtshipCamera(src,vid,FPS)

partial_ROI = [0 0 2048 2048];
% src = getselectedsource(vid);
% camera parameters
src.Brightness = 29;
src.Exposure  = 1.5648;
src.FrameRate = FPS;
src.Gain = 1.752491; 
src.Gamma = 1.5338; 
src.Shutter = 11.6188; 
% vid.FramesPerTrigger = inf;
vid.ROIPosition = partial_ROI;
disp('Purple BACK camera initialized')