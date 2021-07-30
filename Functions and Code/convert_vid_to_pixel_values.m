function [fly_images, pret_images] = convert_vid_to_pixel_values(video_name, parameters, color)
%
% [fly_images, pret_images] = convert_vid_to_pixel_values(video_name, parameters, color)
% 
% [fly_images, pret_images](frame#).frames(ROI-X, ROI-Y) 
% load and read all the frames of a video (video_name)
% save the grayscale pixel values for the image to fly_images
% add a buffer to the beginning/end to match the total stimulus unit for
% graphing
%
% Inputs:
% 'video_name' [filepath and name for the video to be analyzed]
% Outputs:
% 'fly_images' [structure with the pixel light values for each frame in vid]
% 'pret_images' [buffered structure with the pixel values for each frame in vid]
% 'color' [binary option for color on (1) or off (0)]
%
% ES Dickinson, University of Washington, Dec 2018

% % Create directory for images
% workingDir = 'D:\Evyn Data Files\11.16.18\Fly 1_0\Processed Video';
% % mkdir(workingDir)
% mkdir(workingDir,'Images') %set this to a specific RC stim
% 
% 
parameters;

if ~exist(video_name, 'file')
  error('Video file "%s" does not exist', video_name);
end
% try
%   videoobj = VideoReader(filename);
% catch
%   error('File "%s" cannot be read as a video', filename);
% end
% img = readFrame(videoobj);   %no frame2im !

% % Create a VideoReader to use for reading frames from the file.
flyVideo = VideoReader(video_name);
% % video_name = ['11162018_fly1_0 R1C1 Cam-A str-cw-0 sec']; %include cam info, etc. basically the vid name


vid_length = flyVideo.Duration*flyVideo.FrameRate;


blackandwhite(vid_length).frames = struct;
colorscale(vid_length).frames = struct; 
ii = 1;
while hasFrame(flyVideo)
   [X,img] = imread(flyVideo);%readFrame(flyVideo);
   blackandwhite(ii).frames = rgb2gray(img);
   colorscale(ii).frames = img;
   ii = ii+1;
end


% Color or Grayscale
if nargin == 2 %choice made
    switch color
        case 1 %color scale
            fly_images = colorscale;
        case 0 %greyscale
            fly_images = blackandwhite;      
    end
else
    fly_images = blackandwhite;  
end
    
    
% % add buffer of grayscale values at the beginning and end of the video
% start_buffer_length = (parameters.OL_time-parameters.basler_delay)*parameters.Basler_fps;
% end_buffer_length = (parameters.OL_time-parameters.basler_length+parameters.basler_delay)*parameters.Basler_fps;
% start_frames = size(fly_images,2);
% image_size = size(fly_images(1).frames);
% buffer_frame = uint8(zeros(image_size));

% add buffer of grayscale values at the beginning and end of the video
start_buffer_length = 2;
end_buffer_length = 2;
start_frames = size(fly_images,2);
image_size = size(fly_images(1).frames);
buffer_frame = uint8(zeros(image_size));

% create a new structure with buffered start and stop for plotting
for i = 1:start_buffer_length
pret_images(i).data = buffer_frame;
end
strt = start_buffer_length+1;
stp = start_frames+start_buffer_length-1;
for i = strt:stp
pret_images(i).data = fly_images(i-start_buffer_length+1).frames;
end
strt = stp+1;
stp = strt+end_buffer_length-1;
for i = strt:stp
pret_images(i).data = buffer_frame;
end



% figure;
% imshow(FLY_IMAGE(300).frames)
end




