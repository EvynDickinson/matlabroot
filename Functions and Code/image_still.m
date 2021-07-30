


function hf = image_still(VideoName, timestamp, zoom)
% hf = displayimage(VideoName, timestamp, zoom)
% Min input is VideoName, complete with address to find it
% Defaults: 
% timestamp = 0.5;
% zoom = 1;


% timestamp = 0.5;
% zoom = 1.8;

% pull the data from the video:
flyVideo = VideoReader(VideoName);
flyVideo.CurrentTime = timestamp;
picture = readFrame(flyVideo);

% constrast increased
k = imadjust(rgb2gray(picture));

% % increase contrast
J = imsharpen(k,'Radius',2,'Amount',1);
% % increase contrast
% J = k;

I = imresize(J, zoom);
postPic = (I);


Width = size(postPic, 2);
Height = size(postPic, 1);

hf = figure;
% set(hf,'position',[39, 234, Width, Height]);

% set(hf,'position',[-655, 72, Width, Height]);
set(hf, 'Name', VideoName(end-45:end))

% Show the image!
imshow(postPic);
% set(hf,'position',[-1200, 234, Width, Height]);



% VideoName = 'D:\Evyn Data Files\11.1.19\Fly 3_0\Raw Video\11012019_fly3_0 R1C14 Cam-B str-ccw-0.72 sec.avi';










