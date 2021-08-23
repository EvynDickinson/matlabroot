
%get base folder pathway
baseFolder = getCloudPath;

params.well_1 = 'Plant';
params.well_2 = 'Empty';
params.well_3 = 'Empty';
params.well_4 = 'Empty';

params.genotype = select_cross; % choose the fly genotypes
params.protocol = select_protocol;  % choose experiment protocol

% move params into structure:
params.date = dirName;
params.expID = videoNames;
params.num = num;


writeExptoExcel(params) % save the data to the experiement list

% resave with same name: 
save([baseFolder dirName '\' videoNames 'dataMat'])


%% 
baseFolder = getCloudPath;

for i = 1:num.vids
    filename = [baseFolder dirName '\analysis\' videoNames '_' num2str(i) ' data'];
    load(filename);
    
    params.well_1 = 'Plant';
    
    save(filename, 'videoData', 'params')
end


%% Fix the temp log information being added...
clear
% manually load the exp data file
% Fix the MAIN file
[~,baseFolder] = getCloudPath(1);

tempLogStart = [(1:num.vids)', tempLogStart];
save([baseFolder, '\', videoNames, 'dataMat.mat'])

% load each vid file & adjust the temp log:
tempLog = readmatrix([baseFolder '\' videoNames '_RampLog']);
roi = [tempLogStart(:,3), tempLogEnd(:,3)];
filepath = [baseFolder, '\analysis\' videoNames '_'];

for i = 1:num.vids
    load([filepath num2str(i) ' data.mat']);
    
    % update file
    loc = tempLog(:,1)>=roi(i,1) & tempLog(:,1)<=roi(i,2);
    videoData.tempLog = tempLog(loc,:);
    
    % save file
    save([filepath num2str(i) ' data.mat'], 'params', 'videoData')
end

%% Video testing


Vq = interpn(tframe,3);
fig = getfig('',1);
hold on
for ii = 1:size(Vq,3)
    imagesc(Vq(:,:,ii))
    pause(0.05)
end


%% Crop out the regions around each food well to compare across trials...

% Demo to have the user click and draw a circle over an image, then blacken outside the circle and crop out the circular portion into a new image.
clc;    % Clear the command window.
fprintf('Beginning to run %s.m ...\n', mfilename);
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing
fontSize = 15;
% Get image.
originalImage = imread('peppers.png');
[rows, columns, numberOfColorChannels] = size(originalImage);
subplot(2, 2, 1);
imshow(originalImage);
axis('on', 'image');
title('Original Image', 'FontSize', fontSize);
% Maximize the window to make it easier to draw.
g = gcf;
g.WindowState = 'maximized';
% Ask user to draw a circle:
uiwait(helpdlg('Please click and drag out a circle.'));
h.Radius = 0;
while h.Radius == 0
	h = drawcircle('Color','k','FaceAlpha',0.4)
	if h.Radius == 0
		uiwait(helpdlg('You double-clicked.  You need to single click, then drag, then single click again.'));
	end
end
% Get coordinates of the circle.
angles = linspace(0, 2*pi, 10000);
x = cos(angles) * h.Radius + h.Center(1);
y = sin(angles) * h.Radius + h.Center(2);
% Show circle over image.
subplot(2, 2, 2);
imshow(originalImage);
axis('on', 'image');
hold on;
plot(x, y, 'r-', 'LineWidth', 2);
title('Original image with circle mask overlaid', 'FontSize', fontSize);
% Get a mask of the circle
mask = poly2mask(x, y, rows, columns);
subplot(2, 2, 3);
imshow(mask);
axis('on', 'image');
title('Circle Mask', 'FontSize', fontSize);
% Mask the image with the circle.
if numberOfColorChannels == 1
	maskedImage = originalImage; % Initialize with the entire image.
	maskedImage(~circleImage) = 0; % Zero image outside the circle mask.
else
	% Mask the image.
	maskedImage = bsxfun(@times, originalImage, cast(mask, class(originalImage)));
end
% Crop the image to the bounding box.
props = regionprops(mask, 'BoundingBox');
maskedImage = imcrop(maskedImage, props.BoundingBox);
% Display it in the lower right plot.
subplot(2, 2, 4);
imshow(maskedImage, []);
% Change imshow to image() if you don't have the Image Processing Toolbox.
title('Image masked with the circle', 'FontSize', fontSize);
fprintf('Done running %s.m ...\n', mfilename);


figure;
mask = ~data(1).quadMask(1).mask;
imshow(mask);


for vid = 1:nvids
    Img = data(1).frame(:,:,vid); % Image
    
    
    imagesc(Img)
    
%     temp(vid) = mean(trialData(vid).videoData.tempLog(:,2)); %avg temp during video
    
    for roi = 1:4
        
        
        imgAdj = Img; % start with full image
        mask = data(1).quadMask(roi).mask;
        imgAdj(mask) = 0; % crop out frame
        %show masked image
        imagesc(imgAdj)
        
        props = regionprops(~mask, 'BoundingBox');
        maskedImage = imcrop(imgAdj, props.BoundingBox);
        % show cropped image
        imagesc(maskedImage);
        
        
        y(vid,roi) = sum(sum(imgAdj));
    end
end


props = regionprops(mask, 'BoundingBox');
maskedImage = imcrop(Img, props.BoundingBox);
imshow(maskedImage, []);


% Crop the image to the bounding box.
props = regionprops(mask, 'BoundingBox');
maskedImage = imcrop(maskedImage, props.BoundingBox);
% Display it in the lower right plot.
subplot(2, 2, 4);
imshow(maskedImage, []);






