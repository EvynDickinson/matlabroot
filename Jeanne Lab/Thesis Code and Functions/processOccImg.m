

function image = processOccImg(image, GPU_yn, PixelThreshold, ClusterThreshold, ~)
% image = processOccImg(image, GPU_yn, PixelThreshold, ClusterThreshold, binaryoption)
% Image processing for each frame on GPU or CPU
% 
% ES Dickinson, Yale University, Aug 2021


% default image processing values:
P_thresh = 60;
Clust = 5;

% Image processing:
switch nargin
    case 1 
        GPU_yn = true; %default to GPU image processing
        ClusterThreshold = Clust;
        PixelThreshold = P_thresh;
    case 2
        ClusterThreshold = Clust;
        PixelThreshold = P_thresh;
    case 3
        ClusterThreshold = Clust;
end



% where is the image going to be processed
if GPU_yn
    IMG = gpuArray(image);
else
    IMG = image;
end

if nargin == 5 % return the sharpened but no binarized image
    IMG = image;
    I = rgb2gray(IMG); % convert to greyscale
    I = imadjust(I);
    image = imsharpen(I,'Radius',20,'Amount',3);
else % return the binarized image
    I = rgb2gray(IMG); % convert to greyscale
    I = imadjust(I);
    I = imsharpen(I,'Radius',20,'Amount',3);
    I = I>PixelThreshold; %color threshold
    image = imopen(I,strel('disk',ClusterThreshold)); %denoise
end


end