

function image = processOccImg(image, PixelThreshold, ClusterThreshold)
% image = processOccImg(image, PixelThreshold, ClusterThreshold)
% Image processing for each frame
% 
% ES Dickinson, Yale University, Aug 2021

% Image processing:
imGPUoriginal = gpuArray(image);
imGPUgray = rgb2gray(imGPUoriginal); % convert to greyscale
imflyGPU = imGPUgray>PixelThreshold; %color threshold            
image = imopen(imflyGPU,strel('disk',ClusterThreshold)); %denoise

end