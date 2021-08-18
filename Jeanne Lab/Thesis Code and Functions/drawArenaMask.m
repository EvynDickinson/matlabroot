
function mask = drawArenaMask(image, numPoly)
% 
% mask = drawArenaMask(image, numPoly)
% 
% ES Dickinson, Yale University, Aug 2021
m = size(image,1);
n = size(image,2);

% ARENA MASK:
f = figure;
imshow(image)
title('Draw circle around the arena')
roi = drawcircle;    
centre = roi.Center;
radius = roi.Radius;
uiwait(f)
% Define coordinates and radius
x1 = centre(1);
y1 = centre(2);
% Generate grid with binary mask representing the circle. 
[xx,yy] = ndgrid((1:m)-y1,(1:n)-x1);
allMasks(:,:,1) = (xx.^2 + yy.^2>radius^2);

% FOOD WELL MASKS
if nargin>1 %only do if option selected
    for ii = 1:numPoly    
        f = figure;
        imshow(image)
        title('Outline the food well')
        roi = drawpolygon;
        BW = poly2mask(roi.Position(:,1),roi.Position(:,2),m,n);
        allMasks(:,:,ii+1) = BW;
        uiwait(f)
    end
end
mask = sum(allMasks,3)>0;

%     imshow(mask)
% 
%     imshowpair(allMasks(:,:,1), BW, 'montage')


end
