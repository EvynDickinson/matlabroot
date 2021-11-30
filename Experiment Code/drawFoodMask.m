
function [X, Y, mask] = drawFoodMask(Img,WellName,X,Y)
% [X, Y, mask] = drawFoodMask(Img,WellName,X,Y)
% eliminate data points from 

    figure;
    imshow(Img)
    title(['Outline the ' WellName ' well'])
    roi = drawpolygon;
    mask = roi.Position;
    
    Xdim = size(X);
    Ydim = size(Y);

    %resize the data:
    X = reshape(X,[numel(X),1]);
    Y = reshape(Y,[numel(Y),1]);
    
    % Find points within the masked region and turn to NaN
    [in,on] = inpolygon(X,Y, mask(:,1),mask(:,2));   % Logical Matrix
    inon = in | on;                                    % Combine ‘in’ And ‘on’
    X(inon) = NaN;
    Y(inon) = NaN;

    % Resize the data to OG structure:
    X = reshape(X,Xdim);
    Y = reshape(Y,Ydim);
end
                                            