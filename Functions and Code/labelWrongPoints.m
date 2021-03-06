
function pts = labelWrongPoints(figHandle,n)

if nargin < 2
    n = Inf;
    pts = zeros(2, 0);
else
    pts = zeros(2, n);
end
% imshow(image);     % display image
xold = 0;
yold = 0;
k = 0;
hold on;           % and keep it there while we plot
while 1
    [xi, yi, but] = ginput(1);      % get a point
    if ~isequal(but, 1)             % stop if not button 1
        break
    end
    k = k + 1;
    pts(1,k) = xi;
    pts(2,k) = yi;
      if xold
          scatter([xold xi], [yold yi], 75, 'r','filled');  % draw as we go
          scatter([xold xi], [yold yi], 40, 'y','filled');  % draw as we go
      else
          scatter(xi, yi, 75, 'r','filled'); % first point on its own
          scatter(xi, yi, 60, 40, 'y','filled'); % first point on its own
      end
      if isequal(k, n)
          break
      end
      xold = xi;
      yold = yi;
end
hold off;
if k < size(pts,2)
    pts = pts(:, 1:k);
end

close(figHandle)
end
