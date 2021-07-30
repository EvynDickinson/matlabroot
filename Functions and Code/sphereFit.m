function [Center,Radius] = sphereFit(X)
% this fits a sphere to a collection of data using a closed form for the
% solution (opposed to using an array the size of the data set). 
% Minimizes Sum((x-xc)^2+(y-yc)^2+(z-zc)^2-r^2)^2
% x,y,z are the data, xc,yc,zc are the sphere's center, and r is the radius

% Assumes that points are not in a singular configuration, real numbers, ...
% if you have coplanar data, use a circle fit with svd for determining the
% plane, recommended Circle Fit (Pratt method), by Nikolai Chernov
% http://www.mathworks.com/matlabcentral/fileexchange/22643

% Input:
% X: n x 3 matrix of cartesian data
% Outputs:
% Center: Center of sphere 
% Radius: Radius of sphere
% Author:
% Alan Jennings, University of Dayton

A = [mean(X(:,1).*(X(:,1)-mean(X(:,1)))), ...
    2*mean(X(:,1).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,1).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    mean(X(:,2).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,2).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    0, ...
    mean(X(:,3).*(X(:,3)-mean(X(:,3))))];
A = A+A.';
B = [mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,1)-mean(X(:,1))));...
    mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,2)-mean(X(:,2))));...
    mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,3)-mean(X(:,3))))];
Center = (A\B).';
Radius = sqrt(mean(sum([X(:,1)-Center(1),X(:,2)-Center(2),X(:,3)-Center(3)].^2,2)));

disp(Radius)


% 
% %%
% % Test the 'fit' of the circle: 
% %check the percentage of points close to the euclidian distance from the
% %center: 
% 
% Center  = mean(X);
% r = 1.5; % 'known' radius
% thresh = r * 0.03; %3 percent variation allowed
% min_val = r-thresh;
% max_val = r+thresh;
% n = size(X,1);
% 
% % dist = pdist2(Center,X)'; %euclidian distance from 'Center'
% % onball = sum(dist>min_val & dist<max_val)/n; %percent of points on the sphere
% 
% % for those not close: try manual fit?
% % use the center point as a starting location then iterate through 10,000
% % possibilities between the total?
% 
% % BUILD A CUBE OF TEST POINTS
% % search_range = abs(mean(X)-Center)*2;
% search_range = r*1.5;
% Edge = Center+search_range;
% n = 100; %iterations per axes: 
% x_range = linspace(Center(1)-search_range,Center(1)+search_range,n)';
% y_range = linspace(Center(2)-search_range,Center(2)+search_range,n)';
% z_range = linspace(Center(3)-search_range,Center(3)+search_range,n)';
% % X-Unit:
% xU = repmat(x_range,n^2,1);
% % Y-Unit:
% yU = [];
% MT = ones(n,1);
% for ii = 1:n
%     yU = [yU; y_range(ii)*MT];
% end
% yU = repmat(yU,n,1);
% % Z-Unit:
% zU = [];
% MT = ones(n^2,1);
% for ii = 1:n
%    zU = [zU; z_range(ii)*MT];
% end
% % figure; 
% % scatter3(Center(1), Center(2), Center(3), 150, 'pk', 'filled')
% % hold on
% % scatter3(Center(1), Center(2), Center(3), 50, 'pr', 'filled')
% % % plot sample points
% % scatter3(xU, yU, zU, 20, 'k', 'filled')
% cloud = [xU, yU, zU];
% dist = pdist2(cloud,X); %euclidian distance from 'Center'
% onball = 1-(sum(dist>min_val & dist<max_val,2)/length(X)); %percent of points on the sphere
% alldata = [cloud, onball];
% 
% 
% 
% 
% cmapall = get_color('Green', 'Red', 100);
% C = cmapall(round(onball.*100),:);
% % 
% % figure;
% % scatter3(xU, yU, zU, 20, C, 'filled')
% % hold on
% % scatter3(X(:,1), X(:,2), X(:,3), 50, 'k', 'filled');
% % 
% 
% figure;
% scatter(1:length(xU),onball, 20, C)
% hline(0.1)
% title('Cluster 1')
% 
% % add a cutoff point (10%) for good fit:
% bestof = cloud((onball<0.1),:);
% bestof((bestof(:,3)>0),:)=[];
% 
% figure;
% scatter3(X(:,1), X(:,2), X(:,3), 50, 'k', 'filled');
% hold on
% scatter3(bestof(:,1), bestof(:,2), bestof(:,3), 50, 'g', 'filled');
% title('Cluster 1')
% % new center of ball suggested points:
% 
% 
% %% BUILD A CUBE OF NEW TEST POINTS 
% 
% newpoint = mean(bestof);
% longest_distance = max(pdist2(newpoint,bestof));
% 
% 
% search_range = longest_distance*1.5;
% Edge = newpoint+search_range;
% 
% n = 10; %iterations per axes: 
% x_range = linspace(newpoint(1)-search_range,newpoint(1)+search_range,n)';
% y_range = linspace(newpoint(2)-search_range,newpoint(2)+search_range,n)';
% z_range = linspace(newpoint(3)-search_range,newpoint(3)+search_range,n)';
% % X-Unit:
% xU = repmat(x_range,n^2,1);
% % Y-Unit:
% yU = [];
% MT = ones(n,1);
% for ii = 1:n
%     yU = [yU; y_range(ii)*MT];
% end
% yU = repmat(yU,n,1);
% % Z-Unit:
% zU = [];
% MT = ones(n^2,1);
% for ii = 1:n
%    zU = [zU; z_range(ii)*MT];
% end
% 
% % NEW THRESHOLD VARIABLES
% thresh = r * 0.01; %3 percent variation allowed
% min_val = r-thresh;
% max_val = r+thresh;
% % ------------------
% cloud = [xU, yU, zU];
% dist = pdist2(cloud,X); %euclidian distance from 'Center'
% onball = 1-(sum(dist>min_val & dist<max_val,2)/length(X)); %percent of points on the sphere
% alldata = [cloud, onball];
% 
% 
% figure;
% scatter(1:length(xU),onball, 20, 'k')
% hline(0.1)
% title('cluster 2')
% 
% %  WRK
% % add a cutoff point (10%) for good fit:
% [~,idx] = sort(onball);
% bestof = cloud(idx(1:10),:); % best 10 data points
% finalCenter = median(bestof);
% 
% 
% figure;
% scatter3(X(:,1), X(:,2), X(:,3), 50, 'k', 'filled');
% hold on
% scatter3(bestof(:,1), bestof(:,2), bestof(:,3), 50, 'g', 'filled');
% scatter3(finalCenter(1), finalCenter(2), finalCenter(3), 100, 'pm', 'filled')
% 
% 
% 
% %% Check ball position:
% 
% fig = getfig('',1);
% scatter3(X(:,1), X(:,2), X(:,3));
% hold on
% axis tight
% box off
% set(gca,'visible','off')
% [x,y,z] = sphere;
% s = surf(r*x+finalCenter(1), r*y+finalCenter(2), r*z+finalCenter(3));
% set(s, 'FaceColor', Color('grey'))
% alpha 0.8
% axis vis3d % sets the aspect ratio for 3d rotation







