
function [center, fig, err] = DLC_MANUALsphereFit(X, r)
% [center, fig, err] = DLC_MANUALsphereFit(X, r)
% Input:
% X -> nx3 matrix of points on the predicted sphere
% r -> predicted/known radius size
% Output: 
% Center -> cartesian coordinates for the center of the sphere
% err -> percentage of input points that fall within 3% of surface
% 
% ES Dickinson & S Walling-Bell
% University of Washington, 2020

fig = getfig('Sphere fit',1);

%%
% Test the 'fit' of the circle: 
%check the percentage of points close to the euclidian distance from the
%center: 
% r = 1.5; % 'known' radius

Center  = mean(X);
thresh = r * 0.03; %3 percent variation allowed
min_val = r-thresh;
max_val = r+thresh;
n = size(X,1); %number of data points given to function

% ---- BUILD A CUBE OF TEST POINTS ---
search_range = r*1.5; % range around avg center point to search
Edge = Center+search_range; %extrenes to test
n = 100; %iterations per axes: total points = n^3
x_range = linspace(Center(1)-search_range,Center(1)+search_range,n)';
y_range = linspace(Center(2)-search_range,Center(2)+search_range,n)';
z_range = linspace(Center(3)-search_range,Center(3)+search_range,n)';
% X-Unit:
xU = repmat(x_range,n^2,1);
% Y-Unit:
yU = [];
MT = ones(n,1);
for ii = 1:n
    yU = [yU; y_range(ii)*MT];
end
yU = repmat(yU,n,1);
% Z-Unit:
zU = [];
MT = ones(n^2,1);
for ii = 1:n
   zU = [zU; z_range(ii)*MT];
end
% figure; 
% scatter3(Center(1), Center(2), Center(3), 150, 'pk', 'filled')
% hold on
% scatter3(Center(1), Center(2), Center(3), 50, 'pr', 'filled')
% % plot sample points
% scatter3(xU, yU, zU, 20, 'k', 'filled')
cloud = [xU, yU, zU];
dist = pdist2(cloud,X); %euclidian distance from 'Center'
onball = 1-(sum(dist>min_val & dist<max_val,2)/length(X)); %percent of points on the sphere
% 
% a = (dist-r).^2;
% b = mean(a,2);
% c = sqrt(b);
% % RMSE of the thousands of test points...
% test = sqrt(mean((dist-r).^2,2));
% % loc = 
% % figure; 
% % plot(1:n^3,test)



% get colors for plotting axes:
% cmapall = get_color('Green', 'Red', 100);
% C = cmapall(round(onball.*100),:);
C = 'k';

% 
% figure;
% scatter3(xU, yU, zU, 20, C, 'filled')
% hold on
% scatter3(X(:,1), X(:,2), X(:,3), 50, 'k', 'filled');
% 

subplot(2,2,1)
scatter(1:length(xU),onball, 20, C)
hline(0.1) 
xlabel('center test point')
ylabel('bad points (%)')
title('Cluster 1 with cutoff')

% add a cutoff point (10%) for good fit:
loc = (onball<0.1);
% if loc<2
%     loc = (onball<0.3);
%     warndlg('No versions with 90% cutoff. Using 70%')
% end
bestof = cloud(loc,:);
bestof((bestof(:,3)>0),:)=[];
if sum(sum(bestof))==0
    warndlg('NO DATA FIT')
    return
end

subplot(2,2,3)
scatter3(X(:,1), X(:,2), X(:,3), 50, 'k', 'filled');
hold on
scatter3(bestof(:,1), bestof(:,2), bestof(:,3), 50, 'g', 'filled');
title('Cluster locations')
xlabel('x')
ylabel('y')
zlabel('z')


%% Second iteration: FIND THE BEST OF THE BEST:
newpoint = mean(bestof);
longest_distance = max(pdist2(newpoint,bestof));

search_range = longest_distance*1.5;
Edge = newpoint+search_range;

n = 10; %iterations per axes: smaller number than previously
x_range = linspace(newpoint(1)-search_range,newpoint(1)+search_range,n)';
y_range = linspace(newpoint(2)-search_range,newpoint(2)+search_range,n)';
z_range = linspace(newpoint(3)-search_range,newpoint(3)+search_range,n)';
% X-Unit:
xU = repmat(x_range,n^2,1);
% Y-Unit:
yU = [];
MT = ones(n,1);
for ii = 1:n
    yU = [yU; y_range(ii)*MT];
end
yU = repmat(yU,n,1);
% Z-Unit:
zU = [];
MT = ones(n^2,1);
for ii = 1:n
   zU = [zU; z_range(ii)*MT];
end

% NEW THRESHOLD VARIABLES
thresh = r * 0.03; %1 percent variation allowed
min_val2 = r-thresh;
max_val2 = r+thresh;
% ------------------
cloud = [xU, yU, zU];
dist = pdist2(cloud,X); %euclidian distance from 'Center'
onball = 1-(sum(dist>min_val2 & dist<max_val2,2)/length(X)); %percent of points on the sphere

subplot(2,2,2)
scatter(1:length(xU),onball, 20, 'k')
title('cluster 2')
xlabel('test point')
ylabel('bad points (%)')

% % add a cutoff point (10#) for good fit:
% [~,idx] = sort(onball);
% bestof = cloud(idx(1:10),:); % best 10 data points
% finalCenter = median(bestof);
% add a cutoff point (10#) for good fit:
[~,idx] = sort(onball);
bestof = cloud(idx(1),:); % best 10 data points
finalCenter = (bestof);



subplot(2,2,3)
hold on
% scatter3(bestof(:,1), bestof(:,2), bestof(:,3), 50, 'b', 'filled');
scatter3(finalCenter(1), finalCenter(2), finalCenter(3), 100, 'pm', 'filled')



% Find the number of points that were within the threshold for the final
% center point selected:
dist = pdist2(finalCenter,X); %euclidian distance from 'Center'
err = (sum(dist>min_val & dist<max_val,2)/length(X))*100; %percent of points on the sphere

center = finalCenter;

% %% Check ball position:
% 
% fig = getfig('',1);
subplot(2,2,4)
scatter3(X(:,1), X(:,2), X(:,3));
hold on
axis tight
box off
set(gca,'visible','off')
[x,y,z] = sphere;
s = surf(r*x+center(1), r*y+center(2), r*z+center(3));
set(s, 'FaceColor', Color('grey'))
alpha 0.8
axis vis3d % sets the aspect ratio for 3d rotation



end