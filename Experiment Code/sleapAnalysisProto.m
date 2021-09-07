
% https://sleap.ai/notebooks/Analysis_examples.html <-- python code for data example

filepath = 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\Tracking\labels.v000.analysis.h5';
occupancy_matrix = h5read(filepath,'/track_occupancy');
tracks_matrix = h5read(filepath,'/tracks');


h5disp(filepath)


%% Basic sorting of the data...
% demo only -- load in the camera image for the frame and overlay...
demovid = 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\NonlinearCooling_1.avi';
movieInfo = VideoReader(demovid);
% mat locations: 
ML.frame = 1;   % frames are 1st dim
ML.node = 2;    % 
ML.xy = 3;      % xy coordinates
ML.fly = 4;     % tracks for 'a' fly

% all data for head tracked location
headdata = squeeze(tracks_matrix(:,1,:,:));

frame = 5;

% tracked points
x = squeeze(headdata(frame, 1, :));
y = squeeze(headdata(frame, 2, :));
x(isnan(x)) = []; % remove empty tracks
y(isnan(y)) = [];
% image
img = read(movieInfo,frame);

% overlay:
fig = getfig; set(fig, 'color', 'k')
hold on
imagesc(img)
scatter(x,y, 30, 'y')
axis tight; axis square
set(gca, 'visible', 'off')

% save_figure(fig, 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\Tracking\demo tracking image','-png');


%% 
pts = readPoints(img,4); % get locations of the wells from the image

x = squeeze(headdata(:, 1, :)); % frame, instance
y = squeeze(headdata(:, 2, :));

X = reshape(x,numel(x),1);
Y = reshape(y,numel(y),1);
loc = isnan(X);
X(loc) = [];
Y(loc) = [];
% squash all points together over the video:


% what if we do a min radius from the coordinates of the food bowl?  
fig = getfig; set(fig, 'color', 'k');
hist2d(X,Y, 'probability', 'tile')
axis tight; axis square
set(gca, 'visible', 'off')
c = colorbar;
c.Color = [1,1,1];
hold on
scatter(pts(1,:),pts(2,:), 75, 'r', 'filled') % could automate the color on this

% save_figure(fig, 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\Tracking\demo histogram vid 1','-png');

%% Pull temperature information for each frame

% get temp log
tempLogName = 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\NonlinearCooling_RampLog.csv';
tempLog = readmatrix(tempLogName);

% what are the temp time points for this video start and stop: 
% -- TODO load in the appropriate matlab file for this video...harvest code
% from VidAnalysisStep1...


% https://sleap.ai/notebooks/Training_and_inference_using_Google_Drive.html

%% Find flies within ROI of each well

% how many points within a specific radius of the well?
radii = 200; %distance must be less than this number to count for a well ROI
nInst = length(X);
% loop for all wells
for well = 1:4
    b = pts(:,well); % well 1 points
    euDist = sqrt((b(1)-X).^2 + (b(2)-Y).^2);
    data(well).well_loc = (euDist<=radii);
    data(well).well_count = sum((euDist<=radii));
    
    % occupation probability
    data(well).occ_prob = data(well).well_count/nInst;
    data(well).nInst = nInst;
end







% % visual confirmation that the selected points are near the well:
% fig = getfig; set(fig, 'color', 'k');
% hist2d(X,Y, 'probability', 'tile')
% axis tight; axis square
% set(gca, 'visible', 'off')
% c = colorbar;
% c.Color = [1,1,1];
% hold on
% scatter(X(well_loc),Y(well_loc), 50, 'y') % could automate the color on this
% scatter(pts(1,:),pts(2,:), 75, 'r', 'filled') % could automate the color on this
% viscircles(b',radii)
% 
% save_figure(fig, 'G:\My Drive\Jeanne Lab\DATA\08.27.2021\Tracking\demo histogram point sel','-png');
% 

%% 











