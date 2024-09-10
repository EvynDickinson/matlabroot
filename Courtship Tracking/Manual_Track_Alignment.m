

%% Load fly data

baseFolder = getDataPath(1,2);
trial = '05.07.2024_C1_F_LRR_25-17_n2_A';
load([baseFolder,trial '/' trial(12:end-2) trial(end) ' timecourse data.mat'])

%% Stitch fly tracks together ...

raw_video_dir =  getDataPath(2,2);
vid_date = trial(1:10);


%% 
basePath = getDataPath(1,2);
[file, path] = uigetfile([basePath '*timecourse data.mat']);




% Create a blank image or load your image
figure;
imshow(ones(256, 256));  % Example blank white image

% Create an empty array to store point handles
points = [];

% Add multiple editable (draggable) points
numPoints = 5;  % Number of points you want to add
for i = 1:numPoints
    % Create a draggable point
    points(i) = drawpoint('Position', [randi([50, 200]), randi([50, 200])]);
    
    % Optionally, you can set a callback for when the point is moved
    addlistener(points(i), 'MovingROI', @(src, evt) disp(['Point ' num2str(i) ' moved to: ' num2str(evt.CurrentPosition)]));
end

% Optionally, collect the final positions of all points
finalPositions = cell2mat(arrayfun(@(p) p.Position, points, 'UniformOutput', false)');
disp('Final positions of points:');
disp(finalPositions);
