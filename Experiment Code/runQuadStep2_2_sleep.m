
function results = runQuadStep2_2_sleep(trialPath,arena, expName)

% check if data already exists before remaking it: 
sleep_file = [trialPath expName ' sleeping data v2.mat'];  
if isfile(sleep_file)
    disp(['Skipping ' trialPath expName ' sleep data v2 exists'])
    results = true;
    return
end

% load data processed from quadstep2_2
load([trialPath expName arena ' timecourse data v2.mat'])

% clearvars('-except',initial_vars{:})
nbins = 50;

% Create sleep data for unprocessed files (trial by trial)

try fps = expData.parameters.FPS;
catch 
    disp('No FPS found')
end

% How long does a fly need to be still to count as 'sleep'
min_duration = 5*fps*60; % 5 mins * 3fps*60sec = data point number that must be met 'unmoving'

sleepData = struct;
%preallocate for speed and space
trial_length = size(data.T,1);

% Set axis limits for the selected arena
x = data.centre(1);
y = data.centre(2);
r = data.r;
xlimit = [x-(r+50),x+(r+50)];
ylimit = [y-(r+50),y+50+r];

% find the 'auto bin' lines
xedge = linspace(xlimit(1),xlimit(2),nbins+1);
yedge = linspace(ylimit(1),ylimit(2),nbins+1);

% pull the fly locations during the trial
x_loc = data.x_loc;
y_loc = data.y_loc;

% Testing faster method: 
% Get all frame indices for each fly
[frameIdx, flyIdx] = find(~isnan(x_loc));
X = x_loc(sub2ind(size(x_loc), frameIdx, flyIdx));
Y = y_loc(sub2ind(size(y_loc), frameIdx, flyIdx));

% Bin each point
[~,~,xBin] = histcounts(X, xedge);
[~,~,yBin] = histcounts(Y, yedge);

% Remove out-of-bound bins
valid = xBin > 0 & yBin > 0;
xBin = xBin(valid); 
yBin = yBin(valid); 
frameIdx = frameIdx(valid);

% Use accumarray to build N
binLinearIdx = sub2ind([nbins nbins trial_length], xBin, yBin, frameIdx);
N = accumarray(binLinearIdx, 1, [nbins*nbins*trial_length 1]);
N = reshape(N, [nbins, nbins, trial_length]);

% find grid space that have continuous occupation for more than min_duration
frameCount = (nan(nbins,nbins,trial_length));
frameCount(:,:,1) = N(:,:,1);
for frame = 2:trial_length
    currFrame = N(:,:,frame); % current frame locations
    resetLoc = currFrame==0; % locations that do not have flies and thus need a count reset

    tempCount = frameCount(:,:,frame-1)+currFrame; % add current frames to list
    tempCount(resetLoc) = 0; % reset counts for spots with no flies

    frameCount(:,:,frame) = tempCount; % add current count into the saving structure
end

% ---- Vectorize the data (find the flies that are sleeping....) -----
% create empty matrixes for the x and y positions of the sleeping flies

% add a catch here if there is not data for the number of flies??
if isnan(data.nflies) || data.nflies<=0
    [excelfile, Excel, ~] = load_QuadBowlExperiments;
    xlRow = find(strcmp(expName,excelfile(:,Excel.expID)) & ...
        strcmp(expData.parameters.date,excelfile(:,Excel.date)) & ...
        strcmp(arena,excelfile(:,Excel.arena)));
    data.nflies = excelfile{xlRow,Excel.numflies};
end

sleeping = struct;
[sleeping.X, sleeping.Y, sleeping.all_distance] = deal(nan(trial_length,data.nflies));
sleeping.sleepNum = zeros(trial_length,1);
[sleeping.dist_avg, sleeping.dist_err] = deal(nan(trial_length,1));

c1 = data.wellcenters(1,expData.parameters.foodLoc);
c2 = data.wellcenters(2,expData.parameters.foodLoc);

% assign data by frame
for frame = 1:trial_length
    frame_data = frameCount(:,:,frame) > min_duration;
    binLoc = find(frame_data>0);
    
    % Find the coordinates of the sleeping flies bins from the discretized data
    y_row = ceil(binLoc/nbins);
    x_row = rem(binLoc-1,nbins)+1;
    x_position = (xedge(x_row) + xedge(x_row+1))/2;
    y_position = (yedge(y_row) + yedge(y_row+1))/2;
    
    % add position data to the matrix:
    if ~isempty(binLoc)
        % number of flies sleeping
        sleepNum = length(x_position);
        sleeping.sleepNum(frame) = sleepNum;
        % location of sleeping flies
        sleeping.X(frame,1:sleepNum) = x_position;
        sleeping.Y(frame,1:sleepNum) = y_position;
        % distance to food...
        dX = (x_position-c1);
        dY = (y_position-c2);
        temp_dist = hypot(dX,dY)./data.pix2mm;
        sleeping.all_distance(frame,1:sleepNum) = temp_dist;
        % average distance:
        sleeping.dist_avg(frame) = mean(temp_dist);
        sleeping.dist_err(frame) = std(temp_dist);
    end
end

save(sleep_file,'sleeping','-v7.3'); 
results = 'saved file';
