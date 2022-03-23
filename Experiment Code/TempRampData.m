%% 
% Find the video number that corresponds to a particular temperature range
% for the experiment for visual comparision

%%
% Select data to load:
[baseFolder, folder] = getCloudPath(2);   
 
list_dirs = dir([baseFolder folder, '/*dataMat.mat']); %only matlab files
list_dirs = {list_dirs(:).name};
expNames = cellfun(@(x) x(1:end-11),list_dirs,'UniformOutput',false); %pull root name
expName = expNames{listdlg('ListString', expNames, 'SelectionMode', 'Single')};
expName = expName(1:end-1);clear expNames

% load data
expData = load([baseFolder folder '/' expName ' dataMat.mat']);
tempLog = readmatrix([baseFolder folder '/' expName '_RampLog']);
nvids = expData.parameters.numVids;

%% Extrapolate temperature data

[vidNums, vidFrame, temperature, tempWork] = deal([]);
for vid = 1:nvids
    movieInfo = VideoReader([baseFolder folder '/' expName '_1.avi']); %read in video
    nframes = movieInfo.NumFrames;
    % video numbers
    vidNums = [vidNums; vid*ones(nframes,1)];
    % frame num in video
    vidFrame = [vidFrame; (1:nframes)'];
    % temperature log
    logROI(1) = find(tempLog(:,1)==expData.tempLogStart(vid,3));
    logROI(2) = find(tempLog(:,1)==expData.tempLogEnd(vid,3));
    tempCourse = tempLog(logROI(1):logROI(2),2);
    x = round(linspace(1, nframes, length(tempCourse)));
    fullTempList = interp1(x,tempCourse,1:nframes,'spline');   
    data(vid).tempLog = fullTempList;
    temperature = [temperature; fullTempList']; 
    % temp plate work log
    workCourse = tempLog(logROI(1):logROI(2),4);
    x = round(linspace(1, nframes, length(workCourse)));
    fullWorkList = interp1(x,workCourse,1:nframes,'spline');   
    tempWork = [tempWork; fullWorkList']; 
end

% Time count
time = (linspace(1, (length(temperature)/3)/60, length(temperature)))';

row = 4; col = 1;
figure;
subplot(row,col,1)
plot(time, vidNums)
subplot(row,col,2:4)
plot(time, temperature)

hline(14)
hline(20)
search_temp = 14;

% what videos have which temperatures?? 
idx = find(temperature==search_temp);
vline(time(idx))

% get video & frame numbers for these points:
frameNumbers = vidFrame(idx);
videoNumbers = vidNums(idx);

T = table(videoNumbers, frameNumbers);
disp(T)














