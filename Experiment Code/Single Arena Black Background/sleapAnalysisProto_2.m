% sleapAnalysisProto2

%% Select Date & Experiment to Post-process
clear
close all
clc  
%get base folder pathway
[baseFolder, folder] = getCloudPath(2); 

% Select the complete experiments to process
list_dirs = dir([baseFolder folder, '\*.mat']); %only matlab files
list_dirs = {list_dirs(:).name};
expNames = cellfun(@(x) x(1:end-11),list_dirs,'UniformOutput',false); %pull root name
expName = expNames{listdlg('ListString', expNames, 'SelectionMode', 'Single')};
% Pull fly summary sheet information on selected experiment
[excelfile, Excel, xlFile] = load_FlyBowlExperiments;
XLrow = find(strcmpi(excelfile(:,Excel.date), folder) & ...
           strcmpi(excelfile(:,Excel.expID), expName)); %rows with sel date

% Get analysis file directory
analysisDir = [baseFolder folder '\analysis\'];

% Load relevant data files (.mat, .csv, .h5)
warning off
% load matlab file for experiment
expData = load([baseFolder folder '\' expName 'dataMat.mat']);

% load the temperature log for the experiment
tempLog = readmatrix([baseFolder folder '\' expName '_RampLog']);
nvids = expData.num.vids;

% load number of flies
nflies = excelfile{XLrow,Excel.numflies};

% load pre-processed data
load([analysisDir expName ' timecourse data.mat'])

%re-select well centers if the data isn't preloaded:
if ~exist('wellcenters')
    fprintf('\nSelect wells in this order:\n'); display(wellLabels)
    movieInfo = VideoReader([baseFolder folder '\' expName '_1.avi']); %read in video
    demoImg = rgb2gray(read(movieInfo,1));
    pts = readPoints(demoImg,4); % get locations of the wells from the image
    wellcenters = pts;
end

% add total well occupancy to occupancy structure
occupancy.allwellOcc = sum(occupancy.occ,2);

%% Calculate degree of clustering
% pull info for each video

LDist = [];
for vid = 1:nvids
    flyDistance = [];
    % all data for head tracked location
    headData = squeeze(data(vid).tracks(:,1,:,:));
    nframes = size(headData,1); 
    
    % x-y coordinates of flies for each frame
    x_loc = squeeze(headData(:,1,:));
    y_loc = squeeze(headData(:,2,:));
    
    for frame = 1:nframes
        test = [x_loc(frame,:)',y_loc(frame,:)'];
        D = pdist(test);
        D(isnan(D)) = [];
        flyDistance(frame) = mean(D);
    end
    data(vid).flyDistance = flyDistance;
    LDist = [LDist,flyDistance]; 
end
time = linspace(1,(length(LDist)/3)/60,length(LDist));
occupancy.time = time;
%add inter-fly-distance to the occupancy time course structure
occupancy.IFD = LDist;

% demo the clustering proxy
if nvids>12
    divisor = round(nvids/6);
    vidList = 1:divisor:nvids;
else
    vidList = 1:nvids;
end
nrow = 2; ncol = length(vidList); ii = 0; 
[minidx,minidx] = deal([]);
fig = getfig('',1); set(fig, 'pos',[120 331 1244 368], 'color', 'k');
for vid = vidList
    ii = ii+1;
    movieInfo = VideoReader([baseFolder folder '\' expName '_' num2str(vid) '.avi']); %read in video
    headData = squeeze(data(vid).tracks(:,1,:,:));
    
    % Most clustered:
    [M,minidx(ii)] = min(data(vid).flyDistance);
    img = read(movieInfo,minidx(ii));
    subplot(nrow, ncol, ii)
    imshow(img); hold on
    axis tight square
    set(gca, 'visible', 'off')
    x = squeeze(headData(minidx(ii), 1, :)); x(isnan(x)) = [];
    y = squeeze(headData(minidx(ii), 2, :)); y(isnan(y)) = [];
    scatter(x,y, 10, 'y', 'filled')
    title(num2str(M))
    % overlay 'size bar' for min dist:
    plot([10,10+M], [20,20], 'linewidth', 0.5, 'color', 'r')
    
    % Least clustered:
    [M,maxidx(ii)] = max(data(vid).flyDistance);
    img = read(movieInfo,maxidx(ii));
    subplot(nrow, ncol, ii+length(vidList))
    imshow(img); hold on
    axis tight square
    set(gca, 'visible', 'off')
    x = squeeze(headData(maxidx(ii), 1, :)); x(isnan(x)) = [];
    y = squeeze(headData(maxidx(ii), 2, :)); y(isnan(y)) = [];
    scatter(x,y, 10, 'y', 'filled')
    title(num2str(M))
    % overlay 'size bar' for min dist:
    plot([10,10+M], [20,20], 'linewidth', 0.5, 'color', 'r')
end
labelHandles = findall(gcf, 'type', 'text', 'handlevisibility', 'off');
set(labelHandles,'FontSize', 15, 'color', 'w')

save_figure(fig, [analysisDir expName ' linear clustering demo'], '-png');

clear x y M minidx maxidx img nrow ncol D test flyDistance x_loc y_loc labelHandles headData

%% Movement analysis

%bin frame and check bin occupation changes across frames as proxy for
%movement --> won't give an accurate 'speed' etc. but it might give
%movement.
tic
nbins = 100;
spaceChange = [];

for vid = 1:nvids   
    N = [];
    % all data for head tracked location
    headData = squeeze(data(vid).tracks(:,1,:,:));
    nframes = size(headData,1); 
    
    % x-y coordinates of flies for each frame
    x_loc = squeeze(headData(:,1,:));
    y_loc = squeeze(headData(:,2,:));

    for frame = 1:nframes
        X = x_loc(frame,:); X(isnan(X)) = [];
        Y = y_loc(frame,:); Y(isnan(Y)) = [];
        N(:,:,frame) = histcounts2(X,Y,nbins);
    end
    if vid>1 % add last frame of previous vid to prevent frame loss 
        N = cat(3, lastFrame, N);
    end
    lastFrame =  N(:,:,end);
    binDiff = diff(N,1,3);
    binDiff(binDiff<0) = 0;
    BinChange = squeeze(sum(sum(binDiff,1),2));
    % save pertinent information:
    data(vid).BinChange = BinChange;
    spaceChange = [spaceChange; BinChange];
end
occupancy.movement = spaceChange;
toc

figure
plot(occupancy.time(2:end), occupancy.movement)
hold on
plot(occupancy.time(2:end), moving_average(occupancy.movement,10))



%% Time course feature comparison
nrow = 8; ncol = 1;
sSpan = 180; %1 minute filter length
sbpts(1).idx = 1;
sbpts(2).idx = 2:4;
sbpts(3).idx = 5:6;
sbpts(4).idx = 7:8;
LW = 1.5;



fig = getfig; set(fig, 'pos', [157 86 1232 878])
subplot(nrow,ncol,sbpts(1).idx)
plot(time,moving_average(occupancy.temp,sSpan),'linewidth', LW, 'color', 'w')
ylabel('(\circ)')
title('temperature')

subplot(nrow,ncol,sbpts(2).idx)
hold on
for well = 1:4
    switch expData.params.(['well_' num2str(well)'])
        case 'Yeast'
            kolor = Color('gold');
        case 'Plant'
            kolor = Color('green');
        case 'Empty'
            kolor = Color('grey');
        case 'Plant_827'
            kolor = Color('palegreen');
        case 'Plant_91'
            kolor = Color('Darkgreen');
    end
    plot(time,moving_average(occupancy.occ(:,well),sSpan),'linewidth', LW, 'color', kolor)
end
plot(time,moving_average(occupancy.allwellOcc,sSpan),':','linewidth',LW,'color', Color('slateblue'))
ylabel('occupancy probability')
title('individual well occupancy')
legend(wellLabels)
l = legend(strrep([wellLabels; {'all wells'}],'_','-'));
set(l, 'color', 'k', 'textcolor', 'w','edgecolor', 'k','position',...
    [0.1552 0.6918 0.0963 0.1126]);% [0.8780 0.8119 0.0963 0.1126])%

subplot(nrow,ncol,sbpts(3).idx)
plot(time,moving_average(occupancy.IFD,sSpan),'linewidth', LW, 'color', 'w')
ylabel('pixels')
title('inter-fly-distance')

subplot(nrow,ncol,sbpts(4).idx)
plot(time(2:end),moving_average(occupancy.movement,sSpan),'linewidth', LW, 'color', 'w')
ylabel('(a.u.)')
title('movement')
xlabel('time (min)')

formatFig(fig,true,[nrow, ncol], sbpts);
for ii = 1:3
    subplot(nrow,ncol,sbpts(ii).idx)
    set(gca, 'XColor', 'k')
end

save_figure(fig, [analysisDir expName ' timecourse features'], '-png');

clear sbpts well nrow ncol LW ans ii

%% 


















